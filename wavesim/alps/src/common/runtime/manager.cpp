//
// Created by xuanx004 on 3/2/24.
//

#include "manager.h"

#include <common/async/event.h>
#include <common/async/streams.h>
#include <common/async/thread_pool.h>
#include <common/base/logging.h>
#include <common/device/device_traits.h>
#include <common/device/devices.h>
#include <common/utils/timeout_mpi_barrier.h>

#include <Kokkos_Core.hpp>
#include <mpipp/comm.h>
#include <mpipp/environment.h>

#include <stdexcept>

namespace alps {

class RuntimeManager::Impl
{
 public:
  Impl() = default;

  void reinit_stream_pool(int pool_size = default_stream_pool_size)
  {
    if (pool_size <= 0) {
      throw std::invalid_argument("Stream pool size must be positive");
    }
    stream_pool_.reset();
    stream_pool_.emplace(pool_size);
  }

  void reinit_host_transfer_stream_pool()
  {
    host_transfer_stream_pool_.reset();
    host_transfer_stream_pool_.emplace(2);
  }

  void reinit_event_pool()
  {
    event_pool_.reset();
    event_pool_.emplace();
  }

  void reset_thread_pool(unsigned int          n_threads,
                         std::function<void()> init_func)
  {
    thread_pool_.reset();
    //  call num_threads() to circumvent a nodiscard warning
    if (auto n_th =
          thread_pool_
            .emplace(async::concurrency_t(n_threads), std::move(init_func))
            .num_threads();
        n_th.value() != n_threads) {
      throw std::runtime_error("Failed to initialize thread pool with the "
                               "requested number of threads");
    }
  }

  static constexpr int default_stream_pool_size = 16;

  static constexpr unsigned default_thread_pool_size = 2u;

  std::optional<StreamPool<Kokkos::DefaultExecutionSpace>> stream_pool_;

  std::optional<StreamPool<Kokkos::DefaultExecutionSpace>>
    host_transfer_stream_pool_;

  std::optional<async::thread_pool> thread_pool_;

  std::optional<async::event_pool<Kokkos::DefaultExecutionSpace>> event_pool_;

  std::vector<std::function<void()>> cleanup_functions_;

  ~Impl()
  {
    // Execute cleanup functions in reverse order
    while (!cleanup_functions_.empty()) {
      cleanup_functions_.back()();
      cleanup_functions_.pop_back();
    }

    event_pool_.reset();

    stream_pool_.reset();
    host_transfer_stream_pool_.reset();

    thread_pool_.reset();
  }
};

RuntimeManager::ScopeGuard::ScopeGuard()
  : manager_{RuntimeManager::instance()}
{}

RuntimeManager::ScopeGuard::~ScopeGuard()
{
  manager_.finalize();
}

int RuntimeManager::execute_main(std::function<int(int, char**)> const& func,
                                 int                                    argc,
                                 char**                                 argv)
{
  auto guard = RuntimeManager::ScopeGuard{};
  try {
    return func(argc, argv);
  } catch (std::exception const& e) {
    default_logger()->error("Exception caught: {}", e.what());
    return 1;
  } catch (...) {
    default_logger()->error("Unknown exception caught");
    return 1;
  }
}

void RuntimeManager::init_runtimes()
{
  init_mpi();
  init_kokkos();
}

void RuntimeManager::init_mpi()
{
  mpipp::init(mpipp::threading_modes::multiple);
  default_logger()->set_pattern(get_logger_pattern_mpi());
}

void RuntimeManager::push_cleanup_function(
  std::function<void()> cleanup_function)
{
  get_impl().cleanup_functions_.push_back(std::move(cleanup_function));
}

StreamPool<Kokkos::DefaultExecutionSpace> const& RuntimeManager::stream_pool()
{
  if (!get_impl().stream_pool_.has_value()) {
    get_impl().reinit_stream_pool();
  }
  return get_impl().stream_pool_.value();
}

Kokkos::DefaultExecutionSpace RuntimeManager::device_to_host_stream()
{
  if (!get_impl().host_transfer_stream_pool_.has_value()) {
    get_impl().reinit_host_transfer_stream_pool();
  }
  return get_impl().host_transfer_stream_pool_.value().get_stream(1);
}

Kokkos::DefaultExecutionSpace RuntimeManager::host_to_device_stream()
{
  if (!get_impl().host_transfer_stream_pool_.has_value()) {
    get_impl().reinit_host_transfer_stream_pool();
  }
  return get_impl().host_transfer_stream_pool_.value().get_stream(2);
}

async::event_pool<Kokkos::DefaultExecutionSpace>& RuntimeManager::event_pool()
{
  if (!get_impl().event_pool_.has_value()) {
    get_impl().reinit_event_pool();
  }
  return get_impl().event_pool_.value();
}

memory::PooledMemoryResources const& RuntimeManager::pooled_memory_resources()
{
  if (!Kokkos::is_initialized()) {
    // memory pools are not dependent on Kokkos being initialized, but Kokkos
    // initialization ensures the correct device being selected
    throw std::runtime_error(
      "Device runtime (Kokkos) must be initialized before using pooled "
      "memory resources");
  }
  if (!pooled_mr_.has_value()) {
    pooled_mr_.emplace();
  }
  return pooled_mr_.value();
}

async::thread_pool& RuntimeManager::thread_pool()
{
  if (get_impl().thread_pool_.has_value()) {
    return get_impl().thread_pool_.value();
  }

  auto const id = [] {
    if (Kokkos::is_initialized()) {
      return Kokkos::device_id();
    }
    return -1;
  }();

  if (id < 0 && has_device_v) {
    throw std::runtime_error(
      "Thread pool cannot be initialized without a device ID (likely because "
      "Kokkos is not initialized)");
  }
  get_impl().reset_thread_pool(Impl::default_thread_pool_size,
                               [=] { set_device(id); });
  return get_impl().thread_pool_.value();
}

RuntimeManager::RuntimeManager()
  : instance_{[] {
    // initialize the logger so the logger outlives the runtime manager
    default_logger()->flush(); // a "no-op" to initialize the default logger
    return nullptr;
  }()}
{}

RuntimeManager::~RuntimeManager() = default;

RuntimeManager::Impl& RuntimeManager::get_impl()
{
  if (instance_ != nullptr) return *instance_;

  if (!Kokkos::is_initialized() || Kokkos::is_finalized()) {
    throw std::runtime_error(
      "Device runtime (Kokkos) must be initialized before using the runtime "
      "resources");
  }
  instance_ = std::make_unique<Impl>();
  return *instance_;
}

void RuntimeManager::finalize()
{
  instance_.reset();
  if (pooled_mr_.has_value()) {
    pooled_mr_->release_all();
  }

  // Finalize the distributed environment (MPI) before local env (Kokkos)
  if (mpipp::initialized() && !mpipp::finalized()) {
    // wait for all ranks to reach the barrier
    default_logger()->debug("Finalizing MPI...");
    if (auto success = utils::timeout_mpi_barrier(mpipp::COMM_WORLD(), 15);
        !success) {
      default_logger()->error("Timeout in MPI barrier during finalization.");
    } else {
      mpipp::finalize();
    }
  }

  if (Kokkos::is_initialized() && !Kokkos::is_finalized()) {
    Kokkos::finalize();
  }

  pooled_mr_.reset(); // memory resources are destroyed after Kokkos in case
                      // there are still Kokkos objects using them
}

namespace detail {
/// Get the MPI rank with respect to the local node.
std::pair<int, int> get_mpi_node_rank_and_size()
{
  if (!mpipp::initialized() || mpipp::finalized()) {
    throw std::runtime_error("MPI is not available.");
  }

  auto shared_comm = mpipp::COMM_WORLD().split_shared();
  auto rank        = shared_comm.rank();
  auto size        = shared_comm.size();
  default_logger()->debug(
    "MPI local rank: {}/{} on {}", rank, size, mpipp::processor_name());
  return {rank, size};
}

/// @brief Get the device ID for the current MPI rank.
/// The device is assigned round-robin according to the on-node MPI rank.
int get_gpu_roundrobin()
{
  auto [local_rank, node_size] = get_mpi_node_rank_and_size();
  auto local_num_devices       = get_device_count();

  if (local_num_devices < 1) {
    throw std::runtime_error("No GPU device found.");
  }
  if (local_num_devices < node_size) {
    default_logger()->warn(
      "GPU oversubscription: {0} GPUs for {1} MPI ranks on {2}",
      local_num_devices,
      node_size,
      mpipp::processor_name());
  }

  // assign round-robin using on-node MPI rank
  return local_rank % local_num_devices;
}

auto set_kokkos_arguments()
{
  Kokkos::InitializationSettings kokkos_arguments;

  if (alps::has_device_v) {
    int gpu_id = get_gpu_roundrobin();
    kokkos_arguments.set_device_id(gpu_id);
  }

  return kokkos_arguments;
}
} // namespace detail

void RuntimeManager::init_kokkos()
{
  if (!mpipp::initialized()) {
    throw std::runtime_error("MPI is not initialized.");
  }

  Kokkos::initialize(detail::set_kokkos_arguments());

  auto world = mpipp::COMM_WORLD();
  world.barrier();
  if (world.rank() == 0) {
    if (default_logger()->should_log(spdlog::level::debug)) {
      std::stringstream sbuf;
      Kokkos::print_configuration(sbuf, true);
      default_logger()->debug(sbuf.str());
    }
  }
  if (has_device_v) {
    default_logger()->info("Selects GPU device {1} on {0}",
                           mpipp::processor_name(),
                           Kokkos::device_id());
    world.barrier();
    if (world.rank() == 0) {
      if constexpr (alps::mpi_is_gpu_aware_v) {
        default_logger()->info("Use GPU-aware MPI communication mode.");
      } else {
        default_logger()->info("Use host-transfer MPI communication mode.");
      }
    }
  }
}
} // namespace alps
