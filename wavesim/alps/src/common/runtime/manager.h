//
// Created by xuanx004 on 3/2/24.
//

#pragma once

#include <common/memory/pooled_memory_resource.h>

#include <Kokkos_Core_fwd.hpp>

#include <functional>
#include <memory>
#include <optional>

namespace alps {

// Forward declaration
template<class ExecutionSpace>
class StreamPool;
namespace async {
// Forward declaration
class thread_pool;
template<class ExecutionSpace>
class event_pool;
} // namespace async
  //
namespace memory {
template<class BaseSpace>
class DynamicSizePool;
} // namespace memory

/**
 * @brief A class to manage the global resources associated with the application
 *
 * This class is a singleton that manages the global runtimes and resources,
 * such as Kokkos, MPI, and pooled resources. Cleanup functions can be
 * registered to be called when the context is finalized.
 *
 * The manager can be used in two ways:
 * 1. Explicitly obtain a scope guard as
 *
 * ```c++
 * {
 *   auto guard = RuntimeManager::ScopeGuard();
 *   // Do something with the context
 *   // The context is finalized when the guard goes out of scope
 * }
 * ```
 *
 *    The resources should be manually finalized by calling `finalize`.
 *
 * 2. Use `execute_main` to run a function with the context, e.g.
 *
 * ```c++
 * int some_function(int argc, char** argv) {
 *   std::cout << "Hello, world!" << std::endl;
 * }
 *
 * int main(int argc, char** argv) {
 *   return RuntimeManager::execute_main(some_function, argc, argv);
 * }
 * ```
 *
 * @note The class ensures that the resources are finalized in the correct
 * order.
 * @note There can be only one instance of the context.
 * @note The class is not thread-safe.
 */
class RuntimeManager
{
 public:
  class [[nodiscard]] ScopeGuard
  {
   public:
    ScopeGuard();
    ~ScopeGuard();

    ScopeGuard(ScopeGuard const&)            = delete;
    ScopeGuard& operator=(ScopeGuard const&) = delete;

   private:
    RuntimeManager& manager_;
  };

  /// Finalizes all runtimes and resources
  void finalize();

  /// @brief Returns an instance of the context, throws an exception if the
  /// context is not initialized
  [[nodiscard]] static RuntimeManager& instance()
  {
    static RuntimeManager instance;
    return instance;
  }

  /// @brief Executes a function with the context
  static int execute_main(std::function<int(int, char**)> const& func,
                          int                                    argc,
                          char**                                 argv);

  /// Initialize all runtimes: MPI and Kokkos
  void init_runtimes();

  /// @brief Add a cleanup function, which is called when the context is
  /// finalized
  void push_cleanup_function(std::function<void()> cleanup_function);

  StreamPool<Kokkos::DefaultExecutionSpace> const& stream_pool();

  Kokkos::DefaultExecutionSpace device_to_host_stream();

  Kokkos::DefaultExecutionSpace host_to_device_stream();

  async::thread_pool& thread_pool();

  async::event_pool<Kokkos::DefaultExecutionSpace>& event_pool();

  memory::PooledMemoryResources const& pooled_memory_resources();

  template<class BaseSpace>
  memory::DynamicSizePool<BaseSpace>& memory_pool()
  {
    return pooled_memory_resources().get<BaseSpace>();
  }

  RuntimeManager(const RuntimeManager&) = delete;

  RuntimeManager& operator=(const RuntimeManager&) = delete;

  RuntimeManager(RuntimeManager&&) = delete;

  RuntimeManager& operator=(RuntimeManager&&) = delete;

  ~RuntimeManager();

 private:
  RuntimeManager();

  void init_mpi();

  void init_kokkos();

  class Impl;
  std::unique_ptr<RuntimeManager::Impl> instance_;
  RuntimeManager::Impl&                 get_impl();

  std::optional<memory::PooledMemoryResources> pooled_mr_;
};
} // namespace alps
