#include "transposer_pool.h"

#include <common/base/logging.h>
#include <common/runtime/manager.h>

namespace alps::transpose {

template<typename T, typename ExecSpace>
typename TransposerPool<T, ExecSpace>::transposer_t const&
TransposerPool<T, ExecSpace>::emplace(const mpipp::communicator& comm,
                                      int                        n0,
                                      int                        n1,
                                      int                        max_nz,
                                      TransposerOptions          options)
{
  auto const key = TransposerKey{n0, n1, comm.size(), comm.raw_handle()};
  for (auto const& item : plans) {
    if (item.first == key) return item.second;
  }
  auto& new_item = plans.emplace_back(
    key, create_transposer<T, ExecSpace>(comm, n0, n1, max_nz, options));
  return new_item.second;
}

template<typename T, typename ExecSpace>
TransposerPool<T, ExecSpace>::TransposerPool()
{
  plans.reserve(8);
  RuntimeManager::instance().push_cleanup_function([&] { plans.clear(); });
}

template<typename T, typename ExecSpace>
TransposerPool<T, ExecSpace>::~TransposerPool()
{
  plans.clear();
}

template class TransposerPool<float, Kokkos::DefaultExecutionSpace>;
template class TransposerPool<double, Kokkos::DefaultExecutionSpace>;

#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
template class TransposerPool<float, Kokkos::DefaultHostExecutionSpace>;
template class TransposerPool<double, Kokkos::DefaultHostExecutionSpace>;
#endif
} // namespace alps::transpose
