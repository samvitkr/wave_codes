#include "spectral.h"

#include <decomp/pencil_plan.h>

#include <utility>
#include <variant>

namespace alps::spectral {

SpectralGrid::SpectralGrid(PencilPlan const& partition, double kx0, double ky0)
  : SpectralGrid(partition, kx0, ky0, {}, {})
{}

SpectralGrid::SpectralGrid(
  PencilPlan const&                                  partition,
  double                                             kx0,
  double                                             ky0,
  FFTBackendVariant_t<Kokkos::DefaultExecutionSpace> dev_backend)
  : SpectralGrid(partition, kx0, ky0, std::move(dev_backend), {})
{}

SpectralGrid::SpectralGrid(PencilPlan const& partition,
                           double            kx0,
                           double            ky0,
                           SpectralOptions   options)
  : SpectralGrid(partition, kx0, ky0, {}, std::move(options))
{}

SpectralGrid::SpectralGrid(
  PencilPlan const&                                  partition,
  double                                             kx0,
  double                                             ky0,
  FFTBackendVariant_t<Kokkos::DefaultExecutionSpace> dev_backend,
  SpectralOptions                                    options)
  : pencil{partition}
  , pex{kx0}
  , pey{ky0}
  , dev_backend_{std::move(dev_backend)}
  , options_{std::move(options)}
{}

template<class T, class ExecSpace>
void SpectralGrid::emplace_plan() const
{
  auto& plan_storage =
    get_plan_storage(alps::is_default_execution_space<ExecSpace>(), T{});
  if constexpr (alps::is_default_execution_space_v<ExecSpace>) {
    std::visit(
      [&plan_storage, this](auto&& backend) {
        using Backend = std::decay_t<decltype(backend)>;
        plan_storage =
          std::move(SpectralPlanFactory::create<T, ExecSpace, Backend>(
            pencil, (T)pex, (T)pey, options_));
      },
      dev_backend_);
  } else if constexpr (alps::is_host_execution_space_v<ExecSpace>) {
    std::visit(
      [&plan_storage, this](auto&& backend) {
        using Backend = std::decay_t<decltype(backend)>;
        plan_storage =
          std::move(SpectralPlanFactory::create<T, ExecSpace, Backend>(
            pencil, (T)pex, (T)pey, options_));
      },
      host_backend_);
  }
}

SpectralGrid::~SpectralGrid() = default;

template void SpectralGrid::emplace_plan<double, Kokkos::OpenMP>() const;
template void SpectralGrid::emplace_plan<float, Kokkos::OpenMP>() const;

#if defined(KOKKOS_ENABLE_CUDA)
template void SpectralGrid::emplace_plan<double, Kokkos::Cuda>() const;
template void SpectralGrid::emplace_plan<float, Kokkos::Cuda>() const;
#elif defined(KOKKOS_ENABLE_HIP)
template void SpectralGrid::emplace_plan<double, Kokkos::HIP>() const;
template void SpectralGrid::emplace_plan<float, Kokkos::HIP>() const;
#endif
} // namespace alps::spectral
