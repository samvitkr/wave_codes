#include "spectral_base.h"

#include <Kokkos_Core.hpp>

namespace alps::spectral {
template<class T, class ExecSpace>
SpectralPlanBase<T, ExecSpace>::~SpectralPlanBase() = default;

template<class T, class ExecSpace>
SpectralPlanBase<T, ExecSpace>::SpectralPlanBase(
  PencilPlan const&              partition,
  T                              kx0,
  T                              ky0,
  FFTBackendVariant_t<ExecSpace> backend)
  : pencil{partition}
  , pex{kx0}
  , pey{ky0}
  , backend_type{backend}
  , logger_{get_logger("spectral")}
{}

template<class T, class ExecSpace>
void SpectralPlanBase<T, ExecSpace>::ddx(
  const SpectralPlanBase::view_arg_t&       output,
  const SpectralPlanBase::const_view_arg_t& input,
  const ExecSpace&                          space) const
{
  Kokkos::Tools::pushRegion("ddx");
  this->do_ddx(output, input, space);
  Kokkos::Tools::popRegion();
}

template<class T, class ExecSpace>
void SpectralPlanBase<T, ExecSpace>::d2dx2(
  const SpectralPlanBase::view_arg_t&       output,
  const SpectralPlanBase::const_view_arg_t& input,
  const ExecSpace&                          space) const
{
  Kokkos::Tools::pushRegion("d2dx2");
  this->do_d2dx2(output, input, space);
  Kokkos::Tools::popRegion();
}

template<class T, class ExecSpace>
void SpectralPlanBase<T, ExecSpace>::ddy(
  const SpectralPlanBase::view_arg_t&       output,
  const SpectralPlanBase::const_view_arg_t& input,
  SpectralPostOp                            output_process,
  const ExecSpace&                          space) const
{
  Kokkos::Tools::pushRegion("ddy");
  this->do_ddy(output, input, output_process, space);
  Kokkos::Tools::popRegion();
}

template<class T, class ExecSpace>
void SpectralPlanBase<T, ExecSpace>::d2dy2(
  const SpectralPlanBase::view_arg_t&       output,
  const SpectralPlanBase::const_view_arg_t& input,
  SpectralPostOp                            output_process,
  const ExecSpace&                          space) const
{
  Kokkos::Tools::pushRegion("d2dy2");
  this->do_d2dy2(output, input, output_process, space);
  Kokkos::Tools::popRegion();
}

template<class T, class ExecSpace>
void SpectralPlanBase<T, ExecSpace>::cutoff_xy(
  const SpectralPlanBase::view_arg_t& input,
  int                                 kc_x,
  int                                 kc_y,
  const ExecSpace&                    space) const
{
  Kokkos::Tools::pushRegion("cutoff");
  this->do_cutoff_xy(input, kc_x, kc_y, space);
  Kokkos::Tools::popRegion();
}

template<class T, class ExecSpace>
void SpectralPlanBase<T, ExecSpace>::fft_r2c_xy(
  const SpectralPlanBase::view_arg_t&       output,
  const SpectralPlanBase::const_view_arg_t& input,
  const ExecSpace&                          space) const
{
  Kokkos::Tools::pushRegion("r2c_xy");
  this->do_fft_r2c_xy(output, input, space);
  Kokkos::Tools::popRegion();
}

template<class T, class ExecSpace>
void SpectralPlanBase<T, ExecSpace>::fft_c2r_xy(
  const SpectralPlanBase::view_arg_t& output,
  const SpectralPlanBase::view_arg_t& input,
  bool                                dealias,
  const ExecSpace&                    space) const
{
  Kokkos::Tools::pushRegion("c2r_xy");
  this->do_fft_c2r_xy(output, input, dealias, space);
  Kokkos::Tools::popRegion();
}

// Explicit instantiation
template class SpectralPlanBase<double, Kokkos::DefaultHostExecutionSpace>;
template class SpectralPlanBase<float, Kokkos::DefaultHostExecutionSpace>;
#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
template class SpectralPlanBase<double, Kokkos::DefaultExecutionSpace>;
template class SpectralPlanBase<float, Kokkos::DefaultExecutionSpace>;
#endif
} // namespace alps::spectral
