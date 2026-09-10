#pragma once

#include "spectral_base.h"

#include <common/base/logging_fwd.h>
#include <decomp/pencil_plan.h>
#include <fft/fftplan_fftw.h>

#include <Kokkos_Complex.hpp>
#include <Kokkos_Core.hpp>

namespace alps::spectral {

template<class T>
class SpectralPlan<T, Kokkos::OpenMP, fft::FFTW> final
  : public SpectralPlanBase<T, Kokkos::OpenMP>
{
 private:
  using base_t    = SpectralPlanBase<T, Kokkos::OpenMP>;
  using ExecSpace = Kokkos::OpenMP;
  using typename base_t::const_view_arg_t;
  using typename base_t::view_arg_t;

  using ComplexT = Kokkos::complex<T>;

 public:
  /// Initialize from a PencilPlan and base wavenumbers
  explicit SpectralPlan(PencilPlan const& partition, T kx0, T ky0);

  /// Initialize from another SpectralPlan's PencilPlan and base wavenumbers
  template<class RT, class RExecSpace, class RBackend>
  explicit SpectralPlan(
    const SpectralPlan<RT, RExecSpace, RBackend>& other_plan)
    : SpectralPlan{other_plan.pencil, other_plan.pex, other_plan.pey}
  {}

  using base_t::logger_;
  using base_t::pencil;
  using base_t::pex;
  using base_t::pey;

  void do_ddx(view_arg_t const&       output,
              const_view_arg_t const& input,
              ExecSpace const&        space) const override;
  void do_d2dx2(view_arg_t const&       output,
                const_view_arg_t const& input,
                ExecSpace const&        space) const override;
  void do_ddy(view_arg_t const&       output,
              const_view_arg_t const& input,
              SpectralPostOp          output_process,
              ExecSpace const&        space) const override;
  void do_d2dy2(view_arg_t const&       output,
                const_view_arg_t const& input,
                SpectralPostOp          output_process,
                ExecSpace const&        space) const override;
  void do_cutoff_xy(view_arg_t const& input,
                    int               kc_x,
                    int               kc_y,
                    ExecSpace const&  space) const override;
  void do_fft_r2c_xy(view_arg_t const&       output,
                     const_view_arg_t const& input,
                     ExecSpace const&        space) const override;
  void do_fft_c2r_xy(view_arg_t const& output,
                     view_arg_t const& input,
                     bool              dealias,
                     ExecSpace const&  space) const override;

  Kokkos::LayoutLeft get_r2c_xy_output_layout() const override;

  ~SpectralPlan() override;

 private:
  template<class Op>
  void do_x_derivatives(view_arg_t const&       output,
                        const_view_arg_t const& input,
                        Op                      derivative_func,
                        ExecSpace const&        space) const;
  template<class Op>
  void do_y_derivatives(view_arg_t const&       output,
                        const_view_arg_t const& input,
                        Op                      derivative_func,
                        SpectralPostOp          output_process,
                        ExecSpace const&        space) const;

  void tune_transpose_and_transform_chunk_size();

  using interface      = fft::FFTWInterface<T>;
  using fftw_complex_t = typename interface::complex_t;
  using plan_t         = typename interface::unique_ptr;

  plan_t plan_r2c_x;
  plan_t plan_r2c_y;
  plan_t plan_c2r_x;
  plan_t plan_c2r_y;

  int transpose_and_transform_chunk_size_z{0};
};

} // namespace alps::spectral
