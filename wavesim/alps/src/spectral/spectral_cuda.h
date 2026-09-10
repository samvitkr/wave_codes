#pragma once

#include "spectral_base.h"
#include <common/base/logging_fwd.h>
#include <decomp/pencil_plan.h>
#include <fft/fftplan_cufft.h>

#include <Kokkos_Complex.hpp>
#include <Kokkos_Core.hpp>

namespace alps::spectral {

template<class T>
class SpectralPlan<T, Kokkos::Cuda, fft::CUFFT> final
  : public SpectralPlanBase<T, Kokkos::Cuda>
{
 private:
  using base_t    = SpectralPlanBase<T, Kokkos::Cuda>;
  using ExecSpace = Kokkos::Cuda;
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

  Kokkos::View<ComplexT***, Kokkos::LayoutLeft, Kokkos::CudaSpace>
    buf_x_spectral;
  Kokkos::View<ComplexT***, Kokkos::LayoutLeft, Kokkos::CudaSpace>
                                                            buf_y_spectral;
  Kokkos::View<T***, Kokkos::LayoutLeft, Kokkos::CudaSpace> buf_y_physical;

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
  template<class Tag>
  void do_x_derivatives(view_arg_t const&       output,
                        const_view_arg_t const& input,
                        Tag const& /*tag*/,
                        Kokkos::Cuda const& space) const;

  template<class Tag>
  void do_y_derivatives(view_arg_t const&       output,
                        const_view_arg_t const& input,
                        SpectralPostOp          output_process,
                        Tag const& /*tag*/,
                        Kokkos::Cuda const& space) const;

  template<Pencil Direction>
  fft::CUFFTPlan<T, fft::R2C> const& get_plan(int nz, fft::R2C /*tag*/) const;

  template<Pencil Direction>
  fft::CUFFTPlan<T, fft::C2R> const& get_plan(int nz, fft::C2R /*tag*/) const;

  void setup_plan(int nz, Pencil direction, fft::R2C /*tag*/) const;

  void setup_plan(int nz, Pencil direction, fft::C2R /*tag*/) const;

  mutable std::unordered_map<int, fft::CUFFTPlan<T, fft::R2C>> plans_x_r2c;
  mutable std::unordered_map<int, fft::CUFFTPlan<T, fft::R2C>> plans_y_r2c;
  mutable std::unordered_map<int, fft::CUFFTPlan<T, fft::C2R>> plans_x_c2r;
  mutable std::unordered_map<int, fft::CUFFTPlan<T, fft::C2R>> plans_y_c2r;
};

} // namespace alps::spectral
