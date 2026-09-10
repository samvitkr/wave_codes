//
// Created by xuananqing on 3/26/23.
//
#pragma once

#include "spectral_fwd.h"
#include <common/base/logging_fwd.h>
#include <decomp/pencil_plan.h>
#include <fft/fftplan_fwd.h>

#include <Kokkos_Core_fwd.hpp>

#include <memory>
#include <variant>

// Forward declaration
namespace Kokkos {
struct LayoutLeft;
} // namespace Kokkos

namespace alps::spectral {

template<class ExecSpace>
struct FFTBackendVariant;

template<class ExecSpace>
using FFTBackendVariant_t = typename FFTBackendVariant<ExecSpace>::type;

template<>
struct FFTBackendVariant<Kokkos::OpenMP>
{
  using type = std::variant<::alps::fft::FFTW>;
};

#if defined(KOKKOS_ENABLE_CUDA)
template<>
struct FFTBackendVariant<Kokkos::Cuda>
{
  using type = std::variant<::alps::fft::VKFFT, ::alps::fft::CUFFT>;
};
#elif defined(KOKKOS_ENABLE_HIP)
template<>
struct FFTBackendVariant<Kokkos::HIP>
{
  using type = std::variant<::alps::fft::VKFFT>;
};
#endif

template<class T, class ExecSpace>
class SpectralPlanBase
{
 protected:
  using view_arg_t = Kokkos::View<T***, Kokkos::LayoutLeft, ExecSpace>;
  using const_view_arg_t =
    Kokkos::View<T const***, Kokkos::LayoutLeft, ExecSpace>;

 public:
  void ddx(view_arg_t const&       output,
           const_view_arg_t const& input,
           ExecSpace const&        space) const;

  void d2dx2(view_arg_t const&       output,
             const_view_arg_t const& input,
             ExecSpace const&        space) const;

  void ddy(view_arg_t const&       output,
           const_view_arg_t const& input,
           SpectralPostOp          output_process,
           ExecSpace const&        space) const;

  void d2dy2(view_arg_t const&       output,
             const_view_arg_t const& input,
             SpectralPostOp          output_process,
             ExecSpace const&        space) const;

  void cutoff_xy(view_arg_t const& input,
                 int               kc_x,
                 int               kc_y,
                 ExecSpace const&  space) const;

  void fft_r2c_xy(view_arg_t const&       output,
                  const_view_arg_t const& input,
                  ExecSpace const&        space) const;

  void fft_c2r_xy(view_arg_t const& output,
                  view_arg_t const& input,
                  bool              dealias,
                  ExecSpace const&  space) const;

  virtual void do_ddx(view_arg_t const&       output,
                      const_view_arg_t const& input,
                      ExecSpace const&        space) const = 0;

  virtual void do_d2dx2(view_arg_t const&       output,
                        const_view_arg_t const& input,
                        ExecSpace const&        space) const = 0;

  virtual void do_ddy(view_arg_t const&       output,
                      const_view_arg_t const& input,
                      SpectralPostOp          output_process,
                      ExecSpace const&        space) const = 0;

  virtual void do_d2dy2(view_arg_t const&       output,
                        const_view_arg_t const& input,
                        SpectralPostOp          output_process,
                        ExecSpace const&        space) const = 0;

  virtual void do_cutoff_xy(view_arg_t const& input,
                            int               kc_x,
                            int               kc_y,
                            ExecSpace const&  space) const = 0;

  virtual void do_fft_r2c_xy(view_arg_t const&       output,
                             const_view_arg_t const& input,
                             ExecSpace const&        space) const = 0;

  virtual void do_fft_c2r_xy(view_arg_t const& output,
                             view_arg_t const& input,
                             bool              dealias,
                             ExecSpace const&  space) const = 0;

  virtual Kokkos::LayoutLeft get_r2c_xy_output_layout() const = 0;

  PencilPlan pencil;
  T          pex;
  T          pey;

  FFTBackendVariant_t<ExecSpace> backend_type;

  virtual ~SpectralPlanBase();

  // make this class non-copyable
  SpectralPlanBase(SpectralPlanBase const&)            = delete;
  SpectralPlanBase& operator=(SpectralPlanBase const&) = delete;
  // make this class non-movable
  SpectralPlanBase(SpectralPlanBase&&)            = delete;
  SpectralPlanBase& operator=(SpectralPlanBase&&) = delete;

 protected:
  explicit SpectralPlanBase(PencilPlan const&              partition,
                            T                              kx0,
                            T                              ky0,
                            FFTBackendVariant_t<ExecSpace> backend);

  Logger logger_;
};

struct SpectralOptions
{
  VkFFTTuningControl vkfft_tuning{VkFFTTuningControl::Auto};
};

class SpectralPlanFactory
{
 public:
  template<class T, class ExecSpace, class Backend>
  static std::unique_ptr<SpectralPlanBase<T, ExecSpace>>
  create(PencilPlan const&      partition,
         T                      kx0,
         T                      ky0,
         SpectralOptions const& options = {});
};

} // namespace alps::spectral
