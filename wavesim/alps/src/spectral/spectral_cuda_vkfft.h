#pragma once

#include "spectral_base.h"
#include <common/hash/hash.h>
#include <decomp/pencil_plan.h>
#include <fft/vkfft_structs.h>
#include <spectral/detail/vkfft_tuning.h>
#include <spectral/spectral_fwd.h>

#include <Kokkos_Complex.hpp>
#include <Kokkos_Core.hpp>

namespace alps::spectral {

template<class T>
class SpectralPlan<T, Kokkos::DefaultExecutionSpace, fft::VKFFT> final
  : public SpectralPlanBase<T, Kokkos::DefaultExecutionSpace>
{
 private:
  using base_t    = SpectralPlanBase<T, Kokkos::DefaultExecutionSpace>;
  using ExecSpace = Kokkos::DefaultExecutionSpace;
  using MemSpace  = Kokkos::DefaultExecutionSpace::memory_space;
  using typename base_t::const_view_arg_t;
  using typename base_t::view_arg_t;

  using ComplexT = Kokkos::complex<T>;

 public:
  explicit SpectralPlan(PencilPlan const&      partition,
                        T                      kx0,
                        T                      ky0,
                        SpectralOptions const& options = {});

  /// Initialize from another SpectralPlan's PencilPlan and base wavenumbers.
  /// Preserves the tuning-enabled state of the source plan so that an
  /// explicit Disabled (or Enabled) is not silently reset to Auto.
  template<class RT, class RExecSpace, class RBackend>
  explicit SpectralPlan(
    const SpectralPlan<RT, RExecSpace, RBackend>& other_plan)
    : SpectralPlan{other_plan.pencil, other_plan.pex, other_plan.pey}
  {}

  using base_t::logger_;
  using base_t::pencil;
  using base_t::pex;
  using base_t::pey;

  Kokkos::View<ComplexT**, Kokkos::LayoutLeft, MemSpace> buf_spectral;
  Kokkos::View<T***, Kokkos::LayoutLeft, MemSpace>       buf_y_physical;

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
  template<::alps::Pencil Direction>
  const auto& get_plan(int nz) const;

  template<::alps::Pencil Direction>
  const auto& get_cutoff_plan(int nz, int cutoff_idx) const;

  template<::alps::Pencil Direction>
  const auto& get_convolution_plan(int nz) const;

  template<::alps::Pencil Direction>
  auto emplace_plan(int nz) const;

  template<::alps::Pencil Direction>
  auto emplace_cutoff_plan(int nz, int cutoff_idx) const;

  template<::alps::Pencil Direction>
  auto emplace_convolution_plan(int nz) const;

  template<::alps::Pencil Direction>
  VkFFTConfiguration get_plan_base_config() const;

  void apply_x_convolution(view_arg_t const&       output,
                           const_view_arg_t const& input,
                           ComplexT*               kernel,
                           ExecSpace const&        space) const;

  void apply_y_convolution(view_arg_t const&       output,
                           const_view_arg_t const& input,
                           ComplexT*               kernel,
                           SpectralPostOp          output_process,
                           ExecSpace const&        space) const;

  // VkFFT configuration needs the sizes and addresses of buffers stored to take
  // address
  using VkApp = std::unique_ptr<VkFFTApplication, void (*)(VkFFTApplication*)>;
  struct VKFFTPlanWithStorage
  {
    VkApp plan;
#if defined(KOKKOS_ENABLE_CUDA)
    CUdevice device{};
#elif defined(KOKKOS_ENABLE_HIP)
    hipDevice_t device{};
#endif
    uint64_t buffer_size{};
    uint64_t kernel_size{};
    void*    buffer{nullptr};
  };

  struct VkFFTPlanKey
  {
    int size{0};  // transform size
    int batch{0}; // transform batch
    int pad{0};

    constexpr bool operator==(const VkFFTPlanKey& rhs) const noexcept
    {
      return (size == rhs.size) && (batch == rhs.batch) && (pad == rhs.pad);
    }

    constexpr bool operator!=(const VkFFTPlanKey& rhs) const noexcept
    {
      return !(*this == rhs);
    }
  };

  struct VkFFTPlanKeyHash
  {
    size_t operator()(const VkFFTPlanKey& p) const
    {
      size_t seed = 0;
      ::alps::hash_ext::hash_combine(seed, p.size);
      ::alps::hash_ext::hash_combine(seed, p.batch);
      ::alps::hash_ext::hash_combine(seed, p.pad);
      return seed;
    }
  };

  using plan_map =
    std::unordered_map<VkFFTPlanKey, VKFFTPlanWithStorage, VkFFTPlanKeyHash>;

  mutable plan_map plans;
  mutable plan_map conv_plans;
  mutable plan_map cutoff_plans;

  Kokkos::View<ComplexT* [2], Kokkos::LayoutLeft, MemSpace>
    kernel_ddx; // Kernels for derivatives up to 2nd order
  Kokkos::View<ComplexT* [2], Kokkos::LayoutLeft, MemSpace>
    kernel_ddy; // Kernels for derivatives up to 2nd order
  Kokkos::View<ComplexT*, Kokkos::LayoutLeft, MemSpace>
    kernel_cutoff_x; // Kernels for cutoff in x
  Kokkos::View<ComplexT*, Kokkos::LayoutLeft, MemSpace>
    kernel_cutoff_y; // Kernels for cutoff in y

  detail::VkFFTTuningConfig tuning_config_;
};

} // namespace alps::spectral
