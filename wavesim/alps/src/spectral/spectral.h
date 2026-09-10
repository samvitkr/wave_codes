#pragma once

#include "spectral_base.h"
#include <common/device/device_traits.h>
#include <decomp/pencil_plan.h>

#include <memory>
#include <type_traits>
#include <utility>

namespace alps::spectral {

struct SpectralGrid
{
 private:
  template<class T, class ExecSpace>
  using storage_t = std::unique_ptr<SpectralPlanBase<T, ExecSpace>>;

 public:
  SpectralGrid(PencilPlan const&                                  partition,
               double                                             kx0,
               double                                             ky0,
               FFTBackendVariant_t<Kokkos::DefaultExecutionSpace> dev_backend);

  explicit SpectralGrid(PencilPlan const& partition, double kx0, double ky0);

  explicit SpectralGrid(PencilPlan const& partition,
                        double            kx0,
                        double            ky0,
                        SpectralOptions   options);

  SpectralGrid(PencilPlan const&                                  partition,
               double                                             kx0,
               double                                             ky0,
               FFTBackendVariant_t<Kokkos::DefaultExecutionSpace> dev_backend,
               SpectralOptions                                    options);

  /// Return the corresponding pencil partition
  const auto& partition(Pencil which = Pencil::X) const noexcept
  {
    return which == Pencil::Y ? pencil.y_pencil : pencil.x_pencil;
  }

  const auto& x_pencil() const noexcept { return pencil.x_pencil; }

  const auto& y_pencil() const noexcept { return pencil.y_pencil; }

  const auto& comm() const noexcept { return pencil.comm(); }

  /// Size of one dimension of the local block
  auto extent(int axis, Pencil which = Pencil::X) const noexcept
  {
    return partition(which).extents[axis];
  }

  /// Sizes of the local block in all dimensions
  const auto& extents(Pencil which = Pencil::X) const noexcept
  {
    return partition(which).extents;
  }

  /// Global size
  auto global_extent(int axis, Pencil which = Pencil::X) const noexcept
  {
    return partition(which).global_extents[axis];
  }

  /// Global index offset of the local block for a specified axis
  auto offset(int axis, Pencil which = Pencil::X) const noexcept
  {
    return partition(which).offsets[axis];
  }

  /// Global index offsets of the local block in all dimensions
  const auto& offsets(Pencil which = Pencil::X) const noexcept
  {
    return partition(which).offsets;
  }

  /// Get the spectral plan for a given floating-point type and execution space
  template<class T, class ExecSpace>
  const SpectralPlanBase<T, ExecSpace>& get_plan() const
  {
    auto& plan_storage =
      get_plan_storage(alps::is_default_execution_space<ExecSpace>(), T{});
    if (plan_storage == nullptr) {
      emplace_plan<T, ExecSpace>();
    }
    return *plan_storage;
  }

  /// Get the layout of the spectral space, can serve as output of fft_r2c_xy
  template<class T, class ExecSpace>
  auto get_r2c_xy_output_layout() const
  {
    const auto& plan = get_plan<T, ExecSpace>();
    return plan.get_r2c_xy_output_layout();
  }

  ~SpectralGrid();

  PencilPlan pencil;

  double pex;
  double pey;

 private:
  auto& get_plan_storage(std::true_type /*is_default_execspace*/,
                         double /*T*/) const
  {
    return plan_double;
  }

  auto& get_plan_storage(std::true_type /*is_default_execspace*/,
                         float /*T*/) const
  {
    return plan_float;
  }

  auto& get_plan_storage(std::false_type /*is_default_execspace*/,
                         double /*T*/) const
  {
    return host_plan_double;
  }

  auto& get_plan_storage(std::false_type /*is_default_execspace*/,
                         float /*T*/) const
  {
    return host_plan_float;
  }

  template<class T, class ExecSpace>
  void emplace_plan() const;

  // actual storage for spectral plans of different types and execution spaces
  // (mutable for lazy initialization)
  mutable storage_t<double, Kokkos::DefaultExecutionSpace> plan_double{};
  mutable storage_t<float, Kokkos::DefaultExecutionSpace>  plan_float{};
  mutable storage_t<double, Kokkos::DefaultHostExecutionSpace>
    host_plan_double{};
  mutable storage_t<float, Kokkos::DefaultHostExecutionSpace> host_plan_float{};

  FFTBackendVariant_t<Kokkos::DefaultExecutionSpace>     dev_backend_{};
  FFTBackendVariant_t<Kokkos::DefaultHostExecutionSpace> host_backend_{};

  SpectralOptions options_{};
};

/// @brief Compute the derivative of `input` along the x-axis and write to
/// `output`
template<class ExecSpace, class InType, class OutType>
void ddx(OutType&&           output,
         InType&&            input,
         const SpectralGrid& grid,
         ExecSpace const&    space)
{
  using T          = typename std::decay_t<InType>::non_const_value_type;
  const auto& plan = grid.get_plan<T, ExecSpace>();
  plan.ddx(std::forward<OutType>(output), std::forward<InType>(input), space);
}

/// @brief Compute the second derivative of `input` along the x-axis and write
/// to `output`
template<class ExecSpace, class InType, class OutType>
void d2dx2(OutType&&           output,
           InType&&            input,
           const SpectralGrid& grid,
           ExecSpace const&    space)
{
  using T          = typename std::decay_t<InType>::non_const_value_type;
  const auto& plan = grid.get_plan<T, ExecSpace>();
  plan.d2dx2(std::forward<OutType>(output), std::forward<InType>(input), space);
}

/// @brief Compute the derivative of `input` along the y-axis and write to
/// `output`
template<class ExecSpace, class InType, class OutType>
void ddy(OutType&&           output,
         InType&&            input,
         const SpectralGrid& grid,
         ExecSpace const&    space)
{
  using T          = typename std::decay_t<InType>::non_const_value_type;
  const auto& plan = grid.get_plan<T, ExecSpace>();
  plan.ddy(std::forward<OutType>(output),
           std::forward<InType>(input),
           SpectralPostOp::AssignAfterTranspose,
           space);
}

/// @brief Compute the second derivative of `input` along the y-axis and write
/// to `output`
template<class ExecSpace, class InType, class OutType>
void d2dy2(OutType&&           output,
           InType&&            input,
           const SpectralGrid& grid,
           ExecSpace const&    space)
{
  using T          = typename std::decay_t<InType>::non_const_value_type;
  const auto& plan = grid.get_plan<T, ExecSpace>();
  plan.d2dy2(std::forward<OutType>(output),
             std::forward<InType>(input),
             SpectralPostOp::AssignAfterTranspose,
             space);
}

/// @brief Compute the derivative of `input` along the y-axis and add to
/// `output`
template<class ExecSpace, class InType, class OutType>
void ddy_and_add(OutType&&           output,
                 InType&&            input,
                 const SpectralGrid& grid,
                 ExecSpace const&    space)
{
  using T          = typename std::decay_t<InType>::non_const_value_type;
  const auto& plan = grid.get_plan<T, ExecSpace>();
  plan.ddy(std::forward<OutType>(output),
           std::forward<InType>(input),
           SpectralPostOp::AddAfterTranspose,
           space);
}

/// @brief Compute the second derivative of `input` along the y-axis and add to
/// `output`
template<class ExecSpace, class InType, class OutType>
void d2dy2_and_add(OutType&&           output,
                   InType&&            input,
                   const SpectralGrid& grid,
                   ExecSpace const&    space)
{
  using T          = typename std::decay_t<InType>::non_const_value_type;
  const auto& plan = grid.get_plan<T, ExecSpace>();
  plan.d2dy2(std::forward<OutType>(output),
             std::forward<InType>(input),
             SpectralPostOp::AddAfterTranspose,
             space);
}

/// @brief Low pass filtering `input` in the spectral space, cutting off wave
/// numbers larger than `kc_x` and `kc_y`
template<class ExecSpace, class InType>
void cutoff_xy(InType&&            input,
               int                 kc_x,
               int                 kc_y,
               const SpectralGrid& grid,
               ExecSpace const&    space)
{
  using T          = typename std::decay_t<InType>::non_const_value_type;
  const auto& plan = grid.get_plan<T, ExecSpace>();
  plan.cutoff_xy(std::forward<InType>(input), kc_x, kc_y, space);
}

/// @brief Dealias `input` by using the 2/3 rule
template<class ExecSpace, class InType>
void dealias(InType&& input, const SpectralGrid& grid, ExecSpace const& space)
{
  cutoff_xy(std::forward<InType>(input),
            grid.global_extent(0) / 3,
            grid.global_extent(1) / 3,
            grid,
            space);
}

/// @brief Forward transform `input` from the physical space to the spectral
/// space
/**
 *  Transform is done by first performing R2C FFT along the x-axis, then treat
 * the complex arrays as real arrays and transpose, and finally perform R2C FFT
 * along the y-axis.
 *
 *  The output array layout may be obtained by calling
 * `grid.get_r2c_xy_output_layout<T, ExecSpace>()`.
 */
template<class ExecSpace, class InType, class OutType>
void fft_r2c_xy(OutType&&           output,
                InType&&            input,
                const SpectralGrid& grid,
                ExecSpace const&    space)
{
  using T          = typename std::decay_t<InType>::non_const_value_type;
  const auto& plan = grid.get_plan<T, ExecSpace>();
  plan.fft_r2c_xy(
    std::forward<OutType>(output), std::forward<InType>(input), space);
}

/// @brief Backward transform `input` from the spectral space to the physical
/// space
/**
 *  The reverse of `fft_r2c_xy`.
 */
template<class ExecSpace, class InType, class OutType>
void fft_c2r_xy(OutType&&           output,
                InType&&            input,
                bool                dealias,
                const SpectralGrid& grid,
                ExecSpace const&    space)
{
  using T          = typename std::decay_t<InType>::non_const_value_type;
  const auto& plan = grid.get_plan<T, ExecSpace>();
  plan.fft_c2r_xy(
    std::forward<OutType>(output), std::forward<InType>(input), dealias, space);
}

} // namespace alps::spectral
