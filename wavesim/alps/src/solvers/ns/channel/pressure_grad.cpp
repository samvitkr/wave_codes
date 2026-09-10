#include <solvers/source_terms/pressure_grad.h>

#include <common/container/view_utils.h>
#include <common/kokkos_abstraction/exec_policy.h>
#include <common/runtime/async_utils.h>
#include <solvers/ns/channel/solver_ab2.h>
#include <solvers/ns/channel/solver_ab2cn.h>

namespace alps {
namespace solver {

namespace {
void add_pressure_gradient(Vector3Field<Real***> const& source,
                           FlowField const& /*flow*/,
                           std::array<Real, 2> gradients)
{
  std::vector<decltype(get_next_stream())> streams;

  for (const auto& [flux, gradient] :
       {tie(source.x, gradients[0]), tie(source.y, gradients[1])}) {
    if (std::abs(gradient) < std::numeric_limits<Real>::epsilon()) continue;

    const auto f  = create_inner_view(flux);
    const auto dp = gradient;
    Kokkos::parallel_for(
      LoopPolicy<3>(streams.emplace_back(get_next_stream()),
                    local_begins(f),
                    local_ends(f)),
      KOKKOS_LAMBDA(int i, int j, int k) { f(i, j, k) += dp; });
  }

  fence(std::move(streams));
}
} // anonymous namespace

template<>
void ConstantPressureGradient<ChannelFlowSolverAB2>::add_forces(
  const Vector3Field<Real***>&         fu,
  const Kokkos::DefaultExecutionSpace& space) const
{
  space.fence();
  add_pressure_gradient(fu, solver_.flow_field, gradients_);
}

template<>
void ConstantPressureGradient<ChannelFlowSolverAB2CN>::add_forces(
  const Vector3Field<Real***>&         fu,
  const Kokkos::DefaultExecutionSpace& space) const
{
  space.fence();
  add_pressure_gradient(fu, solver_.flow_field, gradients_);
}

} // namespace solver
} // namespace alps
