#include <solvers/field/flow_field.h>

#include <Kokkos_Core_fwd.hpp>

namespace alps {
namespace solver {

void apply_top_bc(const FlowField&                     flow,
                  const GradientWall&                  bc,
                  const Kokkos::DefaultExecutionSpace& stream);

void apply_bottom_bc(const FlowField&                     flow,
                     const GradientWall&                  bc,
                     const Kokkos::DefaultExecutionSpace& stream);

void apply_top_bc(const FlowField&                     flow,
                  const NoSlipWall&                    bc,
                  const Kokkos::DefaultExecutionSpace& stream);

void apply_bottom_bc(const FlowField&                     flow,
                     const NoSlipWall&                    bc,
                     const Kokkos::DefaultExecutionSpace& stream);

} // namespace solver
} // namespace alps
