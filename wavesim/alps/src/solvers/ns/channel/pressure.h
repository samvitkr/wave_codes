#pragma once

#include <common/container/view_types.h>
#include <common/real_type.h>

namespace alps {
namespace solver {

class FlowField;

template<typename ValueT>
class TridiagonalSolver;

class PressureEqn
{
 public:
  FlowField const& flow_field;
  MDView<Real***>  d, dl, du;

  std::unique_ptr<TridiagonalSolver<Real>> solver;

  explicit PressureEqn(const FlowField& flow_field);

  void initialize() const;

  ~PressureEqn();
};

void solve(const MDView<Real***>& div_u, const PressureEqn& peqn, Real dt);

} // namespace solver
} // namespace alps
