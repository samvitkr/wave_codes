//
// Created by xuanx004 on 12/29/23.
//

#pragma once

#include <solvers/field/flow_field.h>

void apply_translation_to_fix_bottom_mean_u(
  alps::solver::FlowField const&       flow,
  Kokkos::DefaultExecutionSpace const& space);
