#pragma once

#include "field.h"
#include "smooth.h"
#include "solver_rk_stage.h"
#include <common/base/logging_fwd.h>
#include <common/container/view_types.h>
#include <common/container/view_utils.h>
#include <common/program_options/config_table.h>

#include <filesystem>

namespace alps::solver::hos {

class HOSSolver
{
 public:
  using PaCallBackFcn = std::function<
    void(decltype(HOSField::pa), double, Kokkos::DefaultExecutionSpace const&)>;

  HOSSolver(const HOSField& wave_field,
            Real            Fr2_,
            Real            We_,
            int             expansion_order);

  HOSSolver(const HOSField& wave_field, ConfigTable const& config);

  void set_smoother(std::unique_ptr<Smoother> smoother);

  /// advance the solution by dt
  void do_step(Real dt, PaCallBackFcn const& pa_callback = {}) const;

  /// advance the solution by dt using RK4
  void do_step(RK4 const&           integrator,
               Real                 dt,
               PaCallBackFcn const& pa_callback = {}) const;

  /// advance the solution by dt using LowStorageRK2N
  void do_step(LowStorageRK2N const& integrator,
               Real                  dt,
               PaCallBackFcn const&  pa_callback = {}) const;

  /// advance the solution by dt using LowStorageRK2C
  void do_step(LowStorageRK2C const& integrator,
               Real                  dt,
               PaCallBackFcn const&  pa_callback = {}) const;

  void update_rhs(MDView<Real** [2]> const&            dFdt,
                  MDView<Real const** [2]> const&      sol,
                  MDView<Real const**> const&          pa,
                  Kokkos::DefaultExecutionSpace const& space) const;

  void get_ws(MDView<Real**> const&                ws,
              MDView<Real const** [2]> const&      sol,
              Kokkos::DefaultExecutionSpace const& space) const;

  double get_time() const { return solution.time; }

  void set_time(double t) const { solution.time = t; }

  /// maximum angular frequency in the domain estimated using linear wave theory
  Real estimate_max_omega() const;

  void save(std::filesystem::path filename) const; // save solution to file

  void load(std::filesystem::path filename); // read solution

  HOSField const& solution;
  Grid const&     grid;

  MDView<Real***> wvn_hos;

  Real Fr2;
  Real We;
  int  order;

  RKIntegrator integrator_;

  std::unique_ptr<Smoother> smoother_;

 private:
  static std::unique_ptr<Smoother> default_smoother();

  Logger logger;
};

std::ostream& operator<<(std::ostream& os, HOSSolver const& solver);

} // namespace alps::solver::hos
