#pragma once

#include "field.h"
#include <common/program_options/config_table.h>

#include <spectral/spectral_fwd.h>

namespace alps::solver::hos {

class Smoother
{
 public:
  void apply(decltype(HOSField::value)            solution,
             Grid const&                          grid,
             Kokkos::DefaultExecutionSpace const& space) const;

  virtual std::string name() const;

  virtual std::string info() const;

  virtual ~Smoother();

 protected:
  virtual void apply_impl(decltype(HOSField::value)            solution,
                          Grid const&                          grid,
                          Kokkos::DefaultExecutionSpace const& space) const;
};

class LowPassFilter : public Smoother
{
 public:
  LowPassFilter(Real factor);

  std::string name() const override;

  std::string info() const override;

 private:
  void apply_impl(decltype(HOSField::value)            solution,
                  Grid const&                          grid,
                  Kokkos::DefaultExecutionSpace const& space) const override;

  Real factor_;
};

std::unique_ptr<Smoother> make_smoother(const alps::ConfigTable& config);

} // namespace alps::solver::hos
