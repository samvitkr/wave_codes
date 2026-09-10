//
// Created by xuanx004 on 7/12/24.
//

#pragma once

#include <common/program_options/config_table.h>
#include <common/real_type.h>

#include <memory>
#include <string>

namespace alps::solver {

struct ScalarBC
{
  virtual std::string info() const = 0;

  virtual std::unique_ptr<ScalarBC> get_updated_bc(double time);

  virtual void load_from(const ConfigTable& config);

  virtual ~ScalarBC() = default;

  bool is_time_dependent{false};

 protected:
  ScalarBC() = default;

  explicit ScalarBC(bool is_time_dependent_)
    : is_time_dependent(is_time_dependent_)
  {}
};

struct ConstantDirichletBC : public ScalarBC
{
  std::string info() const override;

  void load_from(const ConfigTable& config) override;

  Real value_{0};
};

struct ConstantGradientBC : public ScalarBC
{
  std::string info() const override;

  void load_from(const ConfigTable& config) override;

  Real grad_{0}; // dc/dz
};

struct ConstantFluxBC : public ScalarBC
{
  std::string info() const override;

  void load_from(const ConfigTable& config) override;

  Real flux_{0}; // flux in the +z direction
};

} // namespace alps::solver
