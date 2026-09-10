//
// Created by xuananqing on 7/4/23.
//

#pragma once

#include <common/container/vector_field.h>
#include <solvers/field/bc_types.h>
#include <solvers/mesh/mesh_fwd.h>

#include <common/container/view_types.h>

namespace alps::solver {

class FreeSurfaceFlowField;

struct FreeSurfaceBCOptions
{
  double hyper_viscosity_c{1};
  int    hyper_viscosity_n{3};

  static FreeSurfaceBCOptions parse_from(ConfigTable const& config);
};

struct FreeSurfaceBC : public VelocityBC
{
  std::string info() const override;

  FreeSurfaceBC(TopWaveMesh const& mesh, FreeSurfaceBCOptions const& options_);

  void update_eta_t(TopWaveMesh const&          mesh,
                    MDView<Real const**> const& u_s,
                    MDView<Real const**> const& v_s,
                    MDView<Real const**> const& w_s) const;

  void get_updated_eta(MDView<Real**> const&       eta,
                       MDView<Real const**> const& u_s,
                       MDView<Real const**> const& v_s,
                       MDView<Real const**> const& w_s,
                       double                      nu,
                       Real                        dt,
                       int                         rk_stage,
                       TopWaveMesh const&          mesh) const;

  /**
   * Compute the velocity gradient at the surface: d(J^{-1}u)/dz, d(J^{-1}v)/dz,
   * d(J^{-1}w)/dz
   *
   * dw/dz is assumed to be unchanged across steps
   * @note vec_u is (u),  vec_uz is (d(J^{-1} u)/dzeta)
   */
  void get_surface_uz(MDView<Real** [3]> const&    vec_invJuz,
                      Vector3Field<Real***> const& vec_u,
                      Real                         nu,
                      TopWaveMesh const&           mesh) const;

  /**
   * Given pressure and velocity gradient, estimate the surface velocity
   * J^{-1}\hat{u}=J^{-1}u+Δt*(J^{-1}∇p)
   *
   * @note vec_u is (J^{-1}u)
   */
  void get_surface_uhat(Vector3Field<Real***> const&    vec_invJu,
                        MDView<Real const** [3]> const& vec_invJuz,
                        HaloView<Real const***> const&  pp,
                        Real                            dt,
                        TopWaveMesh const&              mesh) const;

  /**
   * Compute J^{-1}w at the surface given J^{-1}u and J^{-1}v
   */
  static void get_surface_invJw(MDView<Real** [2]> const&            invJw,
                                MDView<Real const** [3]> const&      invJu,
                                MDView<Real const** [3]> const&      invJv,
                                TopWaveMesh const&                   mesh,
                                Kokkos::DefaultExecutionSpace const& stream);

  void get_surface_pressure(MDView<Real**> const&       p_top,
                            MDView<Real const**> const& eta,
                            Real                        Fr2,
                            Real                        RWe,
                            Real                        nu,
                            TopWaveMesh const&          mesh,
                            bool update_stored_p_top) const;

  FreeSurfaceBCOptions options;

  /// Shear stress parallel to the surface in x-direction
  MDView<Real**> tau_1;
  /// Shear stress parallel to the surface in y-direction
  MDView<Real**> tau_2;
  /// Estimated velocity u
  MDView<Real** [3]> ue_;
  /// Estimated velocity v
  MDView<Real** [3]> ve_;
  /// Estimated velocity w
  MDView<Real** [2]> we_;
  /// pressure at the surface at the t_n
  MDView<Real**> p_top_0;
};

} // namespace alps::solver
