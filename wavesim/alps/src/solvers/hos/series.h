#pragma once

#include <common/container/view_types.h>
#include <common/real_type.h>
#include <spectral/spectral_fwd.h>

namespace alps::solver::hos {

/** @brief Calculate the Taylor series coefficients for the surface elevation
 *
 * \[ Z_m = \eta \frac{\eta^m}{m!} \]
 *
 * @note This function is asynchronous, i.e. the coefficients are not guaranteed
 * to finish computing when the function returns
 */
void taylor_series_coeff_async(MDView<Real***>                      zp_hos,
                               MDView<const Real**>                 eta,
                               int                                  npw,
                               Grid const&                          grid,
                               Kokkos::DefaultExecutionSpace const& space);

/** @brief Calculate the expansion of the surface potential, r_hat */
void surface_vp_expansion(MDView<Real***>                      r_hat,
                          MDView<const Real**>                 vps,
                          MDView<const Real***>                zp_hos,
                          MDView<const Real***>                wvn_hos,
                          Grid const&                          grid,
                          Kokkos::DefaultExecutionSpace const& space);

/** @brief Calculate the vertical velocity at the surface */
void surface_w(MDView<Real**>                       ws,
               MDView<Real***>                      r_hat,
               MDView<const Real***>                zp_hos,
               MDView<const Real***>                wvn_hos,
               Grid const&                          grid,
               Kokkos::DefaultExecutionSpace const& space);

void calc_wavenumbers(MDView<Real***> wvn_hos,
                      Real            kx0,
                      Real            ky0,
                      Grid const&     grid);

} // namespace alps::solver::hos
