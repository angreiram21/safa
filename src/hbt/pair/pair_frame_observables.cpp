/**
 * @file pair_frame_observables.cpp
 * @brief HBT Lab/LCMS/PRF pair separation observable implementation.
 */

#include "hbt/pair/pair_frame_observables.h"

#include "hbt/boosts/lab_to_lcms.h"
#include "hbt/boosts/lcms_to_prf.h"

#include <cmath>

namespace hbt {
namespace {

/**
 * @brief Minimal OSL projection of one LCMS relative separation.
 */
struct LCMSOSLProjection {
    double delta_t_fm;  ///< Relative emission time in LCMS, in fm.
    double r_out_fm;    ///< Out component in LCMS, in fm.
    double r_side_fm;   ///< Side component in LCMS, in fm.
    double r_long_fm;   ///< Long component in LCMS, in fm.
};

/**
 * @brief Project one LCMS relative separation onto the HBT OSL basis.
 * @param relative_separation_lcms_fm Relative separation in LCMS.
 * @param particle_a First accepted particle in canonical pair order.
 * @param particle_b Second accepted particle in canonical pair order.
 * @param kinematics Validated kinematics for the same ordered pair.
 * @return Temporal component and signed OSL spatial components.
 *
 * @pre kinematics.kt_gev is finite and non-negative.
 * @pre If kinematics.kt_gev == 0, both accepted particles have pT > 0 and
 *      therefore qT = |pa_T - pb_T| is strictly positive.
 *
 * For kT > 0, out is K_T / kT. For exact kT == 0, qT is calculated only in
 * that branch and defines out. side = ez x out = (-out_y, out_x, 0), so no
 * second normalization is needed. long is the beam z direction.
 */
LCMSOSLProjection project_to_osl(
    const common::FourVector& relative_separation_lcms_fm,
    const Particle& particle_a,
    const Particle& particle_b,
    const PairKinematics& kinematics
) noexcept {
    double out_x = 0.0;
    double out_y = 0.0;

    if (kinematics.kt_gev == 0.0) {
        const double qx_gev =
            particle_a.momentum.x1 - particle_b.momentum.x1;
        const double qy_gev =
            particle_a.momentum.x2 - particle_b.momentum.x2;
        const double qt_gev = std::hypot(qx_gev, qy_gev);
        out_x = qx_gev / qt_gev;
        out_y = qy_gev / qt_gev;
    } else {
        out_x = kinematics.kx_gev / kinematics.kt_gev;
        out_y = kinematics.ky_gev / kinematics.kt_gev;
    }

    const double delta_x_fm = relative_separation_lcms_fm.x1;
    const double delta_y_fm = relative_separation_lcms_fm.x2;

    return {
        relative_separation_lcms_fm.x0,
        delta_x_fm * out_x + delta_y_fm * out_y,
        -delta_x_fm * out_y + delta_y_fm * out_x,
        relative_separation_lcms_fm.x3
    };
}

/**
 * @brief Return the Euclidean radius of three spatial components.
 * @param out_fm Out separation component, in fm.
 * @param side_fm Side separation component, in fm.
 * @param long_fm Long separation component, in fm.
 * @return Spatial radius sqrt(out^2 + side^2 + long^2), in fm.
 */
double spatial_radius(
    double out_fm,
    double side_fm,
    double long_fm
) noexcept {
    const double transverse_radius_fm = std::hypot(out_fm, side_fm);
    return std::hypot(transverse_radius_fm, long_fm);
}

}  // namespace

PairFrameObservables calculate_pair_frame_observables(
    const Particle& particle_a,
    const Particle& particle_b,
    const PairKinematics& kinematics
) noexcept {
    const common::FourVector relative_separation_lab_fm{
        particle_a.position.x0 - particle_b.position.x0,
        particle_a.position.x1 - particle_b.position.x1,
        particle_a.position.x2 - particle_b.position.x2,
        particle_a.position.x3 - particle_b.position.x3
    };

    const LabToLCMSResult lcms = boost_lab_to_lcms(
        relative_separation_lab_fm,
        kinematics
    );
    const LCMSOSLProjection osl = project_to_osl(
        lcms.relative_separation_fm,
        particle_a,
        particle_b,
        kinematics
    );
    const LCMSToPRFResult prf = boost_lcms_to_prf(
        osl.delta_t_fm,
        osl.r_out_fm,
        lcms.beta_out
    );

    const double r_radial_lcms_fm = spatial_radius(
        osl.r_out_fm,
        osl.r_side_fm,
        osl.r_long_fm
    );
    const double r_radial_prf_fm = spatial_radius(
        prf.r_out_prf_fm,
        osl.r_side_fm,
        osl.r_long_fm
    );

    return {
        relative_separation_lab_fm.x0,
        osl.delta_t_fm,
        prf.delta_t_prf_fm,
        osl.r_out_fm,
        prf.r_out_prf_fm,
        osl.r_side_fm,
        osl.r_long_fm,
        r_radial_lcms_fm,
        r_radial_prf_fm
    };
}

}  // namespace hbt
