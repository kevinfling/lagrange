/**
 * lagrange_transfer.h - Impulsive Orbital Transfers (Hohmann, Bi-elliptic, General)
 * 
 * Analytic two-impulse and three-impulse transfer solutions between Keplerian orbits.
 * Integrates with patched conics for interplanetary mission design.
 * 
 * Features:
 *   - Classical Hohmann (coplanar circular)
 *   - Bi-elliptic (3-burn, optimal for large radius ratios)
 *   - Generalized elliptic transfers (arbitrary initial/final eccentricities)
 *   - Plane change combinations (combined inclination/eccentricity maneuvers)
 *   - Delta-v and time-of-flight calculations
 *   - Export to Lambert problems for high-fidelity refinement
 */

#ifndef LAGRANGE_TRANSFER_H
#define LAGRANGE_TRANSFER_H

#include "gravity.h"
#include "patched_conic.h"
#include "body.h"
#include <math.h>

#ifdef __cplusplus
extern "C" {
#endif

/*============================================================================
 * Transfer Orbit Structure
 *===========================================================================*/

typedef struct {
    lg_orbital_elements_t initial;
    lg_orbital_elements_t transfer;    /* Intermediate orbit */
    lg_orbital_elements_t target;
    
    float dv_total;                    /* Total characteristic delta-v (km/s or m/s) */
    float dv_departure;                /* First burn magnitude */
    float dv_arrival;                  /* Second burn magnitude */
    float dv_plane_change;               /* Inclination change component (if combined) */
    
    float tof;                         /* Time of flight (seconds) */
    float phase_angle;                 /* Initial lead angle for rendezvous (rad) */
    
    int n_burns;                       /* 2 for Hohmann, 3 for bi-elliptic */
    float efficiency;                    /* Ratio to ideal (C3 or energy metric) */
} lg_transfer_solution_t;

/*============================================================================
 * Hohmann Transfer (Optimal 2-impulse for coplanar circular)
 * 
 * Most efficient 2-impulse transfer between circular coplanar orbits.
 * Transfer orbit: a = (r1 + r2)/2, e = |r2-r1|/(r2+r1)
 *===========================================================================*/

/* Classic Hohmann: circular -> circular */
static inline lg_transfer_solution_t lg_transfer_hohmann(float mu, 
                                                        float r1, 
                                                        float r2) {
    lg_transfer_solution_t sol = {0};
    sol.n_burns = 2;
    
    /* Circular velocities */
    float v1 = sqrtf(mu / r1);
    float v2 = sqrtf(mu / r2);
    
    /* Transfer orbit parameters */
    float a_xfer = (r1 + r2) * 0.5f;
    float e_xfer = fabsf(r2 - r1) / (r1 + r2);
    
    /* Velocities at periapsis and apoapsis of transfer */
    float vp = sqrtf(mu * (2.0f/r1 - 1.0f/a_xfer));  /* At r1 (peri if r1<r2) */
    float va = sqrtf(mu * (2.0f/r2 - 1.0f/a_xfer));  /* At r2 (apo if r1<r2) */
    
    /* Delta-v calculations */
    if (r2 > r1) {
        /* Ascending transfer (LEO to GEO, LEO to Mars, etc) */
        sol.dv_departure = vp - v1;        /* First burn at periapsis */
        sol.dv_arrival = v2 - va;          /* Second burn at apoapsis */
        
        sol.initial.semi_major_axis = r1;
        sol.target.semi_major_axis = r2;
        sol.transfer.semi_major_axis = a_xfer;
        sol.transfer.eccentricity = e_xfer;
        sol.transfer.arg_of_periapsis = 0.0f;  /* Departure at periapsis */
    } else {
        /* Descending transfer */
        sol.dv_departure = v1 - va;        /* Burn at "apoapsis" (which is r1) */
        sol.dv_arrival = vp - v2;          /* Circularize at r2 */
        
        sol.initial.semi_major_axis = r1;
        sol.target.semi_major_axis = r2;
        sol.transfer.semi_major_axis = a_xfer;
        sol.transfer.eccentricity = e_xfer;
        sol.transfer.arg_of_periapsis = LG_PI;  /* Departure at apoapsis */
    }
    
    sol.dv_total = fabsf(sol.dv_departure) + fabsf(sol.dv_arrival);
    
    /* Time of flight: half period of transfer ellipse */
    sol.tof = LG_PI * sqrtf(a_xfer * a_xfer * a_xfer / mu);
    
    /* Phase angle for rendezvous: target must be pi*(1 - sqrt((1+e)/8)) behind */
    /* For circular Hohmann: lead angle = pi * (1 - 1/sqrt(8)) ≈ 1.72 rad ≈ 98.6 deg */
    sol.phase_angle = LG_PI * (1.0f - sqrtf(powf(2.0f*r1/(r1+r2), 3.0f)));
    
    sol.efficiency = 1.0f;  /* By definition, optimal for 2-impulse circular */
    
    return sol;
}

/* Hohmann with plane change (combined maneuver at node) */
static inline lg_transfer_solution_t lg_transfer_hohmann_plane_change(float mu,
                                                                     float r1,
                                                                     float r2,
                                                                     float delta_i,  /* Inclination change (rad) */
                                                                     bool combine_at_departure) {
    lg_transfer_solution_t sol = lg_transfer_hohmann(mu, r1, r2);
    
    /* Law of cosines for combined maneuver: dv = sqrt(v1^2 + v2^2 - 2*v1*v2*cos(di)) */
    float v_depart = sqrtf(mu/r1);
    float di = delta_i;
    
    if (combine_at_departure) {
        /* Combined burn at departure (more efficient if r1 < r2 usually) */
        float v_xfer_p = sqrtf(mu * (2.0f/r1 - 2.0f/(r1+r2)));
        sol.dv_departure = sqrtf(v_depart*v_depart + v_xfer_p*v_xfer_p - 
                                2.0f*v_depart*v_xfer_p*cosf(di));
        sol.dv_plane_change = v_depart * di;  /* Approximation for small angles */
    } else {
        /* Plane change at apoapsis (where velocity is lower, more efficient for large di) */
        float v_xfer_a = sqrtf(mu * (2.0f/r2 - 2.0f/(r1+r2)));
        float v_target = sqrtf(mu/r2);
        
        /* If ascending, apoapsis is at r2 */
        sol.dv_arrival = sqrtf(v_xfer_a*v_xfer_a + v_target*v_target - 
                              2.0f*v_xfer_a*v_target*cosf(di));
        sol.dv_plane_change = v_xfer_a * di;  /* Lower velocity = cheaper plane change */
    }
    
    sol.dv_total = sol.dv_departure + sol.dv_arrival;
    return sol;
}

/*============================================================================
 * Bi-elliptic Transfer (3-impulse, optimal for r2/r1 > ~11.94)
 * 
 * Strategy: 
 *   1. Burn to high apoapsis rb (> r2)
 *   2. Circularize (or adjust) at rb 
 *   3. Burn to drop to r2
 *   
 * Surprisingly, can be more efficient than Hohmann when r2/r1 > 11.94
 * because burn at apoapsis (very low velocity) is cheap.
 *===========================================================================*/

static inline lg_transfer_solution_t lg_transfer_bielliptic(float mu,
                                                            float r1,
                                                            float r2,
                                                            float rb) {    /* Intermediate apoapsis, must be > max(r1,r2) */
    lg_transfer_solution_t sol = {0};
    sol.n_burns = 3;
    
    if (rb <= fmaxf(r1, r2)) {
        sol.dv_total = 1e20f; /* Invalid */
        return sol;
    }
    
    float a_xfer1 = (r1 + rb) / 2.0f;
    float a_xfer2 = (rb + r2) / 2.0f;
    
    float v1_circ = sqrtf(mu/r1);
    float v2_circ = sqrtf(mu/r2);
    
    /* Velocities on transfer legs */
    float v1_xfer1 = sqrtf(mu * (2.0f/r1 - 1.0f/a_xfer1));
    float vb_xfer1 = sqrtf(mu * (2.0f/rb - 1.0f/a_xfer1));
    float vb_xfer2 = sqrtf(mu * (2.0f/rb - 1.0f/a_xfer2));
    float v2_xfer2 = sqrtf(mu * (2.0f/r2 - 1.0f/a_xfer2));
    
    /* Three burns */
    float dv1 = fabsf(v1_xfer1 - v1_circ);
    float dv2 = fabsf(vb_xfer2 - vb_xfer1);  /* Change from leg1 apoapsis velocity to leg2 apoapsis velocity */
    float dv3 = fabsf(v2_circ - v2_xfer2);
    
    sol.dv_departure = dv1;
    sol.dv_arrival = dv3;
    sol.dv_total = dv1 + dv2 + dv3;
    
    /* TOF: half period of first ellipse + half period of second ellipse */
    sol.tof = LG_PI * (sqrtf(a_xfer1*a_xfer1*a_xfer1/mu) + sqrtf(a_xfer2*a_xfer2*a_xfer2/mu));
    
    return sol;
}

/* Find optimal rb for bi-elliptic given r1, r2 */
static inline float lg_transfer_bielliptic_optimal_rb(float r1, float r2) {
    /* As rb -> infinity, bi-elliptic becomes parabolic transfer */
    /* Optimal is usually very high, but practically limited by SOI/planetary constraints */
    /* Return a reasonable multiple, e.g., 10 * max(r1, r2) or SOI radius */
    return 10.0f * fmaxf(r1, r2);
}

/* Check if bi-elliptic beats Hohmann */
static inline bool lg_transfer_bielliptic_wins(float r1, float r2, float rb_limit) {
    float ratio = fmaxf(r1, r2) / fminf(r1, r2);
    /* Critical ratio is ~11.94 for rb -> infinity */
    /* With finite rb, the threshold is slightly higher */
    return ratio > 11.94f && rb_limit > 5.0f * fmaxf(r1, r2);
}

/*============================================================================
 * Generalized Coplanar Transfer (arbitrary initial/final eccentricity)
 * Uses Lambert's problem from patched_conic.h
 *===========================================================================*/

static inline lg_transfer_solution_t lg_transfer_generalized(float mu,
                                                            const lg_orbital_elements_t* initial,
                                                            const lg_orbital_elements_t* target,
                                                            float tof_hint,  /* Guess for TOF, or 0 for auto */
                                                            bool short_way) {
    lg_transfer_solution_t sol = {0};
    
    /* Convert elements to state vectors at specific true anomalies */
    /* For minimum energy transfer, we need to solve Lambert's problem */
    
    lg_vec3_t r1, v1, r2, v2;
    
    /* Current position at some anomaly (0 for simplicity, or optimized) */
    lg_elements_to_state(initial, mu, &r1, &v1);
    lg_elements_to_state(target, mu, &r2, &v2);
    
    /* If tof_hint is 0, estimate from Hohmann */
    if (tof_hint <= 0.0f) {
        float a1 = initial->semi_major_axis;
        float a2 = target->semi_major_axis;
        float a_xfer = (a1 + a2) * 0.5f;
        tof_hint = LG_PI * sqrtf(a_xfer*a_xfer*a_xfer/mu);
    }
    
    lg_lambert_problem_t prob = {
        .r1 = r1, .r2 = r2, .dt = tof_hint, .mu = mu,
        .short_way = short_way, .max_iter = 50, .tol = 1e-6f
    };
    
    lg_lambert_solution_t lambert = lg_lambert_solve(&prob);
    
    if (lambert.status == 0) {
        sol.dv_departure = lg_vec3_len(lg_vec3_sub(lambert.v1, v1));
        sol.dv_arrival = lg_vec3_len(lg_vec3_sub(v2, lambert.v2));
        sol.dv_total = sol.dv_departure + sol.dv_arrival;
        sol.tof = tof_hint;
        sol.transfer.semi_major_axis = lambert.a;
        sol.transfer.eccentricity = lambert.ecc;
    } else {
        sol.dv_total = 1e20f; /* Failed */
    }
    
    return sol;
}

/*============================================================================
 * Phasing/Rendezvous (matching position after transfer)
 *===========================================================================*/

/* Calculate wait time before departure for optimal phase alignment */
static inline float lg_transfer_phase_wait(float mu, float r1, float r2, 
                                          float target_mean_motion,  /* Of target orbit */
                                          bool target_ahead) {
    /* Phase angle calculated in Hohmann */
    lg_transfer_solution_t hohm = lg_transfer_hohmann(mu, r1, r2);
    float lead_angle = hohm.phase_angle;  /* How far ahead target should be */
    
    /* Current angular separation */
    float n_target = target_mean_motion;
    float n_transfer = LG_PI / hohm.tof;  /* Mean motion of transfer (half period = pi) */
    
    /* Synodic period */
    float synodic = 2.0f * LG_PI / fabsf(n_transfer - n_target);
    
    /* Wait time to achieve lead_angle separation */
    float wait = lead_angle / fabsf(n_target - n_transfer);
    if (!target_ahead) wait = synodic - wait;
    
    float result = fmodf(wait, synodic);
    return result > 0 ? result : result + synodic;
}

/*============================================================================
 * Batch Mission Planning (Porkchop for simple transfers)
 *===========================================================================*/

typedef struct {
    float r1, r2;              /* Initial and final radii (or use ephemeris) */
    float* departure_times;    /* Array of departure epochs */
    float* tofs;               /* Array of times of flight */
    float* dv_grid;            /* Output: total dv for each combination */
    int n_deps;
    int n_tofs;
} lg_transfer_grid_t;

/* Fill a grid with Hohmann delta-v (analytic, fast) */
static inline void lg_transfer_grid_hohmann(const lg_transfer_grid_t* grid, float mu) {
    for (int i = 0; i < grid->n_deps; i++) {
        for (int j = 0; j < grid->n_tofs; j++) {
            /* For Hohmann, TOF is fixed by geometry, but we can calculate dv for given r1,r2 */
            /* Actually for grid we usually vary departure epoch and TOF independently */
            /* Here we just compute the theoretical Hohmann dv as baseline */
            lg_transfer_solution_t sol = lg_transfer_hohmann(mu, grid->r1, grid->r2);
            grid->dv_grid[i * grid->n_tofs + j] = sol.dv_total;
        }
    }
}

/*============================================================================
 * Export to Particle System for Simulation
 *===========================================================================*/

/* Initialize spacecraft and targets in particle system for transfer verification */
static inline void lg_transfer_to_particle_system(const lg_transfer_solution_t* sol,
                                                 const lg_body_t* central_body,
                                                 float mu,
                                                 lg_particle_system_t* ps,
                                                 int* sc_idx) {
    /* Add spacecraft as massless particle (or small mass) */
    int idx = ps->n;  /* Append */
    ps->n++;
    
    /* Initial state at departure */
    lg_elements_to_state(&sol->initial, mu, 
                        &(lg_vec3_t){ps->pos.x[idx], ps->pos.y[idx], ps->pos.z[idx]},
                        &(lg_vec3_t){ps->vel.x[idx], ps->vel.y[idx], ps->vel.z[idx]});
    
    /* Apply departure burn */
    /* ... modify velocity by sol->dv_departure in prograde direction ... */
    
    *sc_idx = idx;
}

/*============================================================================
 * Utility: Interplanetary Transfer Calculator (Earth -> Mars example)
 *===========================================================================*/

/* Quick calculation for Earth-Mars Hohmann window */
static inline lg_transfer_solution_t lg_transfer_earth_mars_hohmann(void) {
    /* Earth: 1 AU, Mars: 1.524 AU */
    float r_earth = 1.496e11f;  /* meters */
    float r_mars = 2.279e11f;
    float mu_sun = 1.32712440018e20f;
    
    return lg_transfer_hohmann(mu_sun, r_earth, r_mars);
}

/* Bi-elliptic via Jupiter (gravity assist style, but as pure 3-burn) */
static inline lg_transfer_solution_t lg_transfer_earth_jupiter_mars_bielliptic(void) {
    float r_e = 1.496e11f;
    float r_j = 7.785e11f;  /* Jupiter semi-major axis */
    float r_m = 2.279e11f;
    float mu_sun = 1.32712440018e20f;
    
    /* Use Jupiter orbit as intermediate rb */
    return lg_transfer_bielliptic(mu_sun, r_e, r_m, r_j);
}

#ifdef __cplusplus
}
#endif

#endif /* LAGRANGE_TRANSFER_H */

