/**
 * lagrange_lambert.h - Lambert Problem Solver
 * 
 * Universal variable solution to the boundary value problem of orbital mechanics.
 * Given two position vectors and time of flight, find the transfer orbit.
 * 
 * This header provides unified type definitions to avoid conflicts between
 * patched_conic.h and propulsion.h.
 */

#ifndef LAGRANGE_LAMBERT_H
#define LAGRANGE_LAMBERT_H

#include "math.h"
#include <stdbool.h>
#include <math.h>

#ifdef __cplusplus
extern "C" {
#endif

/*============================================================================
 * Lambert Problem Types (Unified)
 *===========================================================================*/

typedef struct {
    lg_vec3_t r1, r2;         /* Initial and final position vectors (km or AU) */
    float dt;                 /* Time of flight (seconds) - unified field name */
    float mu;                 /* Gravitational parameter (km^3/s^2 or AU^3/day^2) */
    bool short_way;           /* True for < 180 deg transfer, false for > 180 */
    int max_iter;             /* Default 50 */
    float tol;                /* Convergence tolerance, default 1e-6 */
} lg_lambert_problem_t;

typedef struct {
    lg_vec3_t v1;             /* Departure velocity */
    lg_vec3_t v2;             /* Arrival velocity */
    float a;                  /* Semi-major axis of transfer (-ve for hyperbolic) */
    float p;                  /* Semi-latus rectum */
    float ecc;                /* Eccentricity */
    float nu1, nu2;           /* True anomalies at r1, r2 */
    int iter;                 /* Iteration count */
    int status;               /* 0 = success, 1 = fail, 2 = no solution */
    float x;                  /* Universal variable solution (x = sqrt(a) * dE for ellipse) */
} lg_lambert_solution_t;

/*============================================================================
 * Stumpff Functions (Universal Variable)
 *===========================================================================*/

static inline float lg_stumpff_c2(float psi) {
    if (psi > 1e-6f) return (1.0f - cosf(sqrtf(psi))) / psi;
    if (psi < -1e-6f) return (1.0f - coshf(sqrtf(-psi))) / psi;
    return 0.5f - psi/24.0f + psi*psi/720.0f; /* Taylor expansion */
}

static inline float lg_stumpff_c3(float psi) {
    if (psi > 1e-6f) {
        float s = sqrtf(psi);
        return (s - sinf(s)) / (psi * s);
    }
    if (psi < -1e-6f) {
        float s = sqrtf(-psi);
        return (sinhf(s) - s) / ((-psi) * s);
    }
    return 1.0f/6.0f - psi/120.0f + psi*psi/5040.0f;
}

/*============================================================================
 * Lambert Solver (Izzo's Algorithm - Robust for all conic sections)
 *===========================================================================*/

static inline lg_lambert_solution_t lg_lambert_solve(const lg_lambert_problem_t* prob) {
    lg_lambert_solution_t sol = {0};
    sol.status = 1;
    
    float r1 = lg_vec3_len(prob->r1);
    float r2 = lg_vec3_len(prob->r2);
    float r1_dot_r2 = lg_vec3_dot(prob->r1, prob->r2);
    float c = sqrtf(r1*r1 + r2*r2 - 2.0f*r1_dot_r2); /* Chord */
    float s = (r1 + r2 + c) * 0.5f; /* Semi-perimeter */
    
    float ir1 = 1.0f/r1, ir2 = 1.0f/r2;
    lg_vec3_t r1_hat = lg_vec3_scale(prob->r1, ir1);
    lg_vec3_t r2_hat = lg_vec3_scale(prob->r2, ir2);
    
    /* Cross product for plane normal */
    lg_vec3_t h_hat = lg_vec3_norm(lg_vec3_cross(r1_hat, r2_hat));
    float dnu = atan2f(lg_vec3_len(lg_vec3_cross(r1_hat, r2_hat)), r1_dot_r2);
    if (!prob->short_way) dnu = 2.0f*M_PI - dnu;
    
    /* Initial guess for x (universal variable) */
    float A = sqrtf(r1*r2) * sinf(dnu) / sqrtf(1.0f - cosf(dnu)); /* Lambert parameter */
    if (!prob->short_way) A = -A;
    
    float x = 0.0f; /* Initial guess */
    if (prob->dt * prob->dt < 4.0f * M_PI * M_PI * s * s * s / (27.0f * prob->mu)) {
        /* Elliptic guess */
        float alpha = M_PI;
        float beta = 2.0f * asinf(sqrtf((s-c)/s));
        x = (alpha - beta) / sqrtf(s);
    }
    
    /* Newton iteration on time equation */
    for (int i = 0; i < prob->max_iter; i++) {
        float psi = x * x;
        float c2 = lg_stumpff_c2(psi);
        float c3 = lg_stumpff_c3(psi);
        
        float y = r1 + r2 + A * (psi * c3 - 1.0f) / sqrtf(c2);
        if (y < 0 || A > 0 && y < 0) { /* Reject invalid region */
            x *= 0.5f;
            continue;
        }
        
        float chi = sqrtf(y / c2);
        float dt_calc = (chi*chi*chi * c3 + A * sqrtf(y)) / sqrtf(prob->mu);
        
        float dtdx = (chi*chi*chi * (4.0f - psi * c2/c3) + A * (6.0f*y*psi + psi - 4.0f*y)) / 
                     (4.0f * y * chi);
        
        float err = dt_calc - prob->dt;
        if (fabsf(err) < prob->tol) {
            /* Converged - compute velocities */
            float f = 1.0f - y/r1;
            float g = A * sqrtf(y / prob->mu);
            float g_dot = 1.0f - y/r2;
            
            lg_vec3_t r1_hat_local = lg_vec3_norm(prob->r1);
            lg_vec3_t r2_hat_local = lg_vec3_norm(prob->r2);
            
            /* Velocity at r1 */
            lg_vec3_t v1_temp = lg_vec3_sub(prob->r2, lg_vec3_scale(prob->r1, f));
            sol.v1 = lg_vec3_scale(v1_temp, 1.0f/g);
            
            /* Velocity at r2 */
            lg_vec3_t v2_temp = lg_vec3_sub(lg_vec3_scale(prob->r2, g_dot), prob->r1);
            sol.v2 = lg_vec3_scale(v2_temp, 1.0f/g);
            
            sol.a = y / (2.0f - psi);
            sol.p = y;
            sol.ecc = sqrtf(1.0f - (A*A / (sol.a * (2.0f - psi))));
            sol.iter = i;
            sol.status = 0;
            sol.x = x;
            return sol;
        }
        
        x -= err / dtdx;
    }
    
    return sol;
}

#ifdef __cplusplus
}
#endif

#endif /* LAGRANGE_LAMBERT_H */
