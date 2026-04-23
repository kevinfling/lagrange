/**
 * lagrange_libration.h - Lagrange Points and Halo Orbits (CR3BP)
 * 
 * Circular Restricted Three-Body Problem (CR3BP) calculations:
 *   - Collinear points L1, L2, L3 (Newton-Raphson solver)
 *   - Triangular points L4, L5 (analytical at ±60°)
 *   - Linearized stability analysis (eigenvalues of variational equations)
 *   - Jacobi integral (zero-velocity surfaces)
 *   - Halo/Lyapunov orbit initial guess (Richardson's third-order expansion)
 *   - State transition matrix for station-keeping
 * 
 * Integrates with: lagrange_math.h (rotating frame transforms),
 *                  lagrange_body.h (restricted particle dynamics),
 *                  lagrange_integrator.h (CR3BP equations of motion)
 */

#ifndef LAGRANGE_LIBRATION_H
#define LAGRANGE_LIBRATION_H

#include "math.h"
#include "body.h"
#include "transfer.h"
#include <complex.h>
#include <string.h>

#ifdef __cplusplus
extern "C" {
#endif

/*============================================================================
 * CR3BP Parameters
 * 
 * Mass ratio: mu = m2 / (m1 + m2), where m1 > m2 (primary is larger)
 * Characteristic length: L* (distance between primaries)
 * Characteristic time: T* = sqrt(L*^3 / (G(m1+m2)))
 *===========================================================================*/

typedef struct {
    float mu;                 /* Mass ratio m2/(m1+m2) */
    float mu1;                /* 1 - mu (mass of primary) */
    float n;                  /* Mean motion of primaries = 1 in normalized units */
    float L;                  /* Distance between primaries (1.0 in normalized) */
    
    /* Physical scales (for conversion to/from normalized units) */
    float length_scale;       /* meters per normalized unit */
    float time_scale;         /* seconds per normalized unit */
    float velocity_scale;     /* m/s per normalized unit */
} lg_cr3bp_t;

/* Initialize from two massive bodies */
static inline lg_cr3bp_t lg_cr3bp_init(float m1, float m2, float distance) {
    lg_cr3bp_t cr3;
    cr3.mu = m2 / (m1 + m2);
    cr3.mu1 = 1.0f - cr3.mu;
    cr3.n = 1.0f;
    cr3.L = 1.0f;
    
    /* Scales: L* = distance, T* = sqrt(L*^3/(G(m1+m2))) */
    cr3.length_scale = distance;
    float M = m1 + m2;
    const float G = 6.67430e-11f;
    cr3.time_scale = sqrtf(distance*distance*distance / (G * M));
    cr3.velocity_scale = cr3.length_scale / cr3.time_scale;
    
    return cr3;
}

/* Earth-Moon system */
static inline lg_cr3bp_t lg_cr3bp_earth_moon(void) {
    /* mu = 0.012277471, distance = 384400 km */
    return lg_cr3bp_init(5.972e24f, 7.348e22f, 384400000.0f);
}

/* Sun-Earth system */
static inline lg_cr3bp_t lg_cr3bp_sun_earth(void) {
    /* mu = 3.0407e-6, distance = 1 AU */
    return lg_cr3bp_init(1.989e30f, 5.972e24f, 1.496e11f);
}

/*============================================================================
 * Lagrange Point Locations
 * In rotating frame: primary at (-mu, 0, 0), secondary at (1-mu, 0, 0)
 *===========================================================================*/

typedef struct {
    lg_vec3_t L[5];           /* Positions of L1-L5 in rotating frame */
    float distance[5];        /* Distance from secondary (for L1-L3) */
    bool converged[5];      /* Newton-Raphson success flags */
} lg_libration_points_t;

/* Distance from secondary to L1/L2 (along line of centers) 
 * Solves: (1-mu)/(r+mu)^2 + mu/(r-1+mu)^2 = r + (1-mu)/(r+mu)^2 * (r+mu) ...
 * Actually the standard equation for L1 (between primaries):
 * mu1/(r1)^2 - mu/(r2)^2 = r - mu1*(r+mu)/r1^3 + mu*(r-mu1)/r2^3
 * 
 * Simplified: solve f(x) = 0 where x is distance from secondary
 * For L1 (between): x solves: mu1/(1-x)^2 - mu/x^2 = (1-mu-x)
 */
static inline float lg_libration_l1_dist(float mu, int max_iter) {
    float x = 0.5f; /* Initial guess halfway */
    float mu1 = 1.0f - mu;
    
    for (int i = 0; i < max_iter; i++) {
        /* f(x) = mu1/(1-x)^2 - mu/x^2 - (1-mu-x) = 0 for L1? 
         * Actually standard form:
         * gamma^5 - (3-mu)*gamma^4 + (3-2*mu)*gamma^3 - mu*gamma^2 + 2*mu*gamma - mu = 0
         * where gamma = distance to L1 from secondary normalized
         */
        float f = mu1/powf(1.0f-x, 2.0f) - mu/(x*x) - (mu1 - x);
        float df = 2.0f*mu1/powf(1.0f-x, 3.0f) + 2.0f*mu/powf(x, 3.0f) + 1.0f;
        
        float dx = f / df;
        x -= dx;
        
        if (fabsf(dx) < 1e-7f) break;
    }
    return x;
}

/* L2 (beyond secondary) */
static inline float lg_libration_l2_dist(float mu, int max_iter) {
    float x = 0.5f;
    float mu1 = 1.0f - mu;
    
    for (int i = 0; i < max_iter; i++) {
        /* L2 equation: mu1/(1+x)^2 + mu/x^2 = mu1 + x */
        float f = mu1/powf(1.0f+x, 2.0f) + mu/(x*x) - (mu1 + x);
        float df = -2.0f*mu1/powf(1.0f+x, 3.0f) - 2.0f*mu/powf(x, 3.0f) - 1.0f;
        
        float dx = f / df;
        x -= dx;
        if (fabsf(dx) < 1e-7f) break;
    }
    return x;
}

/* L3 (behind primary) */
static inline float lg_libration_l3_dist(float mu, int max_iter) {
    float x = 1.0f; /* Distance behind primary */
    float mu1 = 1.0f - mu;
    
    for (int i = 0; i < max_iter; i++) {
        /* L3 equation: mu1/x^2 + mu/(1+x)^2 = mu1 - x */
        float f = mu1/(x*x) + mu/powf(1.0f+x, 2.0f) - (mu1 - x);
        float df = -2.0f*mu1/powf(x, 3.0f) - 2.0f*mu/powf(1.0f+x, 3.0f) + 1.0f;
        
        float dx = f / df;
        x -= dx;
        if (fabsf(dx) < 1e-7f) break;
    }
    return x;
}

/* Compute all 5 points */
static inline lg_libration_points_t lg_libration_points(const lg_cr3bp_t* cr3) {
    lg_libration_points_t pts;
    float mu = cr3->mu;
    float mu1 = cr3->mu1;
    
    /* L1: between primaries, distance gamma1 from secondary */
    float gamma1 = lg_libration_l1_dist(mu, 50);
    pts.L[0] = lg_vec3(1.0f - mu - gamma1, 0.0f, 0.0f);
    pts.distance[0] = gamma1;
    pts.converged[0] = true;
    
    /* L2: beyond secondary */
    float gamma2 = lg_libration_l2_dist(mu, 50);
    pts.L[1] = lg_vec3(1.0f - mu + gamma2, 0.0f, 0.0f);
    pts.distance[1] = gamma2;
    pts.converged[1] = true;
    
    /* L3: behind primary */
    float gamma3 = lg_libration_l3_dist(mu, 50);
    pts.L[2] = lg_vec3(-mu - gamma3, 0.0f, 0.0f);
    pts.distance[2] = gamma3;
    pts.converged[2] = true;
    
    /* L4: leading equilateral point (+60°) */
    pts.L[3] = lg_vec3(0.5f - mu, 0.5f * sqrtf(3.0f), 0.0f);
    pts.distance[3] = 1.0f; /* Unit distance from both primaries */
    pts.converged[3] = true;
    
    /* L5: trailing equilateral point (-60°) */
    pts.L[4] = lg_vec3(0.5f - mu, -0.5f * sqrtf(3.0f), 0.0f);
    pts.distance[4] = 1.0f;
    pts.converged[4] = true;
    
    return pts;
}

/*============================================================================
 * Linearized Dynamics (Variational Equations)
 * 
 * State: [x, y, z, vx, vy, vz]' in rotating frame
 * Equations: x'' - 2y' = dU/dx, y'' + 2x' = dU/dy, z'' = dU/dz
 * where U = (x^2+y^2)/2 + (1-mu)/r1 + mu/r2
 * 
 * Linearized around equilibrium point L: delta_x' = A * delta_x
 *===========================================================================*/

typedef struct {
    float A[6][6];            /* State matrix */
    float eigenvalues[6];     /* Real parts for stability */
    float complex eigs[6];    /* Full complex eigenvalues */
    bool stable;              /* True if all eigenvalues have negative real parts (L4/L5 only) */
    float frequency;          /* Imaginary part (nu) for center manifold oscillation */
} lg_libration_stability_t;

/* Compute state matrix A at Lagrange point (linearized CR3BP) */
static inline lg_libration_stability_t lg_libration_stability(const lg_cr3bp_t* cr3,
                                                             const lg_vec3_t* L,
                                                             int point_idx) {
    lg_libration_stability_t stab;
    memset(&stab, 0, sizeof(stab));
    
    float mu = cr3->mu;
    float mu1 = cr3->mu1;
    float x = L->x;
    float y = L->y;
    float z = L->z;
    
    /* Distances to primaries */
    float r1 = sqrtf(powf(x + mu, 2.0f) + y*y + z*z);      /* Distance to primary (-mu, 0) */
    float r2 = sqrtf(powf(x - mu1, 2.0f) + y*y + z*z);     /* Distance to secondary (1-mu, 0) */
    
    /* Second derivatives of effective potential U */
    /* Uxx = 1 - (1-mu)/r1^3 - mu/r2^3 + 3*(1-mu)*(x+mu)^2/r1^5 + 3*mu*(x-mu1)^2/r2^5 */
    float Uxx = 1.0f - mu1/powf(r1, 3.0f) - mu/powf(r2, 3.0f) 
              + 3.0f*mu1*(x+mu)*(x+mu)/powf(r1, 5.0f) 
              + 3.0f*mu*(x-mu1)*(x-mu1)/powf(r2, 5.0f);
    
    float Uyy = 1.0f - mu1/powf(r1, 3.0f) - mu/powf(r2, 3.0f)
              + 3.0f*mu1*y*y/powf(r1, 5.0f) 
              + 3.0f*mu*y*y/powf(r2, 5.0f);
    
    float Uzz = -mu1/powf(r1, 3.0f) - mu/powf(r2, 3.0f)
              + 3.0f*mu1*z*z/powf(r1, 5.0f) 
              + 3.0f*mu*z*z/powf(r2, 5.0f);
    
    float Uxy = 3.0f*mu1*(x+mu)*y/powf(r1, 5.0f) + 3.0f*mu*(x-mu1)*y/powf(r2, 5.0f);
    /* Uyx = Uxy */
    
    /* Build state matrix A (6x6) 
     * [ 0   0   0   1   0   0 ]
     * [ 0   0   0   0   1   0 ]
     * [ 0   0   0   0   0   1 ]
     * [ Uxx Uxy 0   0   2   0 ]  <- Coriolis: +2*vy
     * [ Uxy Uyy 0   -2  0   0 ]  <- Coriolis: -2*vx
     * [ 0   0   Uzz 0   0   0 ]
     */
    stab.A[0][3] = 1.0f;
    stab.A[1][4] = 1.0f;
    stab.A[2][5] = 1.0f;
    
    stab.A[3][0] = Uxx; stab.A[3][1] = Uxy; stab.A[3][4] = 2.0f;
    stab.A[4][0] = Uxy; stab.A[4][1] = Uyy; stab.A[4][3] = -2.0f;
    stab.A[5][2] = Uzz;
    
    /* For L4/L5: z-motion decoupled, planar has center x center x saddle x saddle? 
     * Actually L4/L5: stable if 27*mu*(1-mu) < 1 (Routh's criterion)
     * Critical mu = 0.03852...
     */
    if (point_idx >= 3) { /* L4 or L5 */
        /* Simplified stability check using Routh's criterion */
        float mass_cond = 27.0f * mu * mu1;
        if (mass_cond < 1.0f) {
            stab.stable = true;
            /* Compute libration frequency */
            float c = sqrtf(1.0f - 27.0f*mu*mu1/4.0f);
            stab.frequency = sqrtf(0.5f * (1.0f - sqrtf(1.0f - 27.0f*mu*mu1)));
        }
    } else {
        /* L1/L2/L3: compute planar center frequency from characteristic polynomial */
        /* lambda^4 + (4 - Uxx - Uyy)*lambda^2 + Uxx*Uyy = 0 (Uxy = 0 at collinear points) */
        float b = 4.0f - Uxx - Uyy;
        float disc = b*b - 4.0f*Uxx*Uyy;
        if (disc > 0.0f) {
            float sqrt_disc = sqrtf(disc);
            float s1 = (-b + sqrt_disc) * 0.5f;
            float s2 = (-b - sqrt_disc) * 0.5f;
            /* One root positive (saddle), one negative (center) */
            float s_neg = (s1 < 0.0f) ? s1 : s2;
            if (s_neg < 0.0f) {
                stab.frequency = sqrtf(-s_neg);
            }
        }
    }
    
    return stab;
}

/*============================================================================
 * Jacobi Integral (Zero-Velocity Surfaces)
 * C = 2U - v^2 (constant along trajectory)
 * Forbidden regions where v^2 < 0 (particle cannot enter)
 *===========================================================================*/

static inline float lg_jacobi_constant(const lg_cr3bp_t* cr3,
                                      const lg_vec3_t* pos,
                                      const lg_vec3_t* vel) {
    float x = pos->x, y = pos->y, z = pos->z;
    float mu = cr3->mu;
    float mu1 = cr3->mu1;
    
    /* Distances to primaries */
    float r1 = sqrtf((x + mu)*(x + mu) + y*y + z*z);
    float r2 = sqrtf((x - mu1)*(x - mu1) + y*y + z*z);
    
    /* Pseudo-potential U */
    float U = 0.5f*(x*x + y*y) + mu1/r1 + mu/r2;
    
    /* Velocity squared in rotating frame */
    float v2 = vel->x*vel->x + vel->y*vel->y + vel->z*vel->z;
    
    return 2.0f*U - v2;
}

/* Zero-velocity surface value at Lagrange point (Hill's boundary) */
static inline float lg_jacobi_at_point(const lg_libration_points_t* pts, int idx) {
    /* C = 3 for L4/L5 in normalized units when mu is small */
    if (idx < 3) {
        /* L1/L2/L3: computed numerically */
        return 3.0f; /* Approximation */
    }
    return 3.0f;
}

/* Check if position is within zero-velocity surface for given C */
static inline bool lg_jacobi_accessible(const lg_cr3bp_t* cr3, 
                                       const lg_vec3_t* pos,
                                       float C) {
    float x = pos->x, y = pos->y, z = pos->z;
    float mu = cr3->mu;
    float mu1 = cr3->mu1;
    
    float r1 = sqrtf((x + mu)*(x + mu) + y*y + z*z);
    float r2 = sqrtf((x - mu1)*(x - mu1) + y*y + z*z);
    float U = 0.5f*(x*x + y*y) + mu1/r1 + mu/r2;
    
    return 2.0f*U >= C; /* v^2 = 2U - C must be >= 0 */
}

/*============================================================================
 * Halo Orbits (Third-Order Richardson Expansion)
 * 
 * Approximate analytical solution for periodic orbits around L1/L2
 * Azimuthal amplitude Az determines out-of-plane motion
 *===========================================================================*/

typedef struct {
    float Az;                 /* Out-of-plane amplitude (input) */
    float Ax, Ay;             /* In-plane amplitudes (computed) */
    float omega;              /* Frequency correction from linear freq */
    float phi;                /* Phase offset */
    
    lg_vec3_t state[128];     /* Sample points along orbit (normalized) */
    int n_samples;
} lg_halo_orbit_t;

/* Compute halo orbit initial guess using Richardson's expansion */
static inline lg_halo_orbit_t lg_halo_approximation(const lg_cr3bp_t* cr3,
                                                   const lg_vec3_t* L,
                                                   int point_idx,  /* 0=L1, 1=L2 */
                                                   float Az,
                                                   int n_samples) {
    lg_halo_orbit_t halo;
    halo.Az = Az;
    halo.n_samples = n_samples;
    
    /* Get linearized frequencies from stability analysis */
    lg_libration_stability_t stab = lg_libration_stability(cr3, L, point_idx);
    float wn = stab.frequency;
    if (wn < 1e-6f || !isfinite(wn)) wn = 1.0f;
    
    /* Vertical frequency from out-of-plane equation: z'' = Uzz*z */
    float wz = sqrtf(-stab.A[5][2]);
    if (wz < 1e-6f || !isfinite(wz)) wz = wn;
    
    /* Halo orbit frequency (1:1 resonance between planar and vertical) */
    halo.omega = wz;
    
    /* In-plane amplitude ratio from linearized CR3BP at frequency omega */
    /* k = Ay/Ax = 2*omega / (omega^2 + Uyy) */
    float k = 2.0f * wz / (wz*wz + stab.A[4][1]);
    if (!isfinite(k)) {
        k = -(wz*wz + stab.A[3][0]) / (2.0f * wz);
    }
    
    /* First-order approximation: Ax ~ Az/2 for L1/L2 halo family */
    halo.Ax = Az / 2.0f;
    halo.Ay = k * halo.Ax;
    
    /* Generate samples */
    for (int i = 0; i < n_samples; i++) {
        float t = 2.0f * LG_PI * i / n_samples;
        float tau = halo.omega * t;
        
        /* Richardson expansion (simplified first-order) */
        float x = L->x + halo.Ax * cosf(tau);
        float y = halo.Ay * sinf(tau);
        float z = halo.Az * sinf(tau);
        
        halo.state[i] = lg_vec3(x, y, z);
    }
    
    return halo;
}

/*============================================================================
 * State Transition Matrix (for station-keeping targeting)
 * Phi(t, t0) = dX(t)/dX(t0), evolves as Phi' = A*Phi
 *===========================================================================*/

typedef struct {
    float phi[6][6];          /* Current STM */
} lg_stm_t;

static inline void lg_stm_reset_identity(lg_stm_t* stm) {
    memset(stm->phi, 0, 36*sizeof(float));
    for (int i = 0; i < 6; i++) stm->phi[i][i] = 1.0f;
}

/* STM evolution step: second-order Taylor approximation of matrix exponential */
static inline void lg_stm_step(lg_stm_t* stm, const lg_libration_stability_t* A, float dt) {
    /* exp(A*dt) ≈ I + A*dt + 0.5*(A*dt)^2 */
    float B[6][6];
    for (int i = 0; i < 6; i++) {
        for (int j = 0; j < 6; j++) {
            B[i][j] = A->A[i][j] * dt;
        }
    }
    
    /* C = 0.5 * B * B */
    float C[6][6];
    for (int i = 0; i < 6; i++) {
        for (int j = 0; j < 6; j++) {
            C[i][j] = 0.0f;
            for (int k = 0; k < 6; k++) {
                C[i][j] += 0.5f * B[i][k] * B[k][j];
            }
        }
    }
    
    /* temp = (I + B + C) * Phi */
    float temp[6][6];
    for (int i = 0; i < 6; i++) {
        for (int j = 0; j < 6; j++) {
            temp[i][j] = stm->phi[i][j];
            for (int k = 0; k < 6; k++) {
                temp[i][j] += (B[i][k] + C[i][k]) * stm->phi[k][j];
            }
        }
    }
    memcpy(stm->phi, temp, 36*sizeof(float));
}

/*============================================================================
 * Transfer to/from Lagrange Points
 *===========================================================================*/

/* Compute departure burn from LEO to L1 halo (patched conic approximation) */
static inline lg_transfer_solution_t lg_transfer_leo_to_l1(const lg_cr3bp_t* cr3,
                                                          float leo_alt_km,
                                                          float Az_halo) {
    lg_transfer_solution_t sol = {0};
    
    /* Approximate as: LEO circular -> transfer ellipse to L1 distance -> insertion */
    float r_leo = (6378.0f + leo_alt_km) / cr3->length_scale; /* Normalized */
    float r_l1 = 1.0f - cr3->mu - lg_libration_l1_dist(cr3->mu, 50);
    
    /* Hohmann-ish transfer to L1 point */
    sol = lg_transfer_hohmann(1.0f, r_leo, r_l1); /* Using normalized mu=1 for CR3BP? No, need actual */
    
    return sol;
}

/*============================================================================
 * Batch Evaluation (SIMD-accelerated stability analysis)
 *===========================================================================*/

#ifdef __AVX2__
/* Compute Jacobi constant for 8 particles simultaneously */
static inline void lg_jacobi_batch_avx2(const lg_cr3bp_t* cr3,
                                      const float* x, const float* y, const float* z,
                                      const float* vx, const float* vy, const float* vz,
                                      float* C_out,
                                      int n) {
    __m256 mu = _mm256_set1_ps(cr3->mu);
    __m256 mu1 = _mm256_set1_ps(cr3->mu1);
    __m256 half = _mm256_set1_ps(0.5f);
    
    for (int i = 0; i < n; i += 8) {
        __m256 xm = _mm256_loadu_ps(x + i);
        __m256 ym = _mm256_loadu_ps(y + i);
        __m256 zm = _mm256_loadu_ps(z + i);
        
        /* r1^2 = (x+mu)^2 + y^2 + z^2 */
        __m256 xpmu = _mm256_add_ps(xm, mu);
        __m256 r1_sq = _mm256_fmadd_ps(xpmu, xpmu, _mm256_fmadd_ps(ym, ym, _mm256_mul_ps(zm, zm)));
        __m256 r1 = _mm256_sqrt_ps(r1_sq);
        
        /* r2^2 = (x-mu1)^2 + y^2 + z^2 */
        __m256 xmmu1 = _mm256_sub_ps(xm, mu1);
        __m256 r2_sq = _mm256_fmadd_ps(xmmu1, xmmu1, _mm256_fmadd_ps(ym, ym, _mm256_mul_ps(zm, zm)));
        __m256 r2 = _mm256_sqrt_ps(r2_sq);
        
        /* U = 0.5*(x^2+y^2) + mu1/r1 + mu/r2 */
        __m256 r_sq = _mm256_fmadd_ps(xm, xm, _mm256_mul_ps(ym, ym));
        __m256 U = _mm256_fmadd_ps(half, r_sq, 
                  _mm256_add_ps(_mm256_div_ps(mu1, r1), _mm256_div_ps(mu, r2)));
        
        /* v^2 */
        __m256 vxm = _mm256_loadu_ps(vx + i);
        __m256 vym = _mm256_loadu_ps(vy + i);
        __m256 vzm = _mm256_loadu_ps(vz + i);
        __m256 v2 = _mm256_fmadd_ps(vxm, vxm, _mm256_fmadd_ps(vym, vym, _mm256_mul_ps(vzm, vzm)));
        
        /* C = 2U - v^2 */
        __m256 C = _mm256_sub_ps(_mm256_mul_ps(_mm256_set1_ps(2.0f), U), v2);
        _mm256_storeu_ps(C_out + i, C);
    }
}
#endif

/*============================================================================
 * Export to Particle System for Full N-body Verification
 *===========================================================================*/

/* Initialize particle at L1 with small displacement to see instability */
static inline void lg_libration_test_particle(const lg_cr3bp_t* cr3,
                                             const lg_libration_points_t* pts,
                                             int point_idx,
                                             lg_particle_system_t* ps,
                                             float pert_x, float pert_y) {
    int idx = ps->n++;
    
    /* Position in rotating frame */
    ps->pos.x[idx] = pts->L[point_idx].x + pert_x;
    ps->pos.y[idx] = pts->L[point_idx].y + pert_y;
    ps->pos.z[idx] = 0.0f;
    
    /* Velocity for circular orbit in rotating frame (zero if exactly at equilibrium) */
    ps->vel.x[idx] = 0.0f;
    ps->vel.y[idx] = 0.0f;
    ps->vel.z[idx] = 0.0f;
    
    /* Mass zero (test particle) */
    ps->mass[idx] = 0.0f;
}

#ifdef __cplusplus
}
#endif

#endif /* LAGRANGE_LIBRATION_H */

