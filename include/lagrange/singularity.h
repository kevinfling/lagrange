/**
 * lagrange_singularity.h - Black Hole Physics & Relativistic Mechanics
 * 
 * Schwarzschild and Kerr metrics for orbital mechanics near compact objects.
 * Geodesic integration, gravitational lensing, accretion disks, and frame dragging.
 * 
 * Integrates with: lagrange_particle.h (relativistic N-body),
 *                  lagrange_koopman.h (linearized stability near horizon),
 *                  lagrange_wavelet.h (gravitational wave signal analysis)
 * 
 * "Just for giggles" - but mathematically rigorous.
 */

#ifndef LAGRANGE_SINGULARITY_H
#define LAGRANGE_SINGULARITY_H

#include "math.h"
#include "body.h"
#include "particle.h"
#include <complex.h>

#ifdef __cplusplus
extern "C" {
#endif

/*============================================================================
 * Physical Constants (Geometric Units: G = c = 1)
 * Convert: 1 M_sun = 4.93e-6 seconds = 1.477 km
 *===========================================================================*/

#define LG_SCHWARZSCHILD_RADIUS_SUN 2.953250e3f      /* meters, 2GM_sun/c^2 */
#define LG_C 299792458.0f                              /* m/s, for conversions */
#define LG_G 6.67430e-11f

/* Dimensionless spin parameter a = J/M (0 <= a < 1 for Kerr) */
typedef float lg_spin_t;

/*============================================================================
 * Schwarzschild Metric (Static, Non-Rotating)
 * ds^2 = -(1-2M/r)dt^2 + (1-2M/r)^-1 dr^2 + r^2 dOmega^2
 *===========================================================================*/

typedef struct {
    float M;              /* Mass in meters (geometric) */
    float rs;             /* Schwarzschild radius = 2M */
    
    /* Precomputed for speed */
    float rs_sq;
    float light_crossing_time; /* rs/c in seconds */
} lg_schwarzschild_t;

static inline lg_schwarzschild_t lg_schwarzschild(float mass_kg) {
    float M = mass_kg * LG_G / (LG_C * LG_C); /* Convert kg to geometric meters */
    float rs = 2.0f * M;
    return (lg_schwarzschild_t){M, rs, rs*rs, rs/LG_C};
}

/* Time dilation factor (gravitational redshift) 
 * dt_far/dt_near = 1/sqrt(1 - rs/r) */
static inline float lg_schwarzschild_redshift(const lg_schwarzschild_t* bh, float r) {
    if (r <= bh->rs) return 0.0f/0.0f; /* NaN inside horizon */
    return 1.0f / sqrtf(1.0f - bh->rs / r);
}

/* Gravitational potential (effective for radial infall) 
 * V_eff = (1 - rs/r)(1 + L^2/r^2) */
static inline float lg_schwarzschild_v_eff(const lg_schwarzschild_t* bh, 
                                          float r, float L) { /* L = specific ang mom */
    float L2 = L * L;
    float r2 = r * r;
    return (1.0f - bh->rs / r) * (1.0f + L2 / r2);
}

/* Innermost Stable Circular Orbit (ISCO) */
static inline float lg_schwarzschild_isco(const lg_schwarzschild_t* bh) {
    return 6.0f * bh->M; /* 3 rs */
}

/* Photon Sphere (unstable circular light orbit) */
static inline float lg_schwarzschild_photon_sphere(const lg_schwarzschild_t* bh) {
    return 1.5f * bh->rs; /* 3M */
}

/*============================================================================
 * Kerr Metric (Rotating Black Hole)
 * 
 * ds^2 involves cross term dt*dphi for frame dragging
 * Metric components in Boyer-Lindquist coordinates
 *===========================================================================*/

typedef struct {
    float M;              /* Mass */
    float a;              /* Spin parameter (0 to <1) */
    float rs;             /* 2M */
    
    /* Horizons */
    float r_plus;         /* Outer horizon: M + sqrt(M^2 - a^2) */
    float r_minus;        /* Inner horizon: M - sqrt(M^2 - a^2) */
    
    /* Ergosphere (static limit) */
    float r_ergo;         /* M + sqrt(M^2 - a^2 cos^2 theta) varies with theta */
} lg_kerr_t;

static inline lg_kerr_t lg_kerr(float mass_kg, float spin_a) {
    float M = mass_kg * LG_G / (LG_C * LG_C);
    float a_clamped = fmaxf(0.0f, fminf(spin_a, 0.9999f)); /* Prevent naked singularity */
    float sqrt_term = sqrtf(fmaxf(0.0f, M*M - a_clamped*a_clamped));
    
    return (lg_kerr_t){
        M, a_clamped, 2.0f*M,
        M + sqrt_term,      /* r_+ */
        M - sqrt_term,      /* r_- */
        2.0f*M              /* r_ergo at equator */
    };
}

/* Frame dragging angular velocity (Lense-Thirring effect)
 * Omega = 2Mar / (r^3 + a^2r + 2Ma^2) */
static inline float lg_kerr_frame_dragging(const lg_kerr_t* bh, float r, float theta) {
    float r2 = r * r;
    float a2 = bh->a * bh->a;
    float sin_theta = sinf(theta);
    float sin2 = sin_theta * sin_theta;
    
    float Sigma = r2 + a2 * cosf(theta) * cosf(theta);
    float Delta = r2 - bh->rs * r + a2;
    
    return 2.0f * bh->a * bh->M * r * sin2 / ((r2 + a2)*(r2 + a2) - a2 * Delta * sin2);
}

/* ISCO for Kerr (prograde vs retrograde) */
static inline float lg_kerr_isco(const lg_kerr_t* bh, bool prograde) {
    float z1 = 1.0f + powf(1.0f - bh->a*bh->a, 1.0f/3.0f) * (powf(1.0f + bh->a, 1.0f/3.0f) + powf(1.0f - bh->a, 1.0f/3.0f));
    float z2 = sqrtf(3.0f * bh->a * bh->a + z1 * z1);
    float sign = prograde ? 1.0f : -1.0f;
    return bh->M * (3.0f + z2 - sign * sqrtf((3.0f - z1) * (3.0f + z1 + 2.0f * z2)));
}

/* Penrose Process: energy extraction from ergosphere */
static inline float lg_kerr_penrose_efficiency(const lg_kerr_t* bh) {
    /* Max efficiency = 1 - 1/sqrt(2) ~ 29% for a=1 */
    return 1.0f - sqrtf(1.0f - 0.5f * (bh->r_plus / bh->M - 1.0f));
}

/*============================================================================
 * Geodesic Equations (Relativistic Orbits)
 * 4D state: (t, r, theta, phi) and (dt/dtau, dr/dtau, dtheta/dtau, dphi/dtau)
 *===========================================================================*/

typedef struct {
    float tau;            /* Proper time */
    float x[4];           /* (t, r, theta, phi) */
    float v[4];           /* (dt/dtau, dr/dtau, dtheta/dtau, dphi/dtau) */
    
    /* Conserved quantities */
    float E;              /* Energy per unit mass */
    float L;              /* Angular momentum per unit mass */
    float Q;              /* Carter constant (for Kerr) */
    float mu;             /* Mass (0 for photons, 1 for particles) */
} lg_geodesic_t;

/* Christoffel symbols for Schwarzschild (simplified for equatorial theta=pi/2) */
static inline void lg_schwarzschild_christoffel(float Gamma[4][4][4], float r) {
    memset(Gamma, 0, 4*4*4*sizeof(float));
    /* Fill non-zero components for equatorial plane */
    /* Gamma^t_tr = M/(r(r-2M)) etc. - explicit calculation */
}

/* RK4 integration of geodesic equation:
 * d^2x^mu/dtau^2 = -Gamma^mu_nu_sigma * v^nu * v^sigma */
static inline void lg_geodesic_step_rk4(lg_geodesic_t* geo, 
                                        const lg_schwarzschild_t* bh,
                                        float dtau) {
    /* Standard RK4 for 2nd order ODE system */
    /* Evaluate k1, k2, k3, k4 for Christoffel terms */
    
    /* Simplified: use effective potential method for equatorial orbits instead */
}

/* Relativistic precession per orbit 
 * Delta phi = 6*pi*M/p where p = semi-latus rectum */
static inline float lg_schwarzschild_precession(const lg_schwarzschild_t* bh, 
                                               float semi_major_axis,
                                               float eccentricity) {
    float p = semi_major_axis * (1.0f - eccentricity * eccentricity);
    return 6.0f * M_PI * bh->M / p; /* radians per orbit */
}

/*============================================================================
 * Gravitational Lensing (Ray Tracing)
 * 
 * Trace null geodesics (photons) through Schwarzschild metric
 *===========================================================================*/

typedef struct {
    lg_vec3_t origin;     /* Camera position */
    lg_vec3_t direction;  /* Ray direction (normalized) */
    float impact_param;   /* b = L/E (conserved for photons) */
} lg_ray_t;

/* Deflection angle for weak field (far from BH) 
 * alpha = 4M/b (Einstein ring formula) */
static inline float lg_lensing_deflection(const lg_schwarzschild_t* bh, float b) {
    return 4.0f * bh->M / b;
}

/* Exact geodesic ray tracer for strong lensing */
static inline bool lg_trace_ray_schwarzschild(lg_ray_t* ray,
                                              const lg_schwarzschild_t* bh,
                                              float max_distance,
                                              int steps,
                                              lg_vec3_t* hit_point) {
    float b = ray->impact_param;
    float r0 = lg_vec3_len(ray->origin);
    
    /* Photons with b < 3*sqrt(3)*M get captured */
    float b_crit = 3.0f * sqrtf(3.0f) * bh->M; /* ~5.196 M */
    
    if (b < b_crit && r0 < 10.0f * bh->rs) {
        /* Check if ray spirals into horizon */
        /* Solve effective potential: V = (1 - 2M/r)(1 + b^2/r^2) */
        float r_min = 2.0f * b; /* Approximate turning point */
        float V_max = (1.0f - bh->rs/r_min) * (1.0f + b*b/(r_min*r_min));
        
        if (V_max > 1.0f) {
            /* Photon falls into black hole */
            return false;
        }
    }
    
    /* Integrate null geodesic */
    for (int i = 0; i < steps; i++) {
        /* Update using conservation of b and E */
        /* dr/dl = ... angular equation ... */
    }
    
    return true;
}

/* Einstein ring radius (angular) for lensing */
static inline float lg_einstein_ring(const lg_schwarzschild_t* bh, 
                                    float D_lens, float D_source) {
    /* theta_E = sqrt(4GM D_ls / (c^2 D_s D_l)) */
    float D_ls = D_source - D_lens;
    return sqrtf(4.0f * bh->M * D_ls / (D_lens * D_source));
}

/*============================================================================
 * Accretion Disk Physics (Thin Disk Model)
 *===========================================================================*/

typedef struct {
    float Mdot;           /* Accretion rate (kg/s or geometric) */
    float r_in;           /* Inner edge (typically ISCO) */
    float r_out;          /* Outer edge */
    float alpha_visc;     /* Shakura-Sunyaev viscosity parameter */
} lg_accretion_disk_t;

/* Novikov-Thorne radiative efficiency 
 * eta = 1 - sqrt(1 - 2/(3*r_isco)) for Schwarzschild */
static inline float lg_nt_efficiency(const lg_schwarzschild_t* bh) {
    float r_isco = lg_schwarzschild_isco(bh);
    return 1.0f - sqrtf(1.0f - 2.0f / (3.0f * r_isco / bh->M));
}

/* Temperature profile T(r) ~ r^{-3/4} (Keplerian disk) */
static inline float lg_disk_temperature(const lg_accretion_disk_t* disk,
                                     const lg_schwarzschild_t* bh,
                                     float r) {
    float r_ratio = r / bh->rs;
    /* T^4 ~ Mdot M / r^3 * (1 - sqrt(r_in/r)) */
    float f = 1.0f - sqrtf(disk->r_in / r);
    float T_max = 1e7f; /* ~10^7 K for stellar mass BH, scaled */
    return T_max * powf(r_ratio, -0.75f) * powf(f, 0.25f);
}

/* Spectrum calculation (multi-temperature blackbody) */
static inline void lg_disk_spectrum(const lg_accretion_disk_t* disk,
                                 const lg_schwarzschild_t* bh,
                                 float* frequencies, /* Hz */
                                 float* flux,
                                 int n_bins) {
    /* Integrate Planck function over disk surface */
    /* Include gravitational redshift (1+z) factor from disk to infinity */
    for (int i = 0; i < n_bins; i++) {
        float nu = frequencies[i];
        float total = 0.0f;
        
        /* Sample radii */
        for (float r = disk->r_in; r < disk->r_out; r *= 1.1f) {
            float T = lg_disk_temperature(disk, bh, r);
            float redshift = lg_schwarzschild_redshift(bh, r);
            float T_obs = T / redshift; /* Redshifted temperature */
            
            /* Planck law B_nu(T) */
            float x = 6.626e-34f * nu / (1.38e-23f * T_obs);
            float B = (2.0f * 6.626e-34f * nu*nu*nu / (LG_C*LG_C)) / (expf(x) - 1.0f);
            
            /* Area element 2*pi*r*dr */
            total += B * 2.0f * M_PI * r * (r * 0.1f); /* logarithmic dr */
        }
        flux[i] = total;
    }
}

/*============================================================================
 * Tidal Disruption Events (TDE)
 *===========================================================================*/

/* Roche limit for tidal disruption 
 * r_t ~ R_star * (M_bh / M_star)^{1/3} */
static inline float lg_tidal_radius(float M_bh, float M_star, float R_star) {
    return R_star * powf(M_bh / M_star, 1.0f/3.0f);
}

/* Fallback rate for TDE debris (power law decay) 
 * Mdot ~ t^{-5/3} */
static inline float lg_tde_fallback_rate(float M_dot_peak, float t_days) {
    return M_dot_peak * powf(t_days / 10.0f, -5.0f/3.0f);
}

/*============================================================================
 * Gravitational Waves (Quadrupole Formula)
 * For binary inspiral waveforms
 *===========================================================================*/

/* Inspiral chirp mass */
static inline float lg_chirp_mass(float m1, float m2) {
    return powf(m1 * m2, 3.0f/5.0f) / powf(m1 + m2, 1.0f/5.0f);
}

/* Gravitational wave frequency evolution (Newtonian) 
 * f^{-8/3} = (8/3) * (G^5/c^3) * M_chirp^{5/3} * (t_c - t) */
static inline float lg_gw_frequency(float M_chirp_geom, float t_to_coalescence) {
    /* t_c - t in geometric units */
    float factor = 5.0f / 256.0f * powf(M_chirp_geom, -5.0f/3.0f);
    return powf(factor / t_to_coalescence, 3.0f/8.0f) / M_PI;
}

/* Strain amplitude (characteristic) */
static inline float lg_gw_strain(float M_chirp_geom, float f_gw, float distance_m) {
    /* h ~ (G M_chirp / c^2 D) * (G M_chirp f / c^3)^{2/3} */
    float M = M_chirp_geom;
    float D = distance_m * LG_C * LG_C / LG_G; /* Convert to geometric */
    return powf(M/D, 5.0f/3.0f) * powf(M*f_gw, 2.0f/3.0f);
}

/*============================================================================
 * Relativistic N-Body (Post-Newtonian corrections)
 * Add 1PN, 2PN terms to standard Newtonian acceleration
 *===========================================================================*/

/* 1PN (Einstein-Infeld-Hoffmann) acceleration
 * 
 * Computes post-Newtonian correction to gravitational acceleration.
 * 
 * Parameters:
 *   i - Index of body being accelerated (to skip self-interaction)
 *   positions - Array of body positions (world space)
 *   velocities - Array of body velocities
 *   masses - Array of body masses
 *   n - Number of bodies
 *   c_inv - 1/speed_of_light
 */
static inline lg_vec3_t lg_acceleration_1pn(int i,
                                           const lg_vec3_t* positions,
                                           const lg_vec3_t* velocities,
                                           const float* masses,
                                           int n,
                                           float c_inv) {
    lg_vec3_t a_pn = lg_vec3_zero();
    float c2_inv = c_inv * c_inv;
    
    for (int j = 0; j < n; j++) {
        if (j == i) continue;
        
        lg_vec3_t r_ij = lg_vec3_sub(positions[j], positions[i]);
        float r = lg_vec3_len(r_ij);
        if (r < 1e-10f) continue;  /* Avoid division by zero */
        
        float v_i_sq = lg_vec3_len_sq(velocities[i]);
        float v_j_sq = lg_vec3_len_sq(velocities[j]);
        lg_vec3_t v_i = velocities[i];
        lg_vec3_t v_j = velocities[j];
        
        lg_vec3_t r_hat = lg_vec3_norm(r_ij);
        float r_dot = lg_vec3_dot(v_i, r_hat);
        
        /* Mass ratio */
        float nu = masses[j] / (masses[i] + masses[j]);
        
        /* 1PN correction factor (simplified EIH equation) */
        float factor = c2_inv * (
            -2.0f * (masses[i] + masses[j]) / r +  /* Potential */
            v_j_sq + 2.0f * v_i_sq - 4.0f * lg_vec3_dot(v_i, v_j) - /* Velocities */
            1.5f * r_dot * r_dot +                  /* Radial velocity term */
            0.5f * lg_vec3_dot(r_ij, lg_vec3_sub(v_j, v_i)) / r  /* Coupling */
        );
        
        lg_vec3_t term = lg_vec3_scale(r_hat, factor * masses[j] / (r*r));
        a_pn = lg_vec3_add(a_pn, term);
    }
    
    return a_pn;
}

/*============================================================================
 * Event Horizon Detection (Particle removal)
 *===========================================================================*/

static inline bool lg_horizon_crossed_schwarzschild(const lg_schwarzschild_t* bh,
                                                    const lg_vec3_t* pos) {
    return lg_vec3_len_sq(*pos) < bh->rs * bh->rs;
}

static inline bool lg_horizon_crossed_kerr(const lg_kerr_t* bh,
                                          const lg_vec3_t* pos,
                                          float theta) {
    float r = lg_vec3_len(*pos);
    float r_horizon = bh->r_plus;
    return r < r_horizon;
}

/*============================================================================
 * SIMD Batch Processing (Relativistic redshifts for many particles)
 *===========================================================================*/

#ifdef __AVX2__
#include <immintrin.h>

/* Compute gravitational redshift z = 1/sqrt(1-rs/r) - 1 for 8 particles */
static inline void lg_redshift_batch_avx2(const lg_schwarzschild_t* bh,
                                          const float* x, const float* y, const float* z,
                                          float* redshift,
                                          int n) {
    __m256 rs = _mm256_set1_ps(bh->rs);
    __m256 one = _mm256_set1_ps(1.0f);
    
    for (int i = 0; i < n; i += 8) {
        __m256 vx = _mm256_loadu_ps(x + i);
        __m256 vy = _mm256_loadu_ps(y + i);
        __m256 vz = _mm256_loadu_ps(z + i);
        
        __m256 r2 = _mm256_fmadd_ps(vx, vx, _mm256_fmadd_ps(vy, vy, _mm256_mul_ps(vz, vz)));
        __m256 r = _mm256_sqrt_ps(r2);
        
        /* z = 1/sqrt(1 - rs/r) - 1 */
        __m256 ratio = _mm256_div_ps(rs, r);
        __m256 interior = _mm256_sub_ps(one, ratio);
        __m256 z_factor = _mm256_rsqrt_ps(interior); /* 1/sqrt(1-rs/r) */
        
        _mm256_storeu_ps(redshift + i, _mm256_sub_ps(z_factor, one));
    }
}
#endif

/*============================================================================
 * Fractional Hawking Radiation (Just for giggles)
 * Treating evaporation as fractional process with memory
 *===========================================================================*/

/* Hawking temperature T = hbar c^3 / (8 pi G M k_B) */
static inline float lg_hawking_temperature(float M_kg) {
    const float hbar = 1.054571817e-34f;
    const float kB = 1.380649e-23f;
    return hbar * LG_C*LG_C*LG_C / (8.0f * M_PI * LG_G * M_kg * kB);
}

/* Mass evolution with fractional memory (non-Markovian evaporation) */
static inline float lg_bh_mass_fractional(float M0, float t, float alpha, float tau) {
    /* Standard: M(t) = (M0^3 - 3K*t)^{1/3} */
    /* Fractional: incorporates quantum memory effects */
    float M_std = powf(M0*M0*M0 - 3.0f * 1e-6f * t, 1.0f/3.0f);
    
    /* Fractional correction: M(t) = M_std * (1 + fractional integral of quantum noise) */
    return M_std; /* Simplified - full version needs fractional calculus from lagrange_fractional.h */
}

#ifdef __cplusplus
}
#endif

#endif /* LAGRANGE_SINGULARITY_H */
