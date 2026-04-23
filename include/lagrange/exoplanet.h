/**
 * lagrange_exoplanet.h - Procedural Generation from NASA Exoplanet Archive
 * 
 * Statistical synthesis of planetary systems using Gaussian Copula sampling
 * to preserve correlations between mass, period, eccentricity, and stellar host.
 * 
 * Features:
 *   - Gaussian Copula with Cholesky decomposition for rank correlations
 *   - Physical constraint enforcement (Hill stability, resonant chains)
 *   - Marginal distributions: Log-normal mass, Beta eccentricity, Rayleigh inclination
 *   - Resonant chain injection (period ratios 3:2, 2:1 from Kepler data)
 *   - Export to lg_particle_system_t for direct simulation
 * 
 * Usage: Randomized Solar System button → lg_exo_generate_system() → N-body sim
 */

#ifndef LAGRANGE_EXOPLANET_H
#define LAGRANGE_EXOPLANET_H

#include "math.h"
#include "particle.h"
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/*============================================================================
 * Data Structures (mirrors IPAC PS columns)
 *===========================================================================*/

typedef struct {
    char pl_name[32];         /* Planet name (for tracking) */
    float pl_mass;            /* Jupiter masses */
    float pl_radius;          /* Jupiter radii */
    float pl_orbper;          /* Orbital period (days) */
    float pl_orbsmax;         /* Semi-major axis (AU) */
    float pl_orbeccen;        /* Eccentricity [0,1] */
    float pl_orbincl;         /* Inclination (deg) */
    float pl_orblper;        /* Longitude of periastron */
    float pl_eqt;             /* Equilibrium temperature (K) */
    
    /* Host star */
    float st_mass;            /* Solar masses */
    float st_rad;             /* Solar radii */
    float st_teff;            /* Effective temp (K) */
    float st_met;             /* Metallicity [Fe/H] */
    float st_age;             /* Age (Gyr) */
} lg_exoplanet_record_t;

/* Statistical model parameters derived from archive */
typedef struct {
    /* Marginal distribution parameters */
    struct { float mu, sigma; } log_mass;       /* Log-normal for masses */
    struct { float alpha, beta; } eccen;        /* Beta distribution for e */
    struct { float sigma; } incl;               /* Rayleigh for inclinations */
    struct { float min, max; } log_period;      /* Uniform-ish in log for cold planets */
    
    /* Copula correlation matrix (Cholesky factor L where Sigma = L*L^T) */
    /* Variables: [log_star_mass, log_pl_mass, log_period, eccen, metallicity] */
    float L[5][5];            /* Lower triangular Cholesky factor */
    
    /* Multi-planet system statistics */
    float mean_planets_per_star;
    float resonant_fraction;  /* Fraction in mean-motion resonances */
    float hot_jupiter_rate;   /* Occurrence rate of P < 10 day giants */
} lg_exo_model_t;

/*============================================================================
 * Random Number Generation (pcg32 for portability)
 *===========================================================================*/

typedef struct { uint64_t state; uint64_t inc; } lg_pcg32_t;

static inline uint32_t lg_pcg32_random(lg_pcg32_t* rng) {
    uint64_t oldstate = rng->state;
    rng->state = oldstate * 6364136223846793005ULL + (rng->inc|1);
    uint32_t xorshifted = ((oldstate >> 18u) ^ oldstate) >> 27u;
    uint32_t rot = oldstate >> 59u;
    return (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
}

static inline float lg_pcg32_uniform(lg_pcg32_t* rng) {
    return (lg_pcg32_random(rng) >> 8) / 16777216.0f; /* [0,1) with 24-bit prec */
}

/* Box-Muller for normal */
static inline float lg_rand_normal(lg_pcg32_t* rng) {
    float u1 = lg_pcg32_uniform(rng);
    float u2 = lg_pcg32_uniform(rng);
    if (u1 < 1e-7f) u1 = 1e-7f;
    return sqrtf(-2.0f * logf(u1)) * cosf(2.0f * M_PI * u2);
}

/*============================================================================
 * Distributions (CDFs and inverse CDFs for copula transforms)
 *===========================================================================*/

/* Log-normal */
static inline float lg_cdf_lognormal(float x, float mu, float sigma) {
    if (x <= 0) return 0.0f;
    return 0.5f + 0.5f * erff((logf(x) - mu) / (sigma * sqrtf(2.0f)));
}

/* Inverse standard normal CDF (Acklam rational approximation) */
static inline float lg_icdf_normal(float p) {
    if (p <= 0.0f) return -INFINITY;
    if (p >= 1.0f) return INFINITY;

    float a1 = -39.6968302866538f;
    float a2 = 220.946098424521f;
    float a3 = -275.928510446969f;
    float a4 = 138.357751867269f;
    float a5 = -30.6647980661472f;
    float a6 = 2.50662827745924f;

    float b1 = -54.4760987982241f;
    float b2 = 161.585836858041f;
    float b3 = -155.698979859887f;
    float b4 = 66.8013118877197f;
    float b5 = -13.2806811558857f;

    float c1 = -0.00778489400243029f;
    float c2 = -0.322396458441136f;
    float c3 = -2.40075827716184f;
    float c4 = -2.54973253934373f;
    float c5 = 4.37466414146497f;
    float c6 = 2.93816398269878f;

    float d1 = 0.00778469570904146f;
    float d2 = 0.32246712907004f;
    float d3 = 2.445134137143f;
    float d4 = 3.75440866190742f;

    float p_low  = 0.02425f;
    float p_high = 1.0f - p_low;

    float q, r;
    if (p < p_low) {
        q = sqrtf(-2.0f * logf(p));
        return (((((c1*q+c2)*q+c3)*q+c4)*q+c5)*q+c6) /
               ((((d1*q+d2)*q+d3)*q+d4)*q+1.0f);
    } else if (p <= p_high) {
        q = p - 0.5f;
        r = q*q;
        return (((((a1*r+a2)*r+a3)*r+a4)*r+a5)*r+a6)*q /
               (((((b1*r+b2)*r+b3)*r+b4)*r+b5)*r+1.0f);
    } else {
        q = sqrtf(-2.0f * logf(1.0f - p));
        return -(((((c1*q+c2)*q+c3)*q+c4)*q+c5)*q+c6) /
                ((((d1*q+d2)*q+d3)*q+d4)*q+1.0f);
    }
}

static inline float lg_icdf_lognormal(float p, float mu, float sigma) {
    if (p <= 0.0f) return 0.0f;
    if (p >= 1.0f) return INFINITY;
    float z = lg_icdf_normal(p);
    return expf(mu + sigma * z);
}

/* Incomplete beta function I_x(a,b) using Lentz's continued fraction */
static inline double lg_incbeta(double a, double b, double x)
{
    #define LG_INCBETA_MAXITER 300
    #define LG_INCBETA_STOP    1.0e-14
    #define LG_INCBETA_TINY    1.0e-50

    if (x <= 0.0) return 0.0;
    if (x >= 1.0) return 1.0;
    if (a <= 0.0 || b <= 0.0) return 0.0 / 0.0;

    int flipped = 0;
    if (x > (a + 1.0) / (a + b + 2.0)) {
        flipped = 1;
        double t = a; a = b; b = t;
        x = 1.0 - x;
    }

    const double log_beta_ab = lgamma(a) + lgamma(b) - lgamma(a + b);
    const double log_front = a * log(x) + b * log1p(-x) - log_beta_ab;
    const double front = exp(log_front) / a;

    double f = 1.0;
    double c = 1.0;
    double d = 0.0;
    int i;
    for (i = 0; i < LG_INCBETA_MAXITER; ++i) {
        const int m = i / 2;
        double num;
        if (i == 0) {
            num = 1.0;
        } else if (i % 2 == 0) {
            num = (m * (b - m) * x) / ((a + 2.0 * m - 1.0) * (a + 2.0 * m));
        } else {
            num = -((a + m) * (a + b + m) * x) / ((a + 2.0 * m) * (a + 2.0 * m + 1.0));
        }

        d = 1.0 + num * d;
        if (fabs(d) < LG_INCBETA_TINY) d = LG_INCBETA_TINY;
        d = 1.0 / d;

        c = 1.0 + num / c;
        if (fabs(c) < LG_INCBETA_TINY) c = LG_INCBETA_TINY;

        const double cd = c * d;
        f *= cd;

        if (fabs(1.0 - cd) < LG_INCBETA_STOP) {
            break;
        }
    }

    double val = front * (f - 1.0);
    return flipped ? (1.0 - val) : val;

    #undef LG_INCBETA_MAXITER
    #undef LG_INCBETA_STOP
    #undef LG_INCBETA_TINY
}

/* Beta distribution (for eccentricity) */
static inline float lg_cdf_beta(float x, float alpha, float beta) {
    return (float)lg_incbeta((double)alpha, (double)beta, (double)x);
}

static inline float lg_icdf_beta(float p, float alpha, float beta) {
    if (p <= 0.0f) return 0.0f;
    if (p >= 1.0f) return 1.0f;
    if (alpha <= 0.0f || beta <= 0.0f) return 0.0f / 0.0f;

    /* Binary search for x such that lg_cdf_beta(x, alpha, beta) == p */
    double lo = 0.0, hi = 1.0;
    double target = (double)p;
    for (int iter = 0; iter < 60; ++iter) {
        double mid = (lo + hi) * 0.5;
        double cdf = lg_incbeta((double)alpha, (double)beta, mid);
        if (cdf < target) {
            lo = mid;
        } else {
            hi = mid;
        }
        if (hi - lo < 1.0e-12) break;
    }
    return (float)((lo + hi) * 0.5);
}

/* Rayleigh for inclination */
static inline float lg_icdf_rayleigh(float p, float sigma) {
    return sigma * sqrtf(-2.0f * logf(1.0f - p));
}

/*============================================================================
 * Gaussian Copula Sampling
 * 
 * 1. Generate Z ~ N(0,I) independent standard normals
 * 2. Transform to correlated: X = L * Z (L is Cholesky of target correlation)
 * 3. Transform to uniform: U = Phi(X) (standard normal CDF)
 * 4. Transform to target marginals: Y = F^{-1}(U)
 *===========================================================================*/

static inline void lg_copula_sample(lg_pcg32_t* rng, 
                                    const lg_exo_model_t* model,
                                    float* out_vars /* [log_st_mass, log_pl_mass, log_period, ecc, met] */) {
    /* Step 1: Generate independent N(0,1) */
    float z[5];
    for (int i = 0; i < 5; i++) z[i] = lg_rand_normal(rng);
    
    /* Step 2: Apply Cholesky L (lower triangular) to get correlated normals */
    float x[5] = {0};
    for (int i = 0; i < 5; i++) {
        for (int j = 0; j <= i; j++) {
            x[i] += model->L[i][j] * z[j];
        }
    }
    
    /* Step 3: Standard normal CDF approximation (error function) */
    float u[5];
    for (int i = 0; i < 5; i++) {
        u[i] = 0.5f + 0.5f * erff(x[i] / sqrtf(2.0f));
    }
    
    /* Step 4: Inverse transform to marginals */
    out_vars[0] = lg_icdf_lognormal(u[0], 0.0f, 0.2f);      /* Stellar mass ~ 1M_sun */
    out_vars[1] = lg_icdf_lognormal(u[1], -2.5f, 1.0f);     /* Planet mass (mix of super-Earths and giants) */
    out_vars[2] = model->log_period.min + u[2] * (model->log_period.max - model->log_period.min);
    out_vars[3] = lg_icdf_beta(u[3], 1.0f, 3.0f);           /* Eccentricity, typical e~0.3 */
    out_vars[4] = u[4] * 0.5f - 0.25f;                      /* Metallicity [-0.25, 0.25] */
}

/*============================================================================
 * Cholesky Decomposition (for correlation matrix)
 * Compute L such that L*L^T = Sigma (correlation matrix)
 *===========================================================================*/

static inline void lg_cholesky_5x5(float L[5][5], const float Sigma[5][5]) {
    memset(L, 0, 25*sizeof(float));
    for (int i = 0; i < 5; i++) {
        for (int j = 0; j <= i; j++) {
            float sum = 0.0f;
            for (int k = 0; k < j; k++) sum += L[i][k] * L[j][k];
            
            if (i == j) {
                L[i][i] = sqrtf(fmaxf(0.001f, Sigma[i][i] - sum));
            } else {
                L[i][j] = (Sigma[i][j] - sum) / L[j][j];
            }
        }
    }
}

/*============================================================================
 * Physical Constraints & Stability
 *===========================================================================*/

/* Hill stability criterion: minimum separation for circular orbits */
static inline float lg_hill_radius(float a, float e, float m_planet, float m_star) {
    return a * (1.0f - e) * powf(m_planet / (3.0f * m_star), 1.0f/3.0f);
}

/* Check if proposed system is Hill stable */
static inline bool lg_system_hill_stable(const lg_exoplanet_record_t* planets,
                                       int n, float safety_factor) {
    for (int i = 0; i < n-1; i++) {
        float a_in = planets[i].pl_orbsmax;
        float a_out = planets[i+1].pl_orbsmax;
        float rh_in = lg_hill_radius(a_in, planets[i].pl_orbeccen, 
                                    planets[i].pl_mass, planets[i].st_mass);
        float rh_out = lg_hill_radius(a_out, planets[i+1].pl_orbeccen,
                                     planets[i+1].pl_mass, planets[i+1].st_mass);
        
        /* Separation must be > Delta * R_Hill */
        if ((a_out - a_in) < safety_factor * (rh_in + rh_out)) {
            return false;
        }
    }
    return true;
}

/* Inject mean-motion resonances (period ratios) */
static inline void lg_resonant_chain(lg_exoplanet_record_t* planets, int n) {
    /* Common resonances: 3:2, 2:1, 4:3 */
    int resonances[] = {3, 2, 2, 1, 4, 3};
    int n_res = 3;
    
    float period_0 = planets[0].pl_orbper;
    for (int i = 1; i < n && i < n_res*2; i += 2) {
        int idx = (lg_pcg32_random(NULL) >> 30) % n_res; /* Pick random resonance */
        int p = resonances[2*idx];
        int q = resonances[2*idx+1];
        
        /* P_{i+1} / P_i = p/q with small scatter */
        float jitter = 1.0f + 0.01f * (lg_pcg32_uniform(NULL) - 0.5f);
        planets[i].pl_orbper = period_0 * powf((float)p/q, i) * jitter;
        planets[i].pl_orbsmax = powf(planets[i].st_mass * 
                                    (planets[i].pl_orbper/365.25f) * 
                                    (planets[i].pl_orbper/365.25f), 1.0f/3.0f); /* Kepler's 3rd */
    }
}

/*============================================================================
 * System Generation API ("Randomized Solar System" Button)
 *===========================================================================*/

typedef struct {
    lg_exoplanet_record_t star;
    lg_exoplanet_record_t* planets;
    int n_planets;
    bool is_resonant_chain;
} lg_exosystem_t;

/* Initialize model with fitted parameters from IPAC archive */
static inline lg_exo_model_t lg_exo_model_default(void) {
    lg_exo_model_t m = {
        .log_mass = {-2.3f, 1.2f},      /* Super-Earth peak + Jupiter tail */
        .eccen = {1.0f, 3.0f},          /* Beta(1,3) peaks at low e */
        .incl = {2.0f},                 /* Rayleigh sigma ~2 deg for aligned */
        .log_period = {0.3f, 4.0f},     /* 2 days to 54 years */
        .mean_planets_per_star = 2.5f,
        .resonant_fraction = 0.4f,
        .hot_jupiter_rate = 0.01f
    };
    
    /* Correlation matrix (simplified from archive data) 
     * [st_mass, pl_mass, period, ecc, met]
     * Positive: st_mass correlated with pl_mass (giant planets around big stars)
     * Negative: period vs ecc (circularization at short periods) */
    float Sigma[5][5] = {
        {1.0f,  0.3f,  0.1f, -0.1f,  0.2f},
        {0.3f,  1.0f,  0.4f,  0.0f,  0.3f},
        {0.1f,  0.4f,  1.0f, -0.5f,  0.0f},  /* -0.5: short period = low e */
        {-0.1f, 0.0f, -0.5f,  1.0f, -0.1f},
        {0.2f,  0.3f,  0.0f, -0.1f,  1.0f}
    };
    
    lg_cholesky_5x5(m.L, Sigma);
    return m;
}

/* Generate one complete system */
static inline lg_exosystem_t lg_exo_generate_system(lg_pcg32_t* rng, 
                                                  const lg_exo_model_t* model) {
    lg_exosystem_t sys = {0};
    
    /* Sample star first */
    float vars[5];
    lg_copula_sample(rng, model, vars);
    
    sys.star.st_mass = vars[0];
    sys.star.st_rad = powf(vars[0], 0.8f); /* R ~ M^0.8 for main sequence */
    sys.star.st_teff = 5778.0f * powf(vars[0], 0.5f); /* Rough scaling */
    sys.star.st_met = vars[4];
    
    /* Poisson sample number of planets */
    float lam = model->mean_planets_per_star;
    int n_pl = 0;
    float L = expf(-lam);
    float p = 1.0f;
    do {
        n_pl++;
        p *= lg_pcg32_uniform(rng);
    } while (p > L && n_pl < 10);
    sys.n_planets = n_pl - 1;
    
    if (sys.n_planets == 0) return sys;
    
    sys.planets = (lg_exoplanet_record_t*)calloc(sys.n_planets, sizeof(lg_exoplanet_record_t));
    
    /* Decide if resonant chain (compact systems like TRAPPIST-1) */
    sys.is_resonant_chain = (lg_pcg32_uniform(rng) < model->resonant_fraction) && (sys.n_planets >= 3);
    
    /* Generate planets inner to outer */
    for (int i = 0; i < sys.n_planets; i++) {
        lg_copula_sample(rng, model, vars);
        
        sys.planets[i].st_mass = sys.star.st_mass;
        sys.planets[i].pl_mass = vars[1];  /* Jupiter masses */
        sys.planets[i].pl_orbper = expf(vars[2]); /* days */
        sys.planets[i].pl_orbeccen = vars[3];
        sys.planets[i].pl_orbincl = lg_icdf_rayleigh(lg_pcg32_uniform(rng), model->incl.sigma);
        sys.planets[i].pl_orblper = 360.0f * lg_pcg32_uniform(rng);
        
        /* Kepler: a^3/P^2 = M_star (AU^3/year^2 = M_sun) */
        float P_yr = sys.planets[i].pl_orbper / 365.25f;
        sys.planets[i].pl_orbsmax = powf(P_yr * P_yr * sys.star.st_mass, 1.0f/3.0f);
        
        /* Hot Jupiter check: if P < 10 days and M > 0.3 Mjup */
        if (sys.planets[i].pl_orbper < 10.0f && sys.planets[i].pl_mass > 0.3f) {
            /* Force circular orbit (tidal circularization) */
            sys.planets[i].pl_orbeccen *= 0.1f;
        }
        
        /* Radius from mass (crude): R = R_jup * (M/M_jup)^0.5 for giants, 
         * or constant for rocky */
        if (sys.planets[i].pl_mass < 0.05f) {
            sys.planets[i].pl_radius = 0.1f; /* ~1 Earth radii */
        } else {
            sys.planets[i].pl_radius = powf(sys.planets[i].pl_mass, 0.5f);
        }
    }
    
    /* Sort by semi-major axis */
    for (int i = 0; i < sys.n_planets-1; i++) {
        for (int j = i+1; j < sys.n_planets; j++) {
            if (sys.planets[j].pl_orbsmax < sys.planets[i].pl_orbsmax) {
                lg_exoplanet_record_t tmp = sys.planets[i];
                sys.planets[i] = sys.planets[j];
                sys.planets[j] = tmp;
            }
        }
    }
    
    /* Apply resonant chain if flagged */
    if (sys.is_resonant_chain) {
        lg_resonant_chain(sys.planets, sys.n_planets);
    }
    
    /* Stability enforcement: if Hill unstable, scale up separations */
    if (!lg_system_hill_stable(sys.planets, sys.n_planets, 10.0f)) {
        for (int i = 0; i < sys.n_planets; i++) {
            sys.planets[i].pl_orbsmax *= 1.5f; /* Spread them out */
            sys.planets[i].pl_orbper *= powf(1.5f, 1.5f); /* P ~ a^{3/2} */
        }
    }
    
    return sys;
}

/*============================================================================
 * Export to Lagrange Particle System
 * Convert generated catalog to simulation initial conditions
 *===========================================================================*/

static inline void lg_exo_to_particle_system(const lg_exosystem_t* sys,
                                             lg_particle_system_t* ps,
                                             lg_pcg32_t* rng) {
    /* Total particles: 1 star + n planets */
    int n = 1 + sys->n_planets;
    ps->n = n;
    
    /* Allocate SoA (caller must manage memory) */
    /* ... assuming ps->pos.x etc. already allocated ... */
    
    /* Star at index 0 */
    ps->mass[0] = sys->star.st_mass * 1.989e30f / 1e24f; /* Convert to sim units if needed */
    ps->pos.x[0] = 0.0f; ps->pos.y[0] = 0.0f; ps->pos.z[0] = 0.0f;
    ps->vel.x[0] = 0.0f; ps->vel.y[0] = 0.0f; ps->vel.z[0] = 0.0f;
    ps->time[0].dt = 0.001f; /* Star gets smallest dt or special handling */
    
    /* Planets on circular-ish orbits in XY plane with slight inclinations */
    for (int i = 0; i < sys->n_planets; i++) {
        int idx = i + 1;
        const lg_exoplanet_record_t* pl = &sys->planets[i];
        
        /* Position at periastron or random true anomaly */
        float a = pl->pl_orbsmax * 1.496e11f; /* AU to meters */
        float e = pl->pl_orbeccen;
        float inc = pl->pl_orbincl * M_PI / 180.0f;
        float lper = pl->pl_orblper * M_PI / 180.0f;
        
        /* Start at periastron for simplicity */
        float r = a * (1.0f - e);
        float x = r * cosf(lper);
        float y = r * sinf(lper);
        
        /* Rotate by inclination around x-axis */
        float y_rot = y * cosf(inc);
        float z_rot = y * sinf(inc);
        
        ps->pos.x[idx] = x;
        ps->pos.y[idx] = y_rot;
        ps->pos.z[idx] = z_rot;
        
        /* Velocity (perpendicular to radius in orbital plane) */
        float mu = sys->star.st_mass * 1.989e30f * 6.67430e-11f; /* GM */
        float h = sqrtf(mu * a * (1.0f - e*e)); /* Specific angular momentum */
        float vx = -h / r * sinf(lper);
        float vy = h / r * cosf(lper) * cosf(inc);
        float vz = h / r * cosf(lper) * sinf(inc);
        
        ps->vel.x[idx] = vx;
        ps->vel.y[idx] = vy;
        ps->vel.z[idx] = vz;
        
        ps->mass[idx] = pl->pl_mass * 1.898e27f; /* Jupiter masses to kg */
        
        /* Time step based on orbital period (fraction of innermost) */
        ps->time[idx].dt = pl->pl_orbper * 86400.0f / 1000.0f; /* 1/1000 of period in seconds */
    }
}

/* Cleanup */
static inline void lg_exosystem_destroy(lg_exosystem_t* sys) {
    free(sys->planets);
    sys->planets = NULL;
    sys->n_planets = 0;
}

#ifdef __cplusplus
}
#endif

#endif /* LAGRANGE_EXOPLANET_H */
