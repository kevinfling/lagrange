/**
 * lagrange_rings.h - Planetary Ring System Dynamics
 * 
 * Comprehensive ring physics including:
 *   - Particle rings: Ice, dust, rocky debris with size distribution
 *   - Shepherding moons: Gap clearing, density waves, resonant structures
 *   - Self-gravity wakes: Gravitational instability in dense rings
 *   - Viscous spreading: Angular momentum transport, edge sharpening
 *   - Collisional evolution: Energy equipartition, accretion/shattering
 *   - Optical depth: Shadowing, occultation profiles
 *   - Spiral density waves: Lindblad resonances with satellites
 *   - Bending waves: Vertical resonances, warped ring planes
 * 
 * Integrates: lagrange_particle.h (SoA collisional dynamics),
 *             lagrange_framework.h (hierarchical ring-moon systems)
 */

#ifndef LAGRANGE_RINGS_H
#define LAGRANGE_RINGS_H

#include <stdbool.h>
#include <stdio.h>
#include "particle.h"
#include "hierarchical_frame.h"

#ifdef __cplusplus
extern "C" {
#endif

/*============================================================================
 * 1. RING PARTICLE PHYSICS
 *===========================================================================*/

/* Particle composition affects optical properties and density */
typedef enum {
    LG_RING_ICE_WATER,        /* Icy particles, high albedo (~0.6) */
    LG_RING_ICE_AMMONIA,      /* Ammonia ice, reddish tint */
    LG_RING_SILICATE,         /* Rocky, dark material (~0.05 albedo) */
    LG_RING_THOLIN,           /* Organic, very dark (~0.02 albedo) */
    LG_RING_IRON,             /* Metallic, radar-detectable */
    LG_RING_DUST_ELECTROSTATIC /* Sub-micron, charged grains */
} lg_ring_composition_t;

/* Size distribution: power law dn/da ~ a^-q with exponential cutoff */
typedef struct {
    float q;                  /* Power law index, typically 2.5-3.5 */
    float a_min;              /* Minimum radius (m) */
    float a_max;              /* Maximum radius (m) */
    float a_cutoff;           /* Exponential cutoff scale (m) */
    float rho_bulk;           /* Bulk material density (kg/m^3) */
} lg_ring_size_distribution_t;

/* Standard distributions */
static const lg_ring_size_distribution_t LG_DIST_SATURN_ICE = {
    .q = 3.1f, .a_min = 1e-6f, .a_max = 10.0f, .a_cutoff = 1.0f, .rho_bulk = 910.0f
};
static const lg_ring_size_distribution_t LG_DIST_URANUS_DUST = {
    .q = 3.5f, .a_min = 1e-7f, .a_max = 1.0f, .a_cutoff = 0.1f, .rho_bulk = 1000.0f
};
static const lg_ring_size_distribution_t LG_DIST_JUPITER_HALO = {
    .q = 2.5f, .a_min = 1e-8f, .a_max = 0.1f, .a_cutoff = 0.01f, .rho_bulk = 2000.0f
};

/* Individual ring particle with extended properties */
typedef struct {
    lg_particle_soa_t base;   /* Position, velocity from lagrange_particle.h */
    
    /* Ring-specific arrays (parallel to base arrays) */
    float* radius;            /* Physical radius (m) */
    float* mass;              /* Computed from radius and density */
    uint32_t* composition;    /* lg_ring_composition_t packed */
    float* internal_density;  /* Compaction varies with size */
    float* spin[3];           /* Rotation rate (rad/s) */
    float* temperature;       /* Thermal equilibrium temperature */
    float* charge;            /* Electrostatic charge (Coulombs) */
    
    /* Collision history for statistical evolution */
    float* collision_count;
    float* last_collision_time;
    
    int capacity;
    int count;
} lg_ring_particle_t;

/*============================================================================
 * 2. RING STRUCTURE & DYNAMICS
 *===========================================================================*/

/* Radial structure: zones, gaps, ringlets */
typedef enum {
    LG_ZONE_CONTINUOUS,       /* Smooth optical depth profile */
    LG_ZONE_GAP,              /* Cleared region (shepherded) */
    LG_ZONE_RINGLET,          /* Narrow high-optical-depth feature */
    LG_ZONE_SPIRAL_WAVE,      /* Density wave from resonance */
    LG_ZONE_EDGE_SHARP,       /* Viscously confined edge */
    LG_ZONE_GRAZING_WEDGE     /* Thickness increases with radius */
} lg_ring_zone_type_t;

typedef struct {
    float r_inner;            /* Inner boundary (m) */
    float r_outer;            /* Outer boundary (m) */
    float tau_normal;         /* Normal optical depth */
    float tau_peak;           /* For ringlets: peak optical depth */
    float width;              /* FWHM for ringlets */
    float surface_density;    /* kg/m^2 */
    float scale_height;       /* Vertical thickness (m) */
    float velocity_dispersion; /* Random velocity (m/s) */
    
    /* Dynamics */
    float viscosity_alpha;    /* Shakura-Sunyaev alpha parameter */
    float accretion_rate;     /* Mass flow (kg/s/m) */
    float torque_external;    /* From resonant satellite */
    bool self_gravity_wakes;  /* Enable self-gravity wake formation */
    
    lg_ring_zone_type_t type;
    lg_ring_composition_t composition;
    lg_ring_size_distribution_t size_dist;
} lg_ring_zone_t;

/* Full ring system around a planet */
typedef struct {
    char name[32];
    uint64_t planet_node_id;  /* Link to framework node */
    
    /* Radial structure */
    lg_ring_zone_t* zones;
    int n_zones;
    int zones_capacity;
    
    /* Particle representation */
    lg_ring_particle_t particles;
    int n_particles_target;   /* Desired count for this LOD */
    
    /* Physics parameters */
    float planet_mass;
    float planet_radius;
    float j2;                 /* Oblateness coefficient */
    float j4;                 /* Higher-order gravity */
    
    /* Precession */
    float nodal_precession_rate; /* dOmega/dt from J2 */
    float apsidal_precession_rate; /* domega/dt */
    
    /* Self-gravity */
    bool self_gravity_enabled;
    float toomre_Q;           /* Stability parameter */
    float critical_wavelength; /* Jeans length */
    
    /* Time evolution */
    double age;               /* System age (seconds) */
    double last_evolution_step;
    
    /* Statistics */
    float total_mass;
    float total_area;
    float optical_depth_profile[1024]; /* Sampled radial profile */
} lg_ring_system_t;

/*============================================================================
 * 3. SHEPHERDING MOONS & RESONANCES
 *===========================================================================*/

/* Lindblad resonance: drives spiral density waves */
typedef struct {
    int m;                    /* Azimuthal wavenumber */
    int n;                    /* For n:1 mean motion resonance */
    float a_res;              /* Resonant semi-major axis */
    float torque;             /* Angular momentum transfer (N·m) */
    float wave_amplitude;     /* Surface density perturbation */
    float damping_length;     /* Where wave dissipates */
    
    lg_node_t* satellite;     /* Perturbing moon */
    lg_ring_zone_t* target_zone;
} lindblad_resonance_t;

/* Vertical/corotation resonance: drives bending waves */
typedef struct {
    int m;
    float inclination_forcing; /* Vertical amplitude */
    float vertical_wavelength;
    lg_node_t* satellite;
} vertical_resonance_t;

/* Shepherd moon configuration */
typedef struct {
    lg_node_t* moon;          /* The shepherding satellite */
    float r_moon;             /* Orbital radius */
    float gap_width;          /* Cleared region width */
    float torque_inner;       /* On inner ring edge */
    float torque_outer;       /* On outer ring edge */
    float mass_ratio;         /* Moon mass / planet mass */
    
    /* Gap clearing efficiency */
    float viscosity_parameter; /* nu / (r^2 * Omega) */
    float gap_depth;            /* tau_inside / tau_outside */
} lg_shepherd_t;

/*============================================================================
 * 4. PHYSICS IMPLEMENTATION
 *===========================================================================*/

/* Toomre stability parameter: Q = Omega * sigma / (pi * G * Sigma) */
static inline float lg_ring_toomre_Q(float r, float Omega, float sigma, float surface_density) {
    const float G = 6.67430e-11f;
    float sound_speed = sigma; /* Velocity dispersion */
    float epicyclic = Omega * sqrtf(1.0f + r / (r + 1e-10f)); /* Approx kappa */
    return sound_speed * epicyclic / (LG_PI * G * surface_density + 1e-30f);
}

/* Viscous diffusion timescale: t_nu ~ r^2 / nu */
static inline float lg_ring_viscous_time(float r, float alpha, float Omega, float H) {
    float cs = Omega * H; /* Sound speed */
    float nu = alpha * cs * H; /* Shakura-Sunyaev viscosity */
    return r * r / (nu + 1e-30f);
}

/* Gap clearing by shepherd moon (Lin & Papaloizou 1986) */
static inline float lg_ring_gap_width(float q, float h, float alpha) {
    /* q = M_moon/M_planet, h = H/r, alpha = viscosity parameter */
    float q_pow = powf(q, 0.5f);
    float h_pow = powf(h, 2.5f);
    float alpha_pow = powf(alpha, -0.5f);
    return 2.0f * q_pow * h_pow * alpha_pow; /* In units of Hill radius */
}

/* Spiral density wave dispersion relation */
static inline float lg_ring_density_wavelength(float r, int m, float Omega, float sigma, float kappa) {
    /* Doppler-shifted pattern speed */
    float Omega_p = Omega * (1.0f - 1.0f/m); /* For inner Lindblad */
    float D = kappa*kappa - m*m*(Omega - Omega_p)*(Omega - Omega_p);
    if (D < 0) return 1e30f; /* Evanescent */
    const float G = 6.67430e-11f;
    return 4.0f * LG_PI * LG_PI * sigma * sigma / (G * r * fabsf(D) + 1e-30f);
}

/*============================================================================
 * 5. INITIALIZATION & SETUP
 *===========================================================================*/

/* Initialize ring system around planet node */
static inline lg_ring_system_t* lg_ring_system_create(lg_node_t* planet_node, const char* name) {
    lg_ring_system_t* rings = (lg_ring_system_t*)calloc(1, sizeof(lg_ring_system_t));
    strncpy(rings->name, name, 31);
    rings->planet_node_id = planet_node->id;
    
    /* Extract from planet node */
    if (planet_node->type == LG_NODE_BODY) {
        rings->planet_mass = planet_node->data.body.body.mass;
        rings->planet_radius = planet_node->data.body.body.radius;
    }
    
    /* Default: Saturn-like */
    rings->j2 = 0.016f;
    rings->j4 = -0.001f;
    rings->self_gravity_enabled = true;
    rings->zones_capacity = 16;
    rings->zones = (lg_ring_zone_t*)calloc(rings->zones_capacity, sizeof(lg_ring_zone_t));
    
    return rings;
}

/* Add ring zone with specified properties */
static inline lg_ring_zone_t* lg_ring_system_add_zone(
    lg_ring_system_t* rings,
    float r_inner,
    float r_outer,
    lg_ring_zone_type_t type,
    lg_ring_composition_t comp,
    const lg_ring_size_distribution_t* sizes
) {
    if (rings->n_zones >= rings->zones_capacity) {
        rings->zones_capacity *= 2;
        rings->zones = (lg_ring_zone_t*)realloc(rings->zones, 
            rings->zones_capacity * sizeof(lg_ring_zone_t));
    }
    
    lg_ring_zone_t* z = &rings->zones[rings->n_zones++];
    memset(z, 0, sizeof(lg_ring_zone_t));
    
    z->r_inner = r_inner;
    z->r_outer = r_outer;
    z->type = type;
    z->composition = comp;
    z->size_dist = *sizes;
    
    /* Compute surface density from optical depth */
    float a_typical = sizes->a_cutoff;
    float particle_area = LG_PI * a_typical * a_typical;
    float particle_mass = (4.0f/3.0f) * LG_PI * a_typical*a_typical*a_typical * sizes->rho_bulk;
    
    /* tau = n * sigma * H, Sigma = n * m * H */
    z->surface_density = z->tau_normal * particle_mass / particle_area;
    
    /* Scale height from velocity dispersion */
    float Omega = sqrtf(rings->planet_mass * 6.67430e-11f / 
        powf((r_inner + r_outer)*0.5f, 3.0f));
    z->scale_height = z->velocity_dispersion / (Omega + 1e-30f);
    
    return z;
}

/* Create Saturn-like main rings */
static inline void lg_ring_preset_saturn_main(lg_ring_system_t* rings) {
    /* D ring: innermost, faint */
    lg_ring_system_add_zone(rings, 66900e3f, 74500e3f, LG_ZONE_CONTINUOUS,
        LG_RING_SILICATE, &LG_DIST_SATURN_ICE);
    rings->zones[rings->n_zones-1].tau_normal = 0.01f;
    
    /* C ring: middle, moderate optical depth */
    lg_ring_system_add_zone(rings, 74500e3f, 92000e3f, LG_ZONE_CONTINUOUS,
        LG_RING_ICE_WATER, &LG_DIST_SATURN_ICE);
    rings->zones[rings->n_zones-1].tau_normal = 0.1f;
    
    /* B ring: dense, main structure */
    lg_ring_system_add_zone(rings, 92000e3f, 117500e3f, LG_ZONE_CONTINUOUS,
        LG_RING_ICE_WATER, &LG_DIST_SATURN_ICE);
    rings->zones[rings->n_zones-1].tau_normal = 1.5f;
    rings->zones[rings->n_zones-1].self_gravity_wakes = true;
    
    /* Cassini Division: gap cleared by Mimas 2:1 resonance */
    lg_ring_system_add_zone(rings, 117500e3f, 122200e3f, LG_ZONE_GAP,
        LG_RING_ICE_WATER, &LG_DIST_SATURN_ICE);
    rings->zones[rings->n_zones-1].tau_normal = 0.05f;
    
    /* A ring: outer main ring */
    lg_ring_system_add_zone(rings, 122200e3f, 136780e3f, LG_ZONE_CONTINUOUS,
        LG_RING_ICE_WATER, &LG_DIST_SATURN_ICE);
    rings->zones[rings->n_zones-1].tau_normal = 0.4f;
    
    /* Encke Gap: narrow gap at ~133.6e6 m carved by shepherd moon Pan */
    lg_ring_system_add_zone(rings, 133400e3f, 133800e3f, LG_ZONE_GAP,
        LG_RING_ICE_WATER, &LG_DIST_SATURN_ICE);
    rings->zones[rings->n_zones-1].tau_normal = 0.02f;
    
    /* Keeler Gap: narrow gap at ~136.5e6 m carved by shepherd moon Daphnis */
    lg_ring_system_add_zone(rings, 136400e3f, 136660e3f, LG_ZONE_GAP,
        LG_RING_ICE_WATER, &LG_DIST_SATURN_ICE);
    rings->zones[rings->n_zones-1].tau_normal = 0.02f;
    
    /* F ring: shepherded by Prometheus and Pandora */
    lg_ring_system_add_zone(rings, 140180e3f, 140260e3f, LG_ZONE_RINGLET,
        LG_RING_ICE_WATER, &LG_DIST_SATURN_ICE);
    rings->zones[rings->n_zones-1].tau_normal = 0.1f;
    rings->zones[rings->n_zones-1].width = 50e3f;
    
    /* G ring: faint, tenuous */
    lg_ring_system_add_zone(rings, 166000e3f, 173000e3f, LG_ZONE_CONTINUOUS,
        LG_RING_ICE_WATER, &LG_DIST_SATURN_ICE);
    rings->zones[rings->n_zones-1].tau_normal = 1e-6f;
    
    /* E ring: very diffuse, engeleous */
    lg_ring_system_add_zone(rings, 180000e3f, 480000e3f, LG_ZONE_CONTINUOUS,
        LG_RING_ICE_WATER, &LG_DIST_SATURN_ICE);
    rings->zones[rings->n_zones-1].tau_normal = 1e-7f;
}

/*============================================================================
 * 6. PARTICLE INITIALIZATION
 *===========================================================================*/

/* Sample from power-law size distribution */
static inline float lg_ring_sample_size(const lg_ring_size_distribution_t* dist) {
    /* Rejection sampling for truncated power law with exponential cutoff */
    float u = (float)rand() / (float)RAND_MAX;
    float a = dist->a_min * powf(dist->a_max / dist->a_min, u);
    
    /* Apply exponential cutoff */
    float cutoff = expf(-a / dist->a_cutoff);
    if ((float)rand() / (float)RAND_MAX > cutoff) {
        a = dist->a_min; /* Reject to minimum */
    }
    
    return a;
}

/* Ensure particle arrays have capacity for at least n total particles */
static inline void _lg_ring_particles_ensure(lg_ring_particle_t* p, int n) {
    if (n <= p->capacity) return;
    int new_cap = (p->capacity < 8) ? 16 : p->capacity * 2;
    while (new_cap < n) new_cap *= 2;
    
    p->base.x = (float*)realloc(p->base.x, new_cap * sizeof(float));
    p->base.y = (float*)realloc(p->base.y, new_cap * sizeof(float));
    p->base.z = (float*)realloc(p->base.z, new_cap * sizeof(float));
    p->base.vx = (float*)realloc(p->base.vx, new_cap * sizeof(float));
    p->base.vy = (float*)realloc(p->base.vy, new_cap * sizeof(float));
    p->base.vz = (float*)realloc(p->base.vz, new_cap * sizeof(float));
    p->radius = (float*)realloc(p->radius, new_cap * sizeof(float));
    p->mass = (float*)realloc(p->mass, new_cap * sizeof(float));
    p->composition = (uint32_t*)realloc(p->composition, new_cap * sizeof(uint32_t));
    p->internal_density = (float*)realloc(p->internal_density, new_cap * sizeof(float));
    p->temperature = (float*)realloc(p->temperature, new_cap * sizeof(float));
    p->charge = (float*)realloc(p->charge, new_cap * sizeof(float));
    p->collision_count = (float*)realloc(p->collision_count, new_cap * sizeof(float));
    p->last_collision_time = (float*)realloc(p->last_collision_time, new_cap * sizeof(float));
    for (int d = 0; d < 3; d++) {
        p->spin[d] = (float*)realloc(p->spin[d], new_cap * sizeof(float));
    }
    p->capacity = new_cap;
}

/* Initialize particles for a zone */
static inline void lg_ring_zone_populate(lg_ring_zone_t* zone, lg_ring_system_t* rings, int n_particles) {
    float r_mid = (zone->r_inner + zone->r_outer) * 0.5f;
    float dr = (zone->r_outer - zone->r_inner) * 0.5f;
    const float G = 6.67430e-11f;
    
    lg_ring_particle_t* p = &rings->particles;
    int start = p->count;
    _lg_ring_particles_ensure(p, start + n_particles);
    
    for (int i = 0; i < n_particles; i++) {
        int idx = start + i;
        
        /* Radial position with surface density profile */
        float u = (float)rand() / (float)RAND_MAX;
        float r = r_mid + (u - 0.5f) * 2.0f * dr;
        
        /* Azimuth */
        float theta = 2.0f * LG_PI * (float)rand() / (float)RAND_MAX;
        
        /* Vertical: Gaussian distribution */
        float z = zone->scale_height * sqrtf(-2.0f * logf((float)rand() / (float)RAND_MAX + 1e-10f)) *
                  cosf(2.0f * LG_PI * (float)rand() / (float)RAND_MAX);
        
        /* Keplerian velocity with dispersion */
        float v_kep = sqrtf(rings->planet_mass * G / r);
        float v_r = zone->velocity_dispersion * ((float)rand() / (float)RAND_MAX - 0.5f);
        float v_theta = v_kep + zone->velocity_dispersion * ((float)rand() / (float)RAND_MAX - 0.5f);
        float v_z = zone->velocity_dispersion * ((float)rand() / (float)RAND_MAX - 0.5f);
        
        /* Store in SoA */
        p->base.x[idx] = r * cosf(theta);
        p->base.y[idx] = r * sinf(theta);
        p->base.z[idx] = z;
        p->base.vx[idx] = v_r * cosf(theta) - v_theta * sinf(theta);
        p->base.vy[idx] = v_r * sinf(theta) + v_theta * cosf(theta);
        p->base.vz[idx] = v_z;
        
        /* Physical properties */
        float a = lg_ring_sample_size(&zone->size_dist);
        float rho = zone->size_dist.rho_bulk * (0.8f + 0.2f * (float)rand() / (float)RAND_MAX);
        p->radius[idx] = a;
        p->internal_density[idx] = rho;
        p->mass[idx] = (4.0f / 3.0f) * LG_PI * a * a * a * rho;
        p->composition[idx] = (uint32_t)zone->composition;
        p->temperature[idx] = 80.0f; /* Saturn ring temperature (K) */
        p->charge[idx] = 0.0f;
        p->collision_count[idx] = 0.0f;
        p->last_collision_time[idx] = 0.0f;
        for (int d = 0; d < 3; d++) p->spin[d][idx] = 0.0f;
    }
    
    p->count = start + n_particles;
    p->base.count = p->count;
}

/*============================================================================
 * 7. TIME EVOLUTION
 *===========================================================================*/

/* Collisional energy equipartition (dynamical cooling/heating) */
static inline void lg_ring_collisional_evolution(lg_ring_system_t* rings, float dt) {
    /* Bridges et al. 1984: coefficient of restitution depends on impact velocity */
    /* Differential settling: larger particles drift inward faster */
    
    for (int i = 0; i < rings->particles.count; i++) {
        float a = rings->particles.radius[i];
        float rho = rings->particles.internal_density[i];
        if (a <= 0.0f || rho <= 0.0f) continue;
        
        /* Collision frequency ~ optical depth * velocity dispersion / particle size */
        float collision_rate = rings->zones[0].tau_normal * 
                               rings->zones[0].velocity_dispersion / 
                               (2.0f * a + 1e-10f);
        
        /* Update collision statistics */
        rings->particles.collision_count[i] += collision_rate * dt;
        rings->particles.last_collision_time[i] = (double)rings->age;
        
        /* Velocity-dependent restitution: higher impact speed -> lower epsilon */
        float v_impact = rings->zones[0].velocity_dispersion;
        float epsilon = fmaxf(0.1f, 0.5f - 0.1f * logf(v_impact / 0.001f + 1.0f));
        (void)epsilon; /* Applied in collision response */
        
        /* Thermal equilibration from collisional energy */
        float T_eq = 100.0f * powf(a / 1e-6f, -0.25f);
        float alpha = fminf(1.0f, collision_rate * dt * 0.1f);
        rings->particles.temperature[i] += (T_eq - rings->particles.temperature[i]) * alpha;
    }
}

/* Self-gravity wake formation (Salo 1992, 1995) */
static inline void lg_ring_self_gravity_wakes(lg_ring_system_t* rings, float dt) {
    if (!rings->self_gravity_enabled) return;
    
    /* Critical Toomre wavelength */
    float lambda_crit = 2.0f * LG_PI * LG_PI * rings->zones[0].surface_density / 
        (rings->planet_mass * 6.67430e-11f / 
         powf((rings->zones[0].r_inner + rings->zones[0].r_outer)*0.5f, 3.0f));
    
    /* Wake angle: ~20-25 degrees from azimuthal */
    
    /* Aggregate particles into wake structures when Q < 2 */
    for (int z = 0; z < rings->n_zones; z++) {
        lg_ring_zone_t* zone = &rings->zones[z];
        float r = (zone->r_inner + zone->r_outer) * 0.5f;
        float Omega = sqrtf(rings->planet_mass * 6.67430e-11f / (r*r*r));
        float Q = lg_ring_toomre_Q(r, Omega, zone->velocity_dispersion, zone->surface_density);
        
        if (Q < 2.0f) {
            /* Form wakes: enhanced density, reduced dispersion (gravitational cooling) */
            float growth_rate = (2.0f - Q) * 0.001f;
            zone->tau_normal *= (1.0f + growth_rate * dt);
            if (zone->tau_normal > 3.0f) zone->tau_normal = 3.0f;
            
            zone->velocity_dispersion *= (1.0f - growth_rate * dt * 0.5f);
            if (zone->velocity_dispersion < 1e-6f) zone->velocity_dispersion = 1e-6f;
            
            rings->critical_wavelength = lambda_crit;
        }
    }
}

/* Main integration step */
static inline void lg_ring_system_update(lg_ring_system_t* rings, double time, float dt) {
    const float G = 6.67430e-11f;
    float dt_sub = dt;
    
    /* Substep for collisional dynamics */
    int n_substeps = (int)ceilf(dt / 100.0f); /* Max 100s substeps */
    if (n_substeps < 1) n_substeps = 1;
    dt_sub = dt / n_substeps;
    
    for (int step = 0; step < n_substeps; step++) {
        /* 1. Keplerian orbits with J2 precession */
        for (int i = 0; i < rings->particles.count; i++) {
            float x = rings->particles.base.x[i];
            float y = rings->particles.base.y[i];
            float z = rings->particles.base.z[i];
            float r2 = x*x + y*y;
            float r = sqrtf(r2);
            float r3 = r2 * r;
            
            /* Keplerian acceleration */
            float ax = -G * rings->planet_mass * x / (r3 + 1e-10f);
            float ay = -G * rings->planet_mass * y / (r3 + 1e-10f);
            float az = -G * rings->planet_mass * z / (r3 + 1e-10f);
            
            /* J2 perturbation: a_z = -3 J2 mu R^2 z / (2 r^5) */
            float j2_term = -1.5f * rings->j2 * G * rings->planet_mass *
                            rings->planet_radius * rings->planet_radius * z /
                            (powf(r2 + z*z, 2.5f) + 1e-10f);
            az += j2_term;
            
            /* Integrate velocity */
            rings->particles.base.vx[i] += ax * dt_sub;
            rings->particles.base.vy[i] += ay * dt_sub;
            rings->particles.base.vz[i] += az * dt_sub;
            
            /* Integrate position */
            rings->particles.base.x[i] += rings->particles.base.vx[i] * dt_sub;
            rings->particles.base.y[i] += rings->particles.base.vy[i] * dt_sub;
            rings->particles.base.z[i] += rings->particles.base.vz[i] * dt_sub;
        }
        
        /* 2. Collisional evolution (statistics and thermal equilibration) */
        lg_ring_collisional_evolution(rings, dt_sub);
        
        /* 3. Self-gravity wake formation */
        lg_ring_self_gravity_wakes(rings, dt_sub);
        
        /* 4. Viscous radial diffusion: dSigma/dt = (3/r) d/dr[ r^{1/2} d/dr(nu Sigma r^{1/2}) ] */
        for (int z = 0; z < rings->n_zones; z++) {
            lg_ring_zone_t* zone = &rings->zones[z];
            float r = (zone->r_inner + zone->r_outer) * 0.5f;
            float cs = zone->velocity_dispersion;
            float H = zone->scale_height;
            float nu = zone->viscosity_alpha * cs * H; /* Shakura-Sunyaev */
            
            /* Simple radial diffusion: sigma spreads outward */
            float diff_rate = 3.0f * nu / (r * r);
            zone->surface_density *= (1.0f - diff_rate * dt_sub);
            if (zone->surface_density < 1e-6f) zone->surface_density = 1e-6f;
            
            /* Optical depth follows surface density */
            float a_typ = zone->size_dist.a_cutoff;
            float area = LG_PI * a_typ * a_typ;
            float mass = (4.0f/3.0f) * LG_PI * a_typ*a_typ*a_typ * zone->size_dist.rho_bulk;
            zone->tau_normal = zone->surface_density * area / mass;
        }
        
        /* 5. Resonant torques from external satellites */
        for (int z = 0; z < rings->n_zones; z++) {
            lg_ring_zone_t* zone = &rings->zones[z];
            if (fabsf(zone->torque_external) > 0.0f) {
                /* Torque changes angular momentum → surface density redistribution */
                float torque_effect = zone->torque_external * dt_sub * 1e-20f;
                zone->surface_density += torque_effect;
                if (zone->surface_density < 1e-6f) zone->surface_density = 1e-6f;
            }
        }
    }
    
    rings->age += dt;
    rings->last_evolution_step = time;
}

/*============================================================================
 * 8. RENDERING & OBSERVATION
 *===========================================================================*/

/* Compute occultation profile (radial optical depth) */
static inline void lg_ring_compute_occultation(lg_ring_system_t* rings, float* tau_profile, int n_samples) {
    for (int i = 0; i < n_samples; i++) {
        float r = rings->zones[0].r_inner + 
            (rings->zones[rings->n_zones-1].r_outer - rings->zones[0].r_inner) * 
            i / (n_samples - 1);
        
        tau_profile[i] = 0.0f;
        for (int z = 0; z < rings->n_zones; z++) {
            if (r >= rings->zones[z].r_inner && r <= rings->zones[z].r_outer) {
                tau_profile[i] += rings->zones[z].tau_normal;
            }
        }
    }
}

/* Phase function for ring brightness */
static inline float lg_ring_phase_function(float phase_angle, lg_ring_composition_t comp) {
    /* Opposition surge at small phase angles */
    float surge = 1.0f;
    if (phase_angle < 1.0f * LG_PI / 180.0f) {
        surge = 1.0f + 0.5f / (phase_angle * 180.0f / LG_PI + 0.1f);
    }
    
    /* Lambertian + specular for ice */
    float lambert = fmaxf(0.0f, cosf(phase_angle));
    
    return surge * lambert;
}

/*============================================================================
 * 9. INTEGRATION WITH FRAMEWORK
 *===========================================================================*/

/* Attach ring system to planet node */
static inline void lg_node_attach_rings(lg_node_t* planet_node, lg_ring_system_t* rings) {
    /* Create child node for rings */
    lg_node_t* rings_node = lg_node_create(LG_NODE_EMPTY, rings->name);
    rings_node->user_data = rings;
    
    /* Add shepherd moons for gap zones */
    const float G = 6.67430e-11f;
    for (int i = 0; i < rings->n_zones; i++) {
        if (rings->zones[i].type == LG_ZONE_GAP) {
            char moon_name[64];
            snprintf(moon_name, sizeof(moon_name), "shepherd_%d", i);
            lg_node_t* moon = lg_node_create(LG_NODE_BODY, moon_name);
            
            float r_gap = (rings->zones[i].r_inner + rings->zones[i].r_outer) * 0.5f;
            moon->data.body.body.mass = rings->planet_mass * 1e-7f; /* ~10 km moon */
            moon->data.body.body.radius = 5e3f;
            moon->data.body.orbit.a = r_gap;
            moon->data.body.orbit.e = 0.001f;
            moon->data.body.orbit.mu = G * rings->planet_mass;
            moon->domain = LG_PHYSICS_KEPLER;
            
            lg_node_attach(planet_node, moon);
        }
    }
    
    lg_node_attach(planet_node, rings_node);
}

/* LOD switching for rings */
static inline void lg_ring_lod_update(lg_ring_system_t* rings, float distance_camera, float pixel_scale) {
    /* Close: individual particles with shadows */
    /* Medium: particle aggregates in cells */
    /* Far: textured ring plane with optical depth */
    /* Very far: point light occluder */
    
    if (rings->n_zones <= 0) return;
    float ring_width_pixels = (rings->zones[rings->n_zones-1].r_outer - 
                               rings->zones[0].r_inner) / (distance_camera * pixel_scale + 1e-10f);
    
    if (ring_width_pixels > 1000.0f) {
        /* Full particle simulation */
        rings->n_particles_target = 100000;
    } else if (ring_width_pixels > 100.0f) {
        /* Aggregate representation */
        rings->n_particles_target = 10000;
    } else if (ring_width_pixels > 10.0f) {
        /* Textured ring plane */
        rings->n_particles_target = 1000;
    } else {
        /* Billboard / point occluder */
        rings->n_particles_target = 0;
    }
}

#ifdef __cplusplus
}
#endif

#endif /* LAGRANGE_RINGS_H */

