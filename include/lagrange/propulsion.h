/**
 * lagrange_propulsion.h - Comprehensive Astrodynamics: Propulsion & Navigation
 * 
 * Complete mission planning and execution toolkit:
 *   - Thrust curves: Variable Isp, staging events, throttle profiles
 *   - Gravity assists: Powered/unpowered flybys, resonant returns
 *   - Porkchop plots: Launch window analysis with C3/v_inf contours
 *   - Lambert solver: Universal variable algorithm (Izzo/Battin)
 *   - Low-thrust: Equinoctial elements for ion drives & solar sails
 *   - Aerobraking: Atmospheric entry with heating constraints
 *   - Tether physics: Skyhook dynamics, momentum exchange
 * 
 * Integrates: lagrange_patched_conic.h, lagrange_floating_origin.h,
 *             lagrange_particle.h (SoA thrusting bodies)
 */

#ifndef LAGRANGE_PROPULSION_H
#define LAGRANGE_PROPULSION_H

#include "math.h"
#include "gravity.h"
#include "patched_conic.h"
#include "floating_origin.h"
#include "body.h"
#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

/*============================================================================
 * 1. THRUST CURVES & STAGING
 * Realistic rocket dynamics with mass flow and variable efficiency
 *===========================================================================*/

typedef enum {
    LG_PROP_CONSTANT_Isp,     /* Fixed specific impulse */
    LG_PROP_VARIABLE_Isp,     /* Throttle-dependent Isp (e.g., nuclear thermal) */
    LG_PROP_SOLAR_ELECTRIC,   /* Power-limited, Isp~throttle^-1/2 */
    LG_PROP_SAIL,             /* Solar sail (no mass flow, force ~ 1/r^2) */
    LG_PROP_AERO_BRAKE      /* Negative thrust (drag) - see section 6 */
} lg_propulsion_type_t;

typedef struct {
    float thrust_N;           /* Current thrust magnitude */
    float Isp_s;              /* Specific impulse (seconds) */
    float mdot;               /* Mass flow rate (kg/s), derived */
    float throttle;           /* 0.0 to 1.0 */
    
    /* Staging event */
    float dry_mass;           /* Stage dry mass (kg) */
    float prop_mass;          /* Current propellant remaining */
    float stage_time;         /* Time since ignition */
    float max_burn_time;      /* Stage lifetime at full throttle */
    
    /* For variable Isp models */
    float Isp_min, Isp_max;   /* Bounds for variable Isp engines */
    float power_W;            /* Available electrical power (for electric) */
    
    /* Sail parameters (if type == LG_PROP_SAIL) */
    float sail_area_m2;       /* Area */
    float sail_efficiency;    /* 0.0-1.0 (reflectivity factor) */
    
    lg_propulsion_type_t type;
} lg_stage_t;

/* Compute mass flow from thrust and Isp: mdot = thrust / (g0 * Isp) */
static inline void lg_stage_update_mdot(lg_stage_t* stage) {
    const float g0 = 9.80665f;
    stage->mdot = (stage->throttle * stage->thrust_N) / (g0 * stage->Isp_s);
}

/* Solar electric propulsion: optimize Isp for max thrust given power P */
static inline void lg_stage_sepp_optimize(lg_stage_t* stage) {
    /* P = thrust * Isp * g0 / (2 * eta), optimal at specific throttle */
    if (stage->type == LG_PROP_SOLAR_ELECTRIC && stage->power_W > 0) {
        /* For fixed power, T = 2*eta*P/(Isp*g0), maximize T at min Isp? */
        /* Actually T = mdot * Isp * g0, mdot = P/(Isp*g0*eta) for ion */
        /* So T = P*eta/(Isp*g0), inverse relationship - lower Isp = higher thrust */
        float eta = 0.6f; /* thruster efficiency */
        stage->thrust_N = 2.0f * eta * stage->power_W / (stage->Isp_s * 9.80665f);
        lg_stage_update_mdot(stage);
    }
}

/* Staging event: drop dry mass, switch to next stage config */
static inline bool lg_stage_burn(lg_stage_t* stage, float dt) {
    float dm = stage->mdot * dt;
    if (dm > stage->prop_mass) {
        dm = stage->prop_mass;
        stage->prop_mass = 0.0f;
        return false; /* Stage burnout */
    }
    stage->prop_mass -= dm;
    stage->stage_time += dt;
    return true;
}

/* Thrust vector in inertial frame (body frame rotated by attitude) */
typedef struct {
    lg_vec3_t direction;      /* Unit vector in body frame */
    lg_vec3_t thrust_vector;  /* Computed inertial thrust = T * R_body * direction */
    lg_stage_t stage;
} lg_thruster_t;

/* Note: Lambert problem solver is now in lambert.h
 * Types: lg_lambert_problem_t, lg_lambert_solution_t  
 * Function: lg_lambert_solve()
 */

/*============================================================================
 * 3. GRAVITY ASSISTS (Powered and Unpowered)
 *===========================================================================*/

typedef struct {
    lg_vec3_t v_inf_in;       /* Incoming hyperbolic excess */
    lg_vec3_t v_inf_out;      /* Outgoing (targeted) */
    float r_periapsis;        /* Closest approach (planet radii) */
    float turn_angle;         /* Deflection from hyperbola */
    float delta_v;            /* Powered flyby delta-v at periapsis */
    float planet_mu;          /* Planet gravitational parameter */
    float planet_radius;      /* For periapsis altitude check */
    float bending_angle;      /* Total deflection: pi - 2*arcsin(1/e) */
} lg_gravity_assist_t;

/* Compute gravity assist deflection from hyperbolic excess velocities */
static inline lg_gravity_assist_t lg_gravity_assist_calculate(
    const lg_vec3_t* v_inf_in,
    const lg_vec3_t* v_inf_out,
    float planet_mu,
    float planet_radius,
    float r_periapsis_factor  /* Multiple of planet radius, e.g., 1.05 for 5% altitude */
) {
    lg_gravity_assist_t ga = {0};
    ga.v_inf_in = *v_inf_in;
    ga.v_inf_out = *v_inf_out;
    ga.planet_mu = planet_mu;
    ga.planet_radius = planet_radius;
    ga.r_periapsis = r_periapsis_factor * planet_radius;
    
    float v_inf_mag = lg_vec3_len(*v_inf_in);
    float ecc = 1.0f + ga.r_periapsis * v_inf_mag * v_inf_mag / planet_mu;
    ga.bending_angle = M_PI - 2.0f * asinf(1.0f / ecc);
    
    /* Required turn angle vs available */
    float cos_turn = lg_vec3_dot(lg_vec3_norm(*v_inf_in), lg_vec3_norm(*v_inf_out));
    ga.turn_angle = acosf(fmaxf(-1.0f, fminf(1.0f, cos_turn)));
    
    /* Powered flyby delta-v if natural deflection insufficient */
    if (ga.turn_angle > ga.bending_angle) {
        /* Need powered assist - compute required delta-v at periapsis */
        float v_peri = sqrtf(v_inf_mag * v_inf_mag + 2.0f * planet_mu / ga.r_periapsis);
        /* Simplified: delta-v to rotate v_inf vector */
        float delta_turn = ga.turn_angle - ga.bending_angle;
        ga.delta_v = 2.0f * v_peri * sinf(delta_turn * 0.5f);
    } else {
        ga.delta_v = 0.0f; /* Unpowered assist sufficient */
    }
    
    return ga;
}

/* Resonant return: compute next encounter after N orbits */
typedef struct {
    float period_ratio;       /* T_spacecraft / T_planet */
    int orbits_planet;        /* Planet orbits between encounters */
    int orbits_sc;            /* Spacecraft orbits between encounters */
    float phase_change;       /* Change in heliocentric phase angle */
    float synodic_period;     /* Time between encounters */
} lg_resonant_return_t;

static inline lg_resonant_return_t lg_resonant_return(
    float a_transfer,         /* Semi-major axis of transfer orbit (AU or km) */
    float a_planet,         /* Planet semi-major axis */
    float mu_central        /* Central body GM */
) {
    lg_resonant_return_t rr = {0};
    
    float T_sc = 2.0f * M_PI * sqrtf(a_transfer * a_transfer * a_transfer / mu_central);
    float T_planet = 2.0f * M_PI * sqrtf(a_planet * a_planet * a_planet / mu_central);
    
    rr.period_ratio = T_sc / T_planet;
    
    /* Find integer ratio approximation (simplified) */
    for (int n = 1; n <= 10; n++) {
        for (int m = 1; m <= 10; m++) {
            float ratio = (float)n / (float)m;
            if (fabsf(ratio - rr.period_ratio) < 0.01f) {
                rr.orbits_planet = m;
                rr.orbits_sc = n;
                rr.synodic_period = T_planet * m;
                rr.phase_change = 2.0f * M_PI * ((float)n / (float)m - 1.0f);
                return rr;
            }
        }
    }
    
    return rr;
}

/*============================================================================
 * 4. PORKCHOP PLOTS - Launch Window Analysis
 *===========================================================================*/

typedef struct {
    float launch_mjd;         /* Launch date (Modified Julian Date) */
    float arrival_mjd;      /* Arrival date */
    float c3;               /* Characteristic energy (km^2/s^2) */
    float v_inf;            /* Hyperbolic excess at target (km/s) */
    float tof;              /* Time of flight (days) */
    float declination;      /* Launch asymptote declination */
    float dv_total;         /* Total mission delta-v budget */
    int flag;               /* 0=valid, 1=no solution, 2=hyperbolic excess too high */
} lg_porkchop_point_t;

typedef struct {
    lg_porkchop_point_t* grid;     /* 2D grid of launch/arrival combinations */
    int n_launch;                  /* Grid dimension for launch dates */
    int n_arrival;                 /* Grid dimension for arrival dates */
    float c3_min, c3_max;          /* Characteristic energy bounds */
    float v_inf_min, v_inf_max;    /* Arrival excess velocity bounds */
    float tof_min, tof_max;        /* Time of flight bounds */
    float best_c3;                 /* Optimal launch window */
    float best_launch_mjd;
    float best_arrival_mjd;
} lg_porkchop_plot_t;

/* Compute C3 (characteristic energy) for interplanetary transfer */
static inline float lg_characteristic_energy(const lg_vec3_t* v_depart, float v_esc) {
    float v = lg_vec3_len(*v_depart);
    return v * v - v_esc * v_esc; /* C3 = v_inf^2 */
}

/* Generate porkchop plot grid */
static inline void lg_porkchop_generate(
    lg_porkchop_plot_t* plot,
    const lg_orbit_t* departure_orbit,
    const lg_orbit_t* arrival_orbit,
    float mu_central,
    float launch_mjd_start,
    float launch_mjd_end,
    int n_launch,
    float tof_min_days,
    float tof_max_days,
    int n_arrival
) {
    plot->n_launch = n_launch;
    plot->n_arrival = n_arrival;
    plot->grid = (lg_porkchop_point_t*)malloc(n_launch * n_arrival * sizeof(lg_porkchop_point_t));
    
    float d_launch = (launch_mjd_end - launch_mjd_start) / (n_launch - 1);
    float d_tof = (tof_max_days - tof_min_days) / (n_arrival - 1);
    
    plot->c3_min = 1e20f;
    plot->c3_max = 0.0f;
    plot->v_inf_min = 1e20f;
    plot->v_inf_max = 0.0f;
    
    for (int i = 0; i < n_launch; i++) {
        float launch_mjd = launch_mjd_start + i * d_launch;
        
        for (int j = 0; j < n_arrival; j++) {
            float tof_days = tof_min_days + j * d_tof;
            float arrival_mjd = launch_mjd + tof_days;
            float tof_seconds = tof_days * 86400.0f;
            
            lg_porkchop_point_t* pt = &plot->grid[i * n_arrival + j];
            pt->launch_mjd = launch_mjd;
            pt->arrival_mjd = arrival_mjd;
            pt->tof = tof_days;
            
            /* Get planet positions at launch and arrival */
            float M1 = departure_orbit->M0 + departure_orbit->n * (launch_mjd - departure_orbit->epoch) * 86400.0f;
            float M2 = arrival_orbit->M0 + arrival_orbit->n * (arrival_mjd - arrival_orbit->epoch) * 86400.0f;
            
            lg_vec3_t r1 = lg_orbit_position_from_mean_anomaly(departure_orbit, M1);
            lg_vec3_t r2 = lg_orbit_position_from_mean_anomaly(arrival_orbit, M2);
            
            /* Solve Lambert problem */
            lg_lambert_problem_t prob = {
                .r1 = r1, .r2 = r2, .dt = tof_seconds, .mu = mu_central,
                .short_way = true, .max_iter = 50, .tol = 1e-6f
            };
            lg_lambert_solution_t sol = lg_lambert_solve(&prob);
            
            if (sol.status == 0) {
                /* Compute C3 at departure */
                float v_esc1 = sqrtf(2.0f * mu_central / lg_vec3_len(r1));
                pt->c3 = lg_characteristic_energy(&sol.v1, v_esc1);
                pt->v_inf = lg_vec3_len(sol.v2) - sqrtf(2.0f * mu_central / lg_vec3_len(r2));
                pt->flag = 0;
                
                /* Update bounds */
                plot->c3_min = fminf(plot->c3_min, pt->c3);
                plot->c3_max = fmaxf(plot->c3_max, pt->c3);
                plot->v_inf_min = fminf(plot->v_inf_min, pt->v_inf);
                plot->v_inf_max = fmaxf(plot->v_inf_max, pt->v_inf);
                
                /* Track best solution (minimum C3) */
                if (pt->c3 < plot->best_c3) {
                    plot->best_c3 = pt->c3;
                    plot->best_launch_mjd = launch_mjd;
                    plot->best_arrival_mjd = arrival_mjd;
                }
            } else {
                pt->flag = 1; /* No solution */
                pt->c3 = -1.0f;
                pt->v_inf = -1.0f;
            }
        }
    }
}

static inline void lg_porkchop_free(lg_porkchop_plot_t* plot) {
    if (plot->grid) {
        free(plot->grid);
        plot->grid = NULL;
    }
}

/*============================================================================
 * 5. LOW-THRUST TRAJECTORIES - Equinoctial Elements
 *===========================================================================*/

/* Equinoctial elements: nonsingular for near-circular, near-equatorial */
typedef struct {
    float p;                  /* Semi-latus rectum (km) */
    float f;                  /* e*cos(w+Omega) */
    float g;                  /* e*sin(w+Omega) */
    float h;                  /* tan(i/2)*cos(Omega) */
    float k;                  /* tan(i/2)*sin(Omega) */
    float L;                  /* Mean longitude = M + w + Omega */
    float mass;               /* Current spacecraft mass (kg) */
} lg_equinoctial_t;

/* Convert classical to equinoctial elements */
static inline lg_equinoctial_t lg_classical_to_equinoctial(const lg_orbit_t* orb) {
    lg_equinoctial_t eq = {0};
    
    float ecc = orb->e;
    float inc = orb->i;
    float raan = orb->Omega;
    float argp = orb->w;
    
    float sin_raan = sinf(raan);
    float cos_raan = cosf(raan);
    float sin_argp_w = sinf(argp + raan);
    float cos_argp_w = cosf(argp + raan);
    float tan_half_i = tanf(inc * 0.5f);
    
    eq.p = orb->a * (1.0f - ecc * ecc);
    eq.f = ecc * cos_argp_w;
    eq.g = ecc * sin_argp_w;
    eq.h = tan_half_i * cos_raan;
    eq.k = tan_half_i * sin_raan;
    eq.L = orb->M + argp + raan;
    
    return eq;
}

/* Convert equinoctial back to position/velocity */
static inline void lg_equinoctial_to_rv(
    const lg_equinoctial_t* eq,
    float mu,
    lg_vec3_t* r,
    lg_vec3_t* v
) {
    float sin_L = sinf(eq->L);
    float cos_L = cosf(eq->L);
    float alpha2 = eq->h * eq->h - eq->k * eq->k;
    float s2 = 1.0f + eq->h * eq->h + eq->k * eq->k;
    float s = sqrtf(s2);
    
    /* Position in orbital plane */
    float X = eq->p / s * ((1.0f + alpha2) * cos_L + 2.0f * eq->h * eq->k * sin_L);
    float Y = eq->p / s * ((1.0f - alpha2) * sin_L + 2.0f * eq->h * eq->k * cos_L);
    
    /* Velocity in orbital plane */
    float mu_over_p = sqrtf(mu / eq->p);
    float Xdot = -mu_over_p / s * ((1.0f + alpha2) * sin_L - 2.0f * eq->h * eq->k * cos_L + eq->g - eq->f * eq->h * eq->k);
    float Ydot = -mu_over_p / s * ((1.0f - alpha2) * cos_L + 2.0f * eq->h * eq->k * sin_L - eq->f + eq->g * eq->h * eq->k);
    
    /* Rotation to inertial frame */
    float rot[3][3] = {
        {1.0f - 2.0f * eq->k * eq->k / s2, 2.0f * eq->h * eq->k / s2, 0.0f},
        {2.0f * eq->h * eq->k / s2, 1.0f - 2.0f * eq->h * eq->h / s2, 0.0f},
        {0.0f, 0.0f, 1.0f}
    };
    
    r->x = rot[0][0] * X + rot[0][1] * Y;
    r->y = rot[1][0] * X + rot[1][1] * Y;
    r->z = rot[2][0] * X + rot[2][1] * Y;
    
    v->x = rot[0][0] * Xdot + rot[0][1] * Ydot;
    v->y = rot[1][0] * Xdot + rot[1][1] * Ydot;
    v->z = rot[2][0] * Xdot + rot[2][1] * Ydot;
}

/* Gauss variational equations for equinoctial elements under thrust */
typedef struct {
    float d_p; float d_f; float d_g; float d_h; float d_k; float d_L; float d_mass;
} lg_equinoctial_derivs_t;

static inline lg_equinoctial_derivs_t lg_equinoctial_gauss(
    const lg_equinoctial_t* eq,
    float mu,
    const lg_vec3_t* accel_rtn  /* Acceleration in radial-tangential-normal frame */
) {
    lg_equinoctial_derivs_t d = {0};
    
    float r = eq->p / (1.0f + eq->f * cosf(eq->L) + eq->g * sinf(eq->L));
    float v_tangential = sqrtf(mu / eq->p) * (1.0f + eq->f * cosf(eq->L) + eq->g * sinf(eq->L));
    float n = sqrtf(mu / (eq->p * eq->p * eq->p)) * 
              powf(1.0f + eq->f * cosf(eq->L) + eq->g * sinf(eq->L), 2.0f);
    
    float s2 = 1.0f + eq->h * eq->h + eq->k * eq->k;
    float sin_L = sinf(eq->L);
    float cos_L = cosf(eq->L);
    
    /* Transform acceleration from RTN to equinoctial variations */
    float a_r = accel_rtn->x;  /* Radial */
    float a_t = accel_rtn->y;  /* Tangential */
    float a_n = accel_rtn->z;  /* Normal */
    
    d.d_p = 2.0f * sqrtf(eq->p / mu) * eq->p / (1.0f + eq->f * cos_L + eq->g * sin_L) * a_t;
    
    d.d_f = sqrtf(eq->p / mu) * (
        sin_L * a_r + 
        ((1.0f + eq->f * cos_L + eq->g * sin_L) * cos_L + eq->f) / (1.0f + eq->f * cos_L + eq->g * sin_L) * a_t -
        (eq->g * (eq->h * sin_L - eq->k * cos_L)) / s2 * a_n
    );
    
    d.d_g = sqrtf(eq->p / mu) * (
        -cos_L * a_r + 
        ((1.0f + eq->f * cos_L + eq->g * sin_L) * sin_L + eq->g) / (1.0f + eq->f * cos_L + eq->g * sin_L) * a_t +
        (eq->f * (eq->h * sin_L - eq->k * cos_L)) / s2 * a_n
    );
    
    d.d_h = sqrtf(eq->p / mu) * s2 * cos_L / (2.0f * (1.0f + eq->f * cos_L + eq->g * sin_L)) * a_n;
    d.d_k = sqrtf(eq->p / mu) * s2 * sin_L / (2.0f * (1.0f + eq->f * cos_L + eq->g * sin_L)) * a_n;
    
    d.d_L = n + sqrtf(eq->p / mu) * (
        ((eq->h * sin_L - eq->k * cos_L) * (eq->f * sin_L - eq->g * cos_L)) / 
        ((1.0f + eq->f * cos_L + eq->g * sin_L) * s2) * a_n
    );
    
    return d;
}

/* Solar sail acceleration: a = (2*P/c) * (A/m) * cos^2(alpha) * r_hat */
static inline lg_vec3_t lg_sail_acceleration(
    const lg_vec3_t* r,       /* Position from sun */
    float sail_area_m2,
    float sc_mass_kg,
    float reflectivity,       /* 0-1, 1=perfect reflection */
    float cone_angle          /* Sail normal to sun angle */
) {
    float r_mag = lg_vec3_len(*r);
    float solar_const = 1361.0f; /* W/m^2 at 1 AU */
    float c = 299792458.0f;      /* Speed of light */
    float AU = 149597870700.0f;  /* meters */
    
    float dist_au = r_mag / AU;
    float flux = solar_const / (dist_au * dist_au);
    float pressure = (1.0f + reflectivity) * flux / c; /* Radiation pressure */
    
    float accel_mag = pressure * sail_area_m2 / sc_mass_kg;
    float cos_cone = cosf(cone_angle);
    accel_mag *= cos_cone * cos_cone; /* cos^2 law for perfect sail */
    
    lg_vec3_t r_hat = lg_vec3_norm(*r);
    return lg_vec3_scale(r_hat, -accel_mag); /* Away from sun */
}

/*============================================================================
 * 6. AEROBRAKING & ATMOSPHERIC ENTRY
 *===========================================================================*/

typedef struct {
    float Cd;                 /* Drag coefficient */
    float Cl;                 /* Lift coefficient */
    float area_m2;            /* Reference area */
    float mass_kg;            /* Vehicle mass */
    float ballistic_coef;     /* m/(Cd*A) */
    float lift_drag_ratio;    /* Cl/Cd */
} lg_aero_vehicle_t;

/* lg_atmosphere_t and lg_atmosphere_density are defined in gravity.h */

/* Aerodynamic acceleration in planet-relative frame */
static inline lg_vec3_t lg_aero_acceleration(
    const lg_aero_vehicle_t* vehicle,
    const lg_atmosphere_t* atm,
    const lg_vec3_t* r_eci,       /* Position in ECI */
    const lg_vec3_t* v_eci,       /* Velocity in ECI */
    float planet_radius
) {
    float altitude = lg_vec3_len(*r_eci) - planet_radius;
    if (altitude < 0.0f) altitude = 0.0f;
    
    float rho = lg_atmosphere_density(atm, altitude);
    float v_rel_mag = lg_vec3_len(*v_eci); /* Simplified - should subtract rotation */
    
    float dynamic_pressure = 0.5f * rho * v_rel_mag * v_rel_mag;
    float drag_accel = dynamic_pressure * vehicle->Cd * vehicle->area_m2 / vehicle->mass_kg;
    
    /* Drag opposes velocity */
    lg_vec3_t v_hat = lg_vec3_norm(*v_eci);
    lg_vec3_t drag = lg_vec3_scale(v_hat, -drag_accel);
    
    /* Lift perpendicular to velocity and vertical (simplified) */
    lg_vec3_t up = lg_vec3_norm(*r_eci);
    lg_vec3_t lift_dir = lg_vec3_norm(lg_vec3_cross(lg_vec3_cross(v_hat, up), v_hat));
    lg_vec3_t lift = lg_vec3_scale(lift_dir, drag_accel * vehicle->lift_drag_ratio);
    
    return lg_vec3_add(drag, lift);
}

/* Entry corridor: check if periapsis is within aerocapture window */
typedef struct {
    float periapsis_min;      /* Minimum for skip-out */
    float periapsis_max;      /* Maximum for entry/descent */
    float heat_rate_limit;    /* W/cm^2 */
    float decel_limit;        /* g's */
    bool feasible;
} lg_entry_corridor_t;

static inline lg_entry_corridor_t lg_entry_corridor(
    const lg_orbit_t* approach,
    const lg_aero_vehicle_t* vehicle,
    const lg_atmosphere_t* atm,
    float v_entry           /* Entry velocity relative to planet */
) {
    lg_entry_corridor_t ec = {0};
    
    /* Periapsis range for aerocapture */
    float beta = vehicle->ballistic_coef;
    float H = (float)atm->H;
    
    /* Minimum: too high = skip out */
    ec.periapsis_min = H * logf(beta * v_entry * v_entry / (2.0f * H * 1.0f)); /* ~1 m/s^2 decel */
    
    /* Maximum: too low = entry/descent */
    float q_max = 500.0f; /* W/cm^2 typical heat limit */
    float v_circ = sqrtf(approach->mu / ((float)atm->body_radius + 100000.0f));
    ec.periapsis_max = H * logf(vehicle->Cd * atm->rho0 * v_entry * v_entry * v_entry / (2.0f * q_max));
    
    ec.feasible = ec.periapsis_min < ec.periapsis_max;
    
    return ec;
}

/*============================================================================
 * 7. TETHER PHYSICS - Skyhooks & Momentum Exchange
 *===========================================================================*/

typedef struct {
    float length;             /* Tether length (m) */
    float mass;               /* Tether mass (kg) */
    float tip_mass;           /* Mass at tether tip (kg) */
    float tensile_strength;   /* Material tensile strength (Pa) */
    float density;            /* Material density (kg/m^3) */
    float area;               /* Cross-sectional area (m^2) */
    float omega;              /* Rotation rate (rad/s) */
} lg_tether_t;

/* Tapered tether: area varies with tension for constant stress */
static inline float lg_tether_taper_area(
    float r,                  /* Distance from center of rotation */
    float r_tip,              /* Tip distance */
    float area_tip,           /* Area at tip */
    float tensile_strength,
    float density,
    float omega
) {
    /* Tension varies along tether due to centrifugal and gravity gradient */
    float a_centrifugal = omega * omega * r;
    float tension_integral = 0.5f * density * a_centrifugal * (r_tip * r_tip - r * r);
    return area_tip * expf(tension_integral / tensile_strength);
}

/* Skyhook tip velocity: V_tip = V_orbit + omega * L */
static inline float lg_skyhook_tip_velocity(
    float orbit_radius,
    float tether_length,
    float tether_omega,
    float mu_planet
) {
    float v_orbit = sqrtf(mu_planet / orbit_radius);
    float v_rotation = tether_omega * tether_length;
    return v_orbit + v_rotation;
}

/* Payload capture/release delta-v */
typedef struct {
    float capture_dv;         /* Delta-v at capture */
    float release_dv;         /* Delta-v at release */
    float momentum_exchange;  /* Change in tether angular momentum */
    float tip_velocity_pre;   /* Tip velocity before capture */
    float tip_velocity_post;  /* Tip velocity after capture (conservation) */
} lg_tether_exchange_t;

static inline lg_tether_exchange_t lg_tether_payload_exchange(
    const lg_tether_t* tether,
    float payload_mass,
    float v_payload,          /* Payload velocity relative to tether CM */
    bool capture              /* true = capture, false = release */
) {
    lg_tether_exchange_t ex = {0};
    
    float I_tether = tether->mass * tether->length * tether->length / 3.0f + 
                     tether->tip_mass * tether->length * tether->length;
    float L_initial = I_tether * tether->omega;
    
    float m_eff = capture ? tether->tip_mass + payload_mass : tether->tip_mass;
    float I_new = (tether->mass * tether->length * tether->length / 3.0f) + 
                  m_eff * tether->length * tether->length;
    
    /* Conservation of angular momentum */
    float omega_new = L_initial / I_new;
    ex.tip_velocity_pre = tether->omega * tether->length;
    ex.tip_velocity_post = omega_new * tether->length;
    
    ex.momentum_exchange = I_new * omega_new - L_initial;
    
    if (capture) {
        ex.capture_dv = fabsf(v_payload - ex.tip_velocity_pre);
        ex.release_dv = 0.0f;
    } else {
        ex.capture_dv = 0.0f;
        ex.release_dv = fabsf(v_payload - ex.tip_velocity_post);
    }
    
    return ex;
}

/* Rotovator: tether in elliptical orbit, tip touches atmosphere at periapsis */
typedef struct {
    float a;                  /* Semi-major axis */
    float e;                  /* Eccentricity */
    float tether_length;
    float omega_tether;       /* Tether spin rate (usually 2*orbital_mean_motion) */
    float v_tip_surface;      /* Tip velocity relative to surface at periapsis */
    float surface_contact_time; /* Time tip is within atmosphere */
} lg_rotovator_t;

static inline lg_rotovator_t lg_rotovator_design(
    float planet_radius,
    float atmosphere_height,
    float mu,
    float desired_tip_v     /* Desired payload velocity change */
) {
    lg_rotovator_t r = {0};
    
    /* Tip touches atmosphere at periapsis */
    float r_peri = planet_radius + atmosphere_height;
    float v_orbit_peri = sqrtf(2.0f * mu / r_peri); /* Escape velocity at periapsis */
    
    /* Tether length from desired delta-v */
    r.tether_length = desired_tip_v / (2.0f * M_PI); /* Rough: v = omega * L, omega ~ 2*pi/T */
    
    /* Orbit: apoapsis high, periapsis at atmosphere */
    r.e = 0.5f; /* Typical */
    r.a = r_peri / (1.0f - r.e);
    
    r.omega_tether = sqrtf(mu / (r.a * r.a * r.a)); /* Match orbital period */
    r.v_tip_surface = r.omega_tether * r.tether_length;
    
    /* Contact time: fraction of orbit near periapsis */
    float true_anomaly_range = 0.3f; /* ~17 degrees */
    r.surface_contact_time = true_anomaly_range / r.omega_tether;
    
    return r;
}

/*============================================================================
 * 8. INTEGRATION WITH PARTICLE SYSTEM (SoA thrusting bodies)
 *===========================================================================*/

/* Update particle with thrust acceleration */
static inline void lg_particle_apply_thrust(
    lg_particle_soa_t* soa,
    int idx,
    const lg_vec3_t* thrust_acceleration,  /* m/s^2 */
    float dt
) {
    /* Add thrust to velocity */
    lg_vec3_t dv = lg_vec3_scale(*thrust_acceleration, dt);
    soa->vx[idx] += dv.x;
    soa->vy[idx] += dv.y;
    soa->vz[idx] += dv.z;
    
    /* Update mass if propellant is consumed */
    /* (mass flow handled by caller) */
}

/* Low-thrust spiral: integrate equinoctial elements with small thrust */
static inline void lg_low_thrust_spiral(
    lg_equinoctial_t* eq,
    float thrust_accel,       /* Tangential acceleration magnitude */
    float Isp,
    float mu,
    float dt,
    int n_steps
) {
    float dt_sub = dt / n_steps;
    
    for (int i = 0; i < n_steps; i++) {
        /* Current radius from semi-latus rectum and true longitude */
        float r = eq->p / (1.0f + eq->f * cosf(eq->L) + eq->g * sinf(eq->L));
        
        /* Tangential thrust in RTN frame */
        lg_vec3_t accel_rtn = {0.0f, thrust_accel, 0.0f};
        
        /* Gauss variational equations */
        lg_equinoctial_derivs_t d = lg_equinoctial_gauss(eq, mu, &accel_rtn);
        
        /* Simple Euler integration (replace with RK4 for production) */
        eq->p += d.d_p * dt_sub;
        eq->f += d.d_f * dt_sub;
        eq->g += d.d_g * dt_sub;
        eq->h += d.d_h * dt_sub;
        eq->k += d.d_k * dt_sub;
        eq->L += d.d_L * dt_sub;
        
        /* Mass flow */
        float g0 = 9.80665f;
        float mdot = thrust_accel * eq->mass / (Isp * g0);
        eq->mass -= mdot * dt_sub;
    }
}

/*============================================================================
 * 9. UTILITY FUNCTIONS
 *===========================================================================*/

/* Delta-v budget for mission with multiple maneuvers */
static inline float lg_mission_delta_v(int n_maneuvers, const float* dv_each) {
    float total = 0.0f;
    for (int i = 0; i < n_maneuvers; i++) {
        total += dv_each[i];
    }
    return total;
}

/* Tsiolkovsky rocket equation: delta-v = Isp * g0 * ln(m0/mf) */
static inline float lg_tsiolkovsky_dv(float Isp, float m0, float mf) {
    const float g0 = 9.80665f;
    return Isp * g0 * logf(m0 / mf);
}

static inline float lg_tsiolkovsky_mass_ratio(float Isp, float dv) {
    const float g0 = 9.80665f;
    return expf(dv / (Isp * g0));
}

/* Bi-elliptic transfer: sometimes lower delta-v than Hohmann for large ratio */
static inline float lg_bielliptic_dv(
    float r1, float r2, float rb,  /* rb = intermediate apoapsis */
    float mu
) {
    float v1 = sqrtf(mu / r1);
    float v2 = sqrtf(mu / r2);
    
    /* Burn 1: circular to intermediate ellipse */
    float a1 = (r1 + rb) * 0.5f;
    float v_peri1 = sqrtf(mu * (2.0f/r1 - 1.0f/a1));
    float dv1 = v_peri1 - v1;
    
    /* Burn 2: at rb, drop to r2 */
    float a2 = (rb + r2) * 0.5f;
    float v_apo1 = sqrtf(mu * (2.0f/rb - 1.0f/a1));
    float v_apo2 = sqrtf(mu * (2.0f/rb - 1.0f/a2));
    float dv2 = fabsf(v_apo2 - v_apo1);
    
    /* Burn 3: circularize at r2 */
    float v_peri2 = sqrtf(mu * (2.0f/r2 - 1.0f/a2));
    float dv3 = v2 - v_peri2;
    
    return dv1 + dv2 + dv3;
}

/* Patched conic with thrust: propagate through sphere of influence with low thrust */
typedef struct {
    lg_orbit_t orbit_planet;      /* Planet-centered during SOI */
    lg_vec3_t v_inf_departure;    /* Excess velocity leaving planet */
    lg_vec3_t v_inf_arrival;      /* Excess velocity approaching next planet */
    float dv_total;               /* Total delta-v during transfer */
    float time_of_flight;         /* Total transfer time */
} lg_patched_conic_thrust_t;

#ifdef __cplusplus
}
#endif

#endif /* LAGRANGE_PROPULSION_H */

