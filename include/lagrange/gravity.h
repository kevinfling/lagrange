/*
 * Lagrange Physics Library - Gravity and Orbital Mechanics
 * Newtonian gravity and perturbations
 */

#ifndef LAGRANGE_GRAVITY_H
#define LAGRANGE_GRAVITY_H

#include "types.h"
#include "math.h"

#ifdef __cplusplus
extern "C" {
#endif

/*============================================================================
 * Physical Constants
 *===========================================================================*/

#define LG_G 6.67430e-11        /* Gravitational constant (m³/(kg·s²)) */
#define LG_G_AU 1.48818e-34     /* Gravitational constant in AU³/(Msun·day²) */

/* Common gravitational parameters (mu = G * M) */
#define LG_MU_SUN 1.32712440018e20      /* m³/s² */
#define LG_MU_EARTH 3.986004418e14      /* m³/s² */
#define LG_MU_MOON 4.9048695e12         /* m³/s² */
#define LG_MU_MARS 4.282837e13          /* m³/s² */
#define LG_MU_JUPITER 1.26686534e17     /* m³/s² */

/* Body radii (meters) */
#define LG_RADIUS_SUN 6.957e8
#define LG_RADIUS_EARTH 6.371e6
#define LG_RADIUS_MOON 1.737e6
#define LG_RADIUS_MARS 3.3895e6
#define LG_RADIUS_JUPITER 6.9911e7

/* Solar radiation pressure at 1 AU (N/m²) */
#define LG_SRP_1AU 4.560e-6

/* Speed of light (m/s) */
#define LG_SPEED_OF_LIGHT 2.99792458e8

/* Astronomical unit (m) */
#define LG_AU 1.495978707e11

/*============================================================================
 * Point Gravity (Newton's Law)
 * F = G * m1 * m2 / r²  directed along r
 *===========================================================================*/

/* Compute gravitational force on 'body' from 'source' at position 'body_pos' */
static inline lg_vec3_t lg_gravity_point(
    lg_vec3_t body_pos,
    float body_mass,
    lg_vec3_t source_pos,
    double source_mu
) {
    lg_vec3d_t r = lg_vec3d_sub(lg_vec3_to_d(source_pos), lg_vec3_to_d(body_pos));
    double dist_sq = lg_vec3d_len_sq(r);
    double dist = sqrt(dist_sq);
    
    if (dist < 1.0) {
        return lg_vec3_zero(); /* Too close, prevent singularity */
    }
    
    /* F = mu * m / r² in direction of r */
    double force_mag = source_mu * body_mass / dist_sq;
    lg_vec3d_t force = lg_vec3d_scale(lg_vec3d_norm(r), force_mag);
    
    return lg_vec3_from_d(force);
}

/* Compute acceleration (force per unit mass) - more common for orbital mechanics */
static inline lg_vec3_t lg_gravity_accel(
    lg_vec3_t body_pos,
    lg_vec3_t source_pos,
    double source_mu
) {
    lg_vec3d_t r = lg_vec3d_sub(lg_vec3_to_d(source_pos), lg_vec3_to_d(body_pos));
    double dist_sq = lg_vec3d_len_sq(r);
    double dist = sqrt(dist_sq);
    
    if (dist < 1.0) {
        return lg_vec3_zero();
    }
    
    /* a = mu / r² in direction of r */
    double accel_mag = source_mu / dist_sq;
    lg_vec3d_t accel = lg_vec3d_scale(lg_vec3d_norm(r), accel_mag);
    
    return lg_vec3_from_d(accel);
}

/* Compute orbital velocity for circular orbit */
static inline double lg_orbital_velocity_circular(double mu, double radius) {
    return sqrt(mu / radius);
}

/* Compute escape velocity */
static inline double lg_escape_velocity(double mu, double radius) {
    return sqrt(2.0 * mu / radius);
}

/* Compute orbital period (Kepler's 3rd law) */
static inline double lg_orbital_period(double mu, double semi_major_axis) {
    return 2.0 * LG_PI * sqrt(semi_major_axis * semi_major_axis * semi_major_axis / mu);
}

/* Mean motion n = sqrt(mu / a³) */
static inline double lg_mean_motion(double mu, double semi_major_axis) {
    return sqrt(mu / (semi_major_axis * semi_major_axis * semi_major_axis));
}

/* Mean anomaly M = n * dt, wrapped to [0, 2π) */
static inline double lg_mean_anomaly(double mean_motion, double dt) {
    double M = fmod(mean_motion * dt, 2.0 * LG_PI);
    if (M < 0.0) M += 2.0 * LG_PI;
    return M;
}

/* Solve Kepler's equation M = E - e*sin(E) using Newton-Raphson */
static inline double lg_kepler_solve(double M, double eccentricity, double tol, int max_iter) {
    double E = M;
    if (eccentricity > 0.8) {
        E = LG_PI;
    }
    for (int i = 0; i < max_iter; i++) {
        double f = E - eccentricity * sin(E) - M;
        double fp = 1.0 - eccentricity * cos(E);
        double dE = -f / fp;
        E += dE;
        if (fabs(dE) < tol) break;
    }
    return E;
}

/* True anomaly from eccentric anomaly */
static inline double lg_true_anomaly_from_eccentric(double E, double e) {
    double cos_E = cos(E);
    double sin_E = sin(E);
    double cos_nu = (cos_E - e) / (1.0 - e * cos_E);
    double sin_nu = (sqrt(1.0 - e * e) * sin_E) / (1.0 - e * cos_E);
    return atan2(sin_nu, cos_nu);
}

/* Flight path angle: angle between velocity vector and local horizontal */
static inline double lg_flight_path_angle(double e, double true_anomaly) {
    return atan2(e * sin(true_anomaly), 1.0 + e * cos(true_anomaly));
}

/*============================================================================
 * J2/J3/J4 Perturbation (Zonal Harmonics)
 * Correction for non-spherical mass distribution
 *===========================================================================*/

typedef struct {
    double mu;              /* Gravitational parameter */
    double radius;          /* Equatorial radius */
    double j2;              /* J2 coefficient (Earth: ~1.08263e-3) */
    double j3;              /* J3 coefficient (Earth: ~-2.5327e-6) */
    double j4;              /* J4 coefficient (Earth: ~-1.6196e-6) */
} lg_gravity_body_t;

/* Earth parameters */
static const lg_gravity_body_t LG_GRAVITY_EARTH = {
    .mu = LG_MU_EARTH,
    .radius = LG_RADIUS_EARTH,
    .j2 = 1.08263e-3,
    .j3 = -2.5327e-6,
    .j4 = -1.6196e-6
};

/* Compute J2 perturbation acceleration */
static inline lg_vec3_t lg_gravity_j2(
    lg_vec3_t body_pos,
    const lg_gravity_body_t* body
) {
    lg_vec3d_t r = lg_vec3_to_d(body_pos);
    double dist = lg_vec3d_len(r);
    
    if (dist < body->radius) {
        return lg_vec3_zero(); /* Inside body, no J2 effect */
    }
    
    double r_ratio = body->radius / dist;
    double r_ratio_sq = r_ratio * r_ratio;
    double factor = 1.5 * body->j2 * r_ratio_sq * (body->mu / (dist * dist));
    
    /* z²/r² term */
    double z_over_r = r.z / dist;
    double z_term = 5.0 * z_over_r * z_over_r - 1.0;
    
    lg_vec3d_t accel;
    accel.x = factor * (r.x / dist) * z_term;
    accel.y = factor * (r.y / dist) * z_term;
    accel.z = factor * (r.z / dist) * (z_term - 2.0);
    
    return lg_vec3_from_d(accel);
}

/* Compute J3 perturbation acceleration */
static inline lg_vec3_t lg_gravity_j3(
    lg_vec3_t body_pos,
    const lg_gravity_body_t* body
) {
    lg_vec3d_t r = lg_vec3_to_d(body_pos);
    double dist = lg_vec3d_len(r);
    
    if (dist < body->radius) {
        return lg_vec3_zero();
    }
    
    double s = r.z / dist;
    double factor = -0.5 * body->j3 * body->mu * pow(body->radius / dist, 3) / (dist * dist);
    
    lg_vec3d_t accel;
    accel.x = factor * 5.0 * (r.x / dist) * (r.z / dist) * (3.0 - 7.0 * s * s);
    accel.y = factor * 5.0 * (r.y / dist) * (r.z / dist) * (3.0 - 7.0 * s * s);
    accel.z = factor * 3.0 * (1.0 - 10.0 * s * s + (35.0 / 3.0) * s * s * s * s);
    
    return lg_vec3_from_d(accel);
}

/* Compute J4 perturbation acceleration */
static inline lg_vec3_t lg_gravity_j4(
    lg_vec3_t body_pos,
    const lg_gravity_body_t* body
) {
    lg_vec3d_t r = lg_vec3_to_d(body_pos);
    double dist = lg_vec3d_len(r);
    
    if (dist < body->radius) {
        return lg_vec3_zero();
    }
    
    double s = r.z / dist;
    double factor = -0.5 * body->j4 * body->mu * pow(body->radius / dist, 4) / (dist * dist);
    double poly_xy = 3.0 - 42.0 * s * s + 63.0 * s * s * s * s;
    double poly_z  = 15.0 - 70.0 * s * s + 63.0 * s * s * s * s;
    
    lg_vec3d_t accel;
    accel.x = factor * (5.0 / 4.0) * (r.x / dist) * poly_xy;
    accel.y = factor * (5.0 / 4.0) * (r.y / dist) * poly_xy;
    accel.z = factor * (5.0 / 4.0) * (r.z / dist) * poly_z;
    
    return lg_vec3_from_d(accel);
}

/* Combined gravity + J2 */
static inline lg_vec3_t lg_gravity_with_j2(
    lg_vec3_t body_pos,
    const lg_gravity_body_t* body
) {
    lg_vec3_t point = lg_gravity_accel(body_pos, lg_vec3_zero(), body->mu);
    lg_vec3_t j2 = lg_gravity_j2(body_pos, body);
    return lg_vec3_add(point, j2);
}

/* Combined gravity + J2 + J3 + J4 */
static inline lg_vec3_t lg_gravity_with_zonal(
    lg_vec3_t body_pos,
    const lg_gravity_body_t* body
) {
    lg_vec3_t point = lg_gravity_accel(body_pos, lg_vec3_zero(), body->mu);
    lg_vec3_t j2 = lg_gravity_j2(body_pos, body);
    lg_vec3_t j3 = lg_gravity_j3(body_pos, body);
    lg_vec3_t j4 = lg_gravity_j4(body_pos, body);
    lg_vec3_t sum = lg_vec3_add(point, j2);
    sum = lg_vec3_add(sum, j3);
    sum = lg_vec3_add(sum, j4);
    return sum;
}

/*============================================================================
 * Atmospheric Drag
 *===========================================================================*/

typedef struct {
    double rho0;        /* Reference density (kg/m³) */
    double h0;          /* Reference altitude (m) */
    double H;           /* Scale height (m) */
    double body_radius; /* Planet radius (m) */
} lg_atmosphere_t;

/* Standard exponential atmosphere model */
static inline double lg_atmosphere_density(const lg_atmosphere_t* atm, double altitude) {
    return atm->rho0 * exp(-(altitude - atm->h0) / atm->H);
}

/* Compute drag acceleration: a = -0.5 * rho * v² * Cd * A / m * v̂ */
static inline lg_vec3_t lg_gravity_drag(
    lg_vec3_t body_pos,
    lg_vec3_t velocity,
    double cd,
    double area,
    double mass,
    const lg_atmosphere_t* atm
) {
    double r = lg_vec3d_len(lg_vec3_to_d(body_pos));
    double alt = r - atm->body_radius;
    if (alt < 0.0) {
        return lg_vec3_zero(); /* Underground */
    }
    
    double rho = lg_atmosphere_density(atm, alt);
    lg_vec3d_t v = lg_vec3_to_d(velocity);
    double v_mag = lg_vec3d_len(v);
    if (v_mag < 1e-10) {
        return lg_vec3_zero();
    }
    
    double factor = -0.5 * cd * area * rho * v_mag / mass;
    lg_vec3d_t accel = lg_vec3d_scale(lg_vec3d_norm(v), factor * v_mag);
    return lg_vec3_from_d(accel);
}

/*============================================================================
 * Solar Radiation Pressure (Cannonball Model)
 *===========================================================================*/

/* Compute SRP acceleration given vector from Sun to satellite */
static inline lg_vec3_t lg_gravity_srp(
    lg_vec3_t sun_to_sat,     /* Vector from Sun to satellite (m) */
    double cr,                /* Reflectivity coefficient (1.0 = perfect absorption, 2.0 = perfect reflection) */
    double area,              /* Cross-sectional area (m²) */
    double mass               /* Satellite mass (kg) */
) {
    lg_vec3d_t r = lg_vec3_to_d(sun_to_sat);
    double d = lg_vec3d_len(r);
    if (d < 1.0) {
        return lg_vec3_zero();
    }
    
    double pressure = LG_SRP_1AU * (LG_AU * LG_AU) / (d * d);
    double factor = cr * area * pressure / mass;
    lg_vec3d_t accel = lg_vec3d_scale(lg_vec3d_norm(r), factor);
    return lg_vec3_from_d(accel);
}

/*============================================================================
 * Third-Body Perturbation
 * a = mu_third * (d_vec/|d|³ - r_third/|r_third|³)
 *===========================================================================*/

static inline lg_vec3_t lg_gravity_third_body(
    lg_vec3_t sat_pos,        /* Satellite position relative to primary */
    lg_vec3_t third_pos,      /* Third body position relative to primary */
    double third_mu           /* Third body gravitational parameter */
) {
    lg_vec3d_t r_s = lg_vec3_to_d(sat_pos);
    lg_vec3d_t r_t = lg_vec3_to_d(third_pos);
    lg_vec3d_t d_vec = lg_vec3d_sub(r_t, r_s); /* vector from sat to third body */
    double d = lg_vec3d_len(d_vec);
    double r_t_mag = lg_vec3d_len(r_t);
    
    if (d < 1.0 || r_t_mag < 1.0) {
        return lg_vec3_zero();
    }
    
    /* Perturbation = mu_third * (d_vec/d^3 + r_t/r_t^3) */
    lg_vec3d_t term1 = lg_vec3d_scale(lg_vec3d_norm(d_vec), 1.0 / (d * d));
    lg_vec3d_t term2 = lg_vec3d_scale(lg_vec3d_norm(r_t), 1.0 / (r_t_mag * r_t_mag));
    
    lg_vec3d_t accel = lg_vec3d_scale(lg_vec3d_add(term1, term2), third_mu);
    return lg_vec3_from_d(accel);
}

/*============================================================================
 * Sphere of Influence
 * For patched conic approximation (e.g., spacecraft leaving Earth for Mars)
 *===========================================================================*/

/* Compute radius of sphere of influence */
static inline double lg_sphere_of_influence(
    double body_orbit_radius,   /* Distance from body to its parent */
    double body_mass,
    double parent_mass
) {
    return body_orbit_radius * pow(body_mass / parent_mass, 0.4);
}

/* Hill sphere: stable satellite orbits */
static inline double lg_hill_sphere(
    double semi_major_axis,
    double eccentricity,
    double body_mass,
    double primary_mass
) {
    return semi_major_axis * (1.0 - eccentricity) * pow(body_mass / (3.0 * primary_mass), 1.0 / 3.0);
}

/* Earth's SOI radius (~0.93 million km) */
#define LG_EARTH_SOI_RADIUS 9.29e8

/* Synodic period of two orbiting bodies */
static inline double lg_synodic_period(double T1, double T2) {
    return 1.0 / fabs(1.0 / T1 - 1.0 / T2);
}

/*============================================================================
 * Orbital Elements
 * Convert between state vectors (position, velocity) and orbital elements
 *===========================================================================*/

typedef struct {
    double semi_major_axis;     /* a (meters) */
    double eccentricity;        /* e (unitless) */
    double inclination;         /* i (radians) */
    double raan;                /* Ω - right ascension of ascending node (radians) */
    double arg_of_periapsis;    /* ω (radians) */
    double true_anomaly;        /* ν (radians) */
} lg_orbital_elements_t;

/* Compact orbit representation for propagation (used by scene graph & propulsion) */
typedef struct {
    float a;        /* Semi-major axis (m) */
    float e;        /* Eccentricity */
    float i;        /* Inclination (rad) */
    float Omega;    /* RAAN (rad) */
    float w;        /* Argument of periapsis (rad) */
    float M;        /* Current mean anomaly (rad) */
    float M0;       /* Mean anomaly at epoch (rad) */
    float n;        /* Mean motion (rad/s) */
    float mu;       /* Gravitational parameter (m³/s²) */
    double epoch;   /* Reference epoch (seconds) */
} lg_orbit_t;

/* Compute position from mean anomaly using classical orbital elements.
 * Solves Kepler's equation and rotates to inertial frame. */
static inline lg_vec3_t lg_orbit_position_from_mean_anomaly(const lg_orbit_t* orb, float M) {
    /* Solve Kepler's equation for eccentric anomaly E */
    float E = M;
    for (int iter = 0; iter < 20; iter++) {
        float dE = (E - orb->e * sinf(E) - M) / (1.0f - orb->e * cosf(E));
        E -= dE;
        if (fabsf(dE) < 1e-8f) break;
    }
    
    /* True anomaly: nu = 2*atan2(sqrt(1+e)*sin(E/2), sqrt(1-e)*cos(E/2)) */
    float half_E = E * 0.5f;
    float nu = 2.0f * atan2f(sinf(half_E) * sqrtf(1.0f + orb->e), cosf(half_E) * sqrtf(1.0f - orb->e));
    
    float r = orb->a * (1.0f - orb->e * cosf(E));
    
    /* Position in orbital plane */
    float x_orb = r * cosf(nu);
    float y_orb = r * sinf(nu);
    
    /* Rotate to inertial frame (3-1-3 Euler angles: Omega, i, w) */
    float cos_O = cosf(orb->Omega);
    float sin_O = sinf(orb->Omega);
    float cos_i = cosf(orb->i);
    float sin_i = sinf(orb->i);
    float cos_w = cosf(orb->w);
    float sin_w = sinf(orb->w);
    
    float x = (cos_O * cos_w - sin_O * sin_w * cos_i) * x_orb +
              (-cos_O * sin_w - sin_O * cos_w * cos_i) * y_orb;
    float y = (sin_O * cos_w + cos_O * sin_w * cos_i) * x_orb +
              (-sin_O * sin_w + cos_O * cos_w * cos_i) * y_orb;
    float z = (sin_i * sin_w) * x_orb + (sin_i * cos_w) * y_orb;
    
    return lg_vec3(x, y, z);
}

/* Compute orbital elements from state vectors */
static inline lg_orbital_elements_t lg_state_to_elements(
    lg_vec3_t position,
    lg_vec3_t velocity,
    double mu
) {
    lg_orbital_elements_t elem = {0};
    
    lg_vec3d_t r = lg_vec3_to_d(position);
    lg_vec3d_t v = lg_vec3_to_d(velocity);
    
    double r_mag = lg_vec3d_len(r);
    double v_sq = lg_vec3d_len_sq(v);
    
    /* Specific angular momentum */
    lg_vec3d_t h = lg_vec3d_cross(r, v);
    double h_mag = lg_vec3d_len(h);
    
    /* Semi-major axis from vis-viva equation */
    double energy = v_sq / 2.0 - mu / r_mag;
    if (fabs(energy) > 1e-15) {
        elem.semi_major_axis = -mu / (2.0 * energy);
    }
    
    /* Eccentricity vector */
    lg_vec3d_t e_vec = lg_vec3d_sub(
        lg_vec3d_cross(v, h),
        lg_vec3d_scale(r, mu / r_mag)
    );
    e_vec = lg_vec3d_scale(e_vec, 1.0 / mu);
    elem.eccentricity = lg_vec3d_len(e_vec);
    
    /* Inclination */
    elem.inclination = acos(h.z / h_mag);
    
    /* RAAN (if inclination is not zero) */
    double n_mag = sqrt(h.x * h.x + h.y * h.y);
    if (n_mag > 1e-15) {
        elem.raan = atan2(h.x, -h.y);
        if (elem.raan < 0) elem.raan += 2.0 * LG_PI;
    }
    
    /* Argument of periapsis and true anomaly */
    if (elem.eccentricity > 1e-10 && n_mag > 1e-15) {
        double cos_omega = (h.x * e_vec.y - h.y * e_vec.x) / 
                           (n_mag * elem.eccentricity);
        cos_omega = fmax(-1.0, fmin(1.0, cos_omega));
        elem.arg_of_periapsis = acos(cos_omega);
        if (e_vec.z < 0.0) {
            elem.arg_of_periapsis = 2.0 * LG_PI - elem.arg_of_periapsis;
        }
        
        double cos_nu = lg_vec3d_dot(e_vec, r) / 
                        (elem.eccentricity * r_mag);
        cos_nu = fmax(-1.0, fmin(1.0, cos_nu));
        elem.true_anomaly = acos(cos_nu);
        if (lg_vec3d_dot(r, v) < 0.0) {
            elem.true_anomaly = 2.0 * LG_PI - elem.true_anomaly;
        }
    } else if (n_mag > 1e-15) {
        /* Near-circular orbit: use argument of latitude */
        elem.arg_of_periapsis = 0.0;
        double cos_u = (h.y * r.x - h.x * r.y) / (n_mag * r_mag);
        cos_u = fmax(-1.0, fmin(1.0, cos_u));
        elem.true_anomaly = acos(cos_u);
        if (lg_vec3d_dot(r, v) < 0.0) {
            elem.true_anomaly = 2.0 * LG_PI - elem.true_anomaly;
        }
    }
    
    return elem;
}

/* Compute state vectors from orbital elements */
static inline void lg_elements_to_state(
    const lg_orbital_elements_t* elem,
    double mu,
    lg_vec3_t* out_pos,
    lg_vec3_t* out_vel
) {
    /* Position and velocity in orbital plane */
    double e = elem->eccentricity;
    double a = elem->semi_major_axis;
    double nu = elem->true_anomaly;
    
    double cos_nu = cos(nu);
    double sin_nu = sin(nu);
    
    double p = a * (1.0 - e * e); /* Semi-latus rectum */
    double r = p / (1.0 + e * cos_nu);
    
    /* Position in orbital plane */
    lg_vec3d_t pos_plane = {r * cos_nu, r * sin_nu, 0.0};
    
    /* Velocity in orbital plane */
    double mu_over_h = sqrt(mu / p);
    lg_vec3d_t vel_plane = {
        mu_over_h * (-sin_nu),
        mu_over_h * (e + cos_nu),
        0.0
    };
    
    /* Rotation matrix from orbital to inertial frame */
    double cos_O = cos(elem->raan);
    double sin_O = sin(elem->raan);
    double cos_w = cos(elem->arg_of_periapsis);
    double sin_w = sin(elem->arg_of_periapsis);
    double cos_i = cos(elem->inclination);
    double sin_i = sin(elem->inclination);
    
    /* Apply rotations */
    lg_vec3d_t pos, vel;
    
    /* Rotation about z by arg of periapsis */
    double x1 = cos_w * pos_plane.x - sin_w * pos_plane.y;
    double y1 = sin_w * pos_plane.x + cos_w * pos_plane.y;
    double vx1 = cos_w * vel_plane.x - sin_w * vel_plane.y;
    double vy1 = sin_w * vel_plane.x + cos_w * vel_plane.y;
    
    /* Rotation about x by inclination */
    double y2 = cos_i * y1;
    double z2 = sin_i * y1;
    double vy2 = cos_i * vy1;
    double vz2 = sin_i * vy1;
    
    /* Rotation about z by RAAN */
    pos.x = cos_O * x1 - sin_O * y2;
    pos.y = sin_O * x1 + cos_O * y2;
    pos.z = z2;
    
    vel.x = cos_O * vx1 - sin_O * vy2;
    vel.y = sin_O * vx1 + cos_O * vy2;
    vel.z = vz2;
    
    *out_pos = lg_vec3_from_d(pos);
    *out_vel = lg_vec3_from_d(vel);
}

/*============================================================================
 * Gauss's Method for Orbit Determination
 * Angles-only observations from three positions
 *===========================================================================*/

typedef struct {
    lg_vec3_t L;      /* Line-of-sight unit vector */
    lg_vec3_t R;      /* Observer position */
    double t;         /* Observation time */
} lg_optical_obs_t;

/*
 * Gauss's method for orbit determination from three optical observations.
 * Returns true on success, false if observations are degenerate (coplanar).
 * Outputs the satellite state vector at the middle observation time.
 */
static inline double lg_gauss_poly_eval(double r, double a, double b, double c) {
    double r2 = r * r;
    double r3 = r2 * r;
    double r6 = r3 * r3;
    return r6 * r2 + a * r6 + b * r3 + c;
}

static inline int lg_gauss_find_root(double a, double b, double c, double* root) {
    double r = 1e-6;
    double f_prev = lg_gauss_poly_eval(r, a, b, c);
    for (int i = 0; i < 60; i++) {
        double r_next = r * 2.0;
        double f_next = lg_gauss_poly_eval(r_next, a, b, c);
        if (f_prev * f_next < 0.0) {
            double lo = r, hi = r_next;
            double f_lo = f_prev;
            for (int j = 0; j < 80; j++) {
                double mid = (lo + hi) * 0.5;
                double f_mid = lg_gauss_poly_eval(mid, a, b, c);
                if (f_lo * f_mid < 0.0) {
                    hi = mid;
                } else {
                    lo = mid;
                    f_lo = f_mid;
                }
            }
            *root = (lo + hi) * 0.5;
            return 1;
        }
        r = r_next;
        f_prev = f_next;
    }
    return 0;
}

static inline int lg_gauss_orbit_determination(
    const lg_optical_obs_t obs[3],
    double mu,
    lg_vec3_t* out_pos,
    lg_vec3_t* out_vel
) {
    lg_vec3d_t L1 = lg_vec3_to_d(obs[0].L);
    lg_vec3d_t L2 = lg_vec3_to_d(obs[1].L);
    lg_vec3d_t L3 = lg_vec3_to_d(obs[2].L);

    lg_vec3d_t R1 = lg_vec3_to_d(obs[0].R);
    lg_vec3d_t R2 = lg_vec3_to_d(obs[1].R);
    lg_vec3d_t R3 = lg_vec3_to_d(obs[2].R);

    double tau1 = obs[0].t - obs[1].t;
    double tau3 = obs[2].t - obs[1].t;
    double tau  = tau3 - tau1;

    lg_vec3d_t p0 = lg_vec3d_cross(L2, L3);
    lg_vec3d_t p1 = lg_vec3d_cross(L1, L3);
    lg_vec3d_t p2 = lg_vec3d_cross(L1, L2);

    double D0 = lg_vec3d_dot(L1, p0);
    if (fabs(D0) < 1e-15) {
        return 0; /* Degenerate: observations are coplanar */
    }

    double D11 = lg_vec3d_dot(R1, p0);
    double D12 = lg_vec3d_dot(R1, p1);
    double D13 = lg_vec3d_dot(R1, p2);

    double D21 = lg_vec3d_dot(R2, p0);
    double D22 = lg_vec3d_dot(R2, p1);
    double D23 = lg_vec3d_dot(R2, p2);

    double D31 = lg_vec3d_dot(R3, p0);
    double D32 = lg_vec3d_dot(R3, p1);
    double D33 = lg_vec3d_dot(R3, p2);

    double A = (-D12 * (tau3 / tau) + D22 + D32 * (tau1 / tau)) / D0;
    double B = (D12 * (tau3 * tau3 - tau * tau) * (tau3 / tau)
                + D32 * (tau * tau - tau1 * tau1) * (tau1 / tau)) / (6.0 * D0);

    double E = lg_vec3d_dot(R2, L2);
    double R2_sq = lg_vec3d_dot(R2, R2);

    double a_poly = -(A * A + 2.0 * A * E + R2_sq);
    double b_poly = -2.0 * mu * B * (A + E);
    double c_poly = -(mu * mu) * (B * B);

    double r2_star;
    if (!lg_gauss_find_root(a_poly, b_poly, c_poly, &r2_star)) {
        return 0;
    }

    double r2_cubed = r2_star * r2_star * r2_star;

    double num1 = 6.0 * (D31 * (tau1 / tau3) + D21 * (tau / tau3)) * r2_cubed
                + mu * D31 * (tau * tau - tau1 * tau1) * (tau1 / tau3);
    double den1 = 6.0 * r2_cubed + mu * (tau * tau - tau3 * tau3);
    double rho1 = ((num1 / den1) - D11) / D0;

    double rho2 = A + mu * B / r2_cubed;

    double num3 = 6.0 * (D13 * (tau3 / tau1) - D23 * (tau / tau1)) * r2_cubed
                + mu * D13 * (tau * tau - tau3 * tau3) * (tau3 / tau1);
    double den3 = 6.0 * r2_cubed + mu * (tau * tau - tau1 * tau1);
    double rho3 = ((num3 / den3) - D33) / D0;

    if (rho1 < 0.0 || rho2 < 0.0 || rho3 < 0.0 || !isfinite(rho1) || !isfinite(rho2) || !isfinite(rho3)) {
        return 0; /* Negative or non-finite ranges */
    }

    lg_vec3d_t r1_vec = lg_vec3d_add(R1, lg_vec3d_scale(L1, rho1));
    lg_vec3d_t r2_vec = lg_vec3d_add(R2, lg_vec3d_scale(L2, rho2));
    lg_vec3d_t r3_vec = lg_vec3d_add(R3, lg_vec3d_scale(L3, rho3));

    double f1 = 1.0 - 0.5 * (mu / r2_cubed) * tau1 * tau1;
    double f3 = 1.0 - 0.5 * (mu / r2_cubed) * tau3 * tau3;
    double g1 = tau1 - (1.0 / 6.0) * (mu / r2_cubed) * tau1 * tau1 * tau1;
    double g3 = tau3 - (1.0 / 6.0) * (mu / r2_cubed) * tau3 * tau3 * tau3;

    double denom = f1 * g3 - f3 * g1;
    if (fabs(denom) < 1e-15) return 0;

    lg_vec3d_t v2_vec = lg_vec3d_sub(lg_vec3d_scale(r3_vec, f1), lg_vec3d_scale(r1_vec, f3));
    v2_vec = lg_vec3d_scale(v2_vec, 1.0 / denom);

    *out_pos = lg_vec3_from_d(r2_vec);
    *out_vel = lg_vec3_from_d(v2_vec);
    return 1;
}

#ifdef __cplusplus
}
#endif

#endif /* LAGRANGE_GRAVITY_H */
