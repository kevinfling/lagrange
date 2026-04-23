/*
 * Lagrange Tests - Gravity and Orbital Mechanics
 */

#include <stdio.h>
#include <math.h>
#include <stdbool.h>

#define LAGRANGE_IMPLEMENTATION
#include "lagrange.h"

#define EPSILON 1e-3f
#define EPSILON_LOOSE 1e-2f
#define D_EPSILON 1e-6

#define ASSERT_NEAR(a, b, eps) do { \
    float _a = (a), _b = (b); \
    if (fabsf(_a - _b) > (eps)) { \
        printf("FAIL: %s:%d: Expected %f, got %f (diff: %e)\n", \
               __FILE__, __LINE__, _b, _a, fabsf(_a - _b)); \
        return 1; \
    } \
} while(0)

#define ASSERT_NEAR_D(a, b, eps) do { \
    double _a = (a), _b = (b); \
    if (fabs(_a - _b) > (eps)) { \
        printf("FAIL: %s:%d: Expected %lf, got %lf (diff: %e)\n", \
               __FILE__, __LINE__, _b, _a, fabs(_a - _b)); \
        return 1; \
    } \
} while(0)

#define ASSERT_VEC3_NEAR(v, vx, vy, vz, eps) do { \
    ASSERT_NEAR((v).x, (vx), (eps)); \
    ASSERT_NEAR((v).y, (vy), (eps)); \
    ASSERT_NEAR((v).z, (vz), (eps)); \
} while(0)

#define ASSERT_TRUE(cond) do { \
    if (!(cond)) { \
        printf("FAIL: %s:%d: Assertion failed: %s\n", __FILE__, __LINE__, #cond); \
        return 1; \
    } \
} while(0)

/*==========================================================================
 * J2 Perturbation
 *==========================================================================*/

int test_j2_perturbation(void) {
    /* Satellite at equator, 400 km altitude */
    double r = LG_RADIUS_EARTH + 400e3;
    lg_vec3_t pos = { (float)r, 0.0f, 0.0f };
    
    lg_vec3_t accel = lg_gravity_j2(pos, &LG_GRAVITY_EARTH);
    
    /* J2 at equator should have zero z-component and point inward (negative x) */
    ASSERT_NEAR(accel.x, -0.0125f, 0.005f); /* ~1.2 cm/s² inward at 400 km */
    ASSERT_NEAR(accel.y, 0.0f, 1e-6f);
    ASSERT_NEAR(accel.z, 0.0f, 1e-6f);
    
    /* At pole, J2 z-component should be twice the radial magnitude */
    lg_vec3_t pole = { 0.0f, 0.0f, (float)r };
    lg_vec3_t accel_pole = lg_gravity_j2(pole, &LG_GRAVITY_EARTH);
    ASSERT_NEAR(accel_pole.x, 0.0f, 1e-6f);
    ASSERT_NEAR(accel_pole.y, 0.0f, 1e-6f);
    /* J2 z-component at pole opposes central gravity (positive = away from Earth) */
    ASSERT_TRUE(accel_pole.z > 0.0f);
    
    printf("PASS: j2_perturbation\n");
    return 0;
}

/*==========================================================================
 * J3/J4 Perturbation
 *==========================================================================*/

int test_j3_j4_perturbation(void) {
    double r = LG_RADIUS_EARTH + 400e3;
    lg_vec3_t pos = { (float)r, 0.0f, 100e3f }; /* slightly north of equator */
    
    lg_vec3_t j2 = lg_gravity_j2(pos, &LG_GRAVITY_EARTH);
    lg_vec3_t j3 = lg_gravity_j3(pos, &LG_GRAVITY_EARTH);
    lg_vec3_t j4 = lg_gravity_j4(pos, &LG_GRAVITY_EARTH);
    lg_vec3_t total = lg_gravity_with_zonal(pos, &LG_GRAVITY_EARTH);
    
    /* J3/J4 should be much smaller than J2 and point */
    ASSERT_TRUE(fabs(j3.x) < fabs(j2.x) * 0.01f);
    ASSERT_TRUE(fabs(j4.x) < fabs(j2.x) * 0.01f);
    
    /* Total should equal sum of components approximately */
    lg_vec3_t sum = lg_vec3_add(lg_vec3_add(
        lg_gravity_accel(pos, lg_vec3_zero(), LG_MU_EARTH), j2), j3);
    sum = lg_vec3_add(sum, j4);
    ASSERT_NEAR(total.x, sum.x, 1e-6f);
    ASSERT_NEAR(total.y, sum.y, 1e-6f);
    ASSERT_NEAR(total.z, sum.z, 1e-6f);
    
    printf("PASS: j3_j4_perturbation\n");
    return 0;
}

/*==========================================================================
 * Atmospheric Drag
 *==========================================================================*/

int test_atmospheric_drag(void) {
    /* Simple exponential atmosphere */
    lg_atmosphere_t atm = {
        .rho0 = 1.225,
        .h0 = 0.0,
        .H = 8500.0,
        .body_radius = LG_RADIUS_EARTH
    };
    
    lg_vec3_t pos = { (float)(LG_RADIUS_EARTH + 100e3), 0.0f, 0.0f }; /* 100 km altitude */
    lg_vec3_t vel = { 0.0f, 7800.0f, 0.0f }; /* LEO velocity */
    
    lg_vec3_t drag = lg_gravity_drag(pos, vel, 2.2, 10.0, 1000.0, &atm);
    
    /* Drag should oppose velocity (negative y) */
    ASSERT_TRUE(drag.y < 0.0f);
    ASSERT_NEAR(drag.x, 0.0f, 1e-6f);
    ASSERT_NEAR(drag.z, 0.0f, 1e-6f);
    
    /* At 100 km drag is significant (~0.01 m/s²) */
    ASSERT_TRUE(drag.y < -1e-3f);
    
    /* Zero velocity -> zero drag */
    lg_vec3_t zero_vel = { 0.0f, 0.0f, 0.0f };
    lg_vec3_t no_drag = lg_gravity_drag(pos, zero_vel, 2.2, 10.0, 1000.0, &atm);
    ASSERT_VEC3_NEAR(no_drag, 0.0f, 0.0f, 0.0f, 1e-6f);
    
    printf("PASS: atmospheric_drag\n");
    return 0;
}

/*==========================================================================
 * Solar Radiation Pressure
 *==========================================================================*/

int test_solar_radiation_pressure(void) {
    /* Vector from Sun to satellite: satellite is 1 AU from Sun along +x */
    lg_vec3_t sun_to_sat = { (float)LG_AU, 0.0f, 0.0f };
    
    lg_vec3_t srp = lg_gravity_srp(sun_to_sat, 1.3, 10.0, 1000.0);
    
    /* Force points away from Sun (+x direction) */
    ASSERT_TRUE(srp.x > 0.0f);
    ASSERT_NEAR(srp.y, 0.0f, 1e-10f);
    ASSERT_NEAR(srp.z, 0.0f, 1e-10f);
    
    /* At 1 AU for 10 m², 1000 kg, Cr=1.3: */
    /* a = 1.3 * 10 * 4.56e-6 / 1000 ≈ 5.93e-8 m/s² */
    ASSERT_NEAR_D((double)srp.x, 5.93e-8, 1e-9);
    
    printf("PASS: solar_radiation_pressure\n");
    return 0;
}

/*==========================================================================
 * Third-Body Perturbation
 *==========================================================================*/

int test_third_body(void) {
    /* Satellite at LEO, Moon along x-axis at ~384,400 km */
    lg_vec3_t sat = { (float)(LG_RADIUS_EARTH + 400e3), 0.0f, 0.0f };
    lg_vec3_t moon = { 3.844e8f, 0.0f, 0.0f };
    
    lg_vec3_t pert = lg_gravity_third_body(sat, moon, LG_MU_MOON);
    
    /* Perturbation should be small but non-zero, pulling toward Moon (+x) */
    ASSERT_TRUE(pert.x > 0.0f);
    ASSERT_NEAR(pert.y, 0.0f, 1e-10f);
    ASSERT_NEAR(pert.z, 0.0f, 1e-10f);
    
    printf("PASS: third_body\n");
    return 0;
}

/*==========================================================================
 * Kepler Solver
 *==========================================================================*/

int test_kepler_solver(void) {
    double e = 0.5;
    double M = 1.0;
    double E = lg_kepler_solve(M, e, 1e-12, 50);
    
    /* Verify Kepler's equation */
    ASSERT_NEAR_D(E - e * sin(E), M, 1e-10);
    
    /* For e=0, E=M */
    double E_circ = lg_kepler_solve(M, 0.0, 1e-12, 50);
    ASSERT_NEAR_D(E_circ, M, 1e-12);
    
    printf("PASS: kepler_solver\n");
    return 0;
}

/*==========================================================================
 * Mean Motion & Mean Anomaly
 *==========================================================================*/

int test_mean_motion_anomaly(void) {
    double a = 7e6;
    double n = lg_mean_motion(LG_MU_EARTH, a);
    double T = lg_orbital_period(LG_MU_EARTH, a);
    ASSERT_NEAR_D(n, 2.0 * LG_PI / T, 1e-10);
    
    double M = lg_mean_anomaly(n, T / 4.0);
    ASSERT_NEAR_D(M, LG_PI / 2.0, 1e-10);
    
    printf("PASS: mean_motion_anomaly\n");
    return 0;
}

/*==========================================================================
 * Flight Path Angle
 *==========================================================================*/

int test_flight_path_angle(void) {
    /* Circular orbit: flight path angle = 0 everywhere */
    ASSERT_NEAR_D(lg_flight_path_angle(0.0, 0.0), 0.0, 1e-12);
    ASSERT_NEAR_D(lg_flight_path_angle(0.0, LG_PI / 2.0), 0.0, 1e-12);
    
    /* Elliptical orbit: max flight path angle at ν where tan(γ) = e*sin(ν)/(1+e*cos(ν)) */
    double e = 0.5;
    double gamma = lg_flight_path_angle(e, LG_PI / 2.0);
    ASSERT_NEAR_D(gamma, atan2(e, 1.0), 1e-12);
    
    printf("PASS: flight_path_angle\n");
    return 0;
}

/*==========================================================================
 * Hill Sphere
 *==========================================================================*/

int test_hill_sphere(void) {
    /* Earth-Moon system: a ~ 384,400 km, e ~ 0.0549 */
    /* Moon's Hill sphere around Earth: a=384400 km, m=Moon, M=Earth */
    double hill = lg_hill_sphere(3.844e8, 0.0549, 7.342e22, 5.972e24);
    /* Expected ~ 58,000 km */
    ASSERT_TRUE(hill > 5.0e7);
    ASSERT_TRUE(hill < 7.0e7);
    
    printf("PASS: hill_sphere\n");
    return 0;
}

/*==========================================================================
 * Synodic Period
 *==========================================================================*/

int test_synodic_period(void) {
    /* Earth year ~ 365.25 days, Mars year ~ 687 days */
    double T_earth = 365.25;
    double T_mars = 687.0;
    double syn = lg_synodic_period(T_earth, T_mars);
    /* Expected ~ 780 days for Earth-Mars opposition */
    ASSERT_NEAR_D(syn, 779.9, 1.0);
    
    printf("PASS: synodic_period\n");
    return 0;
}

/*==========================================================================
 * Gauss Orbit Determination
 *==========================================================================*/

int test_gauss_orbit_determination(void) {
    /* Known orbit: circular at 7000 km, 30 deg inclination */
    double mu = LG_MU_EARTH;
    double a = 7e6;
    double v_circ = sqrt(mu / a);
    double n_sat = sqrt(mu / (a*a*a));
    double inc_sat = 30.0 * LG_PI / 180.0;
    double raan_sat = 0.0;
    double arg_lat0 = -0.5;
    
    /* Ground observer at Earth's surface */
    double R_earth = LG_RADIUS_EARTH;
    double omega_earth = 7.2921159e-5;
    double obs_lat = 0.3;
    double obs_lon0 = 0.0;
    
    /* Three observation times spanning ~5 minutes */
    double dt = 150.0;
    
    lg_optical_obs_t obs[3];
    for (int i = 0; i < 3; i++) {
        double t = (i - 1) * dt;
        double arg_lat = arg_lat0 + n_sat * t;
        
        lg_vec3d_t sat = {
            a * (cos(arg_lat) * cos(raan_sat) - sin(arg_lat) * cos(inc_sat) * sin(raan_sat)),
            a * (cos(arg_lat) * sin(raan_sat) + sin(arg_lat) * cos(inc_sat) * cos(raan_sat)),
            a * sin(arg_lat) * sin(inc_sat)
        };
        
        double lon = obs_lon0 + omega_earth * t;
        lg_vec3d_t observer = {
            R_earth * cos(obs_lat) * cos(lon),
            R_earth * cos(obs_lat) * sin(lon),
            R_earth * sin(obs_lat)
        };
        
        lg_vec3d_t los = lg_vec3d_sub(sat, observer);
        los = lg_vec3d_norm(los);
        
        obs[i].L = lg_vec3_from_d(los);
        obs[i].R = lg_vec3_from_d(observer);
        obs[i].t = t;
    }
    
    lg_vec3_t recovered_pos, recovered_vel;
    int ok = lg_gauss_orbit_determination(obs, mu, &recovered_pos, &recovered_vel);
    ASSERT_TRUE(ok);
    
    /* True state at t=0 */
    lg_vec3d_t true_sat = {
        a * (cos(arg_lat0) * cos(raan_sat) - sin(arg_lat0) * cos(inc_sat) * sin(raan_sat)),
        a * (cos(arg_lat0) * sin(raan_sat) + sin(arg_lat0) * cos(inc_sat) * cos(raan_sat)),
        a * sin(arg_lat0) * sin(inc_sat)
    };
    lg_vec3d_t true_vel_d = {
        -v_circ * (sin(arg_lat0) * cos(raan_sat) + cos(arg_lat0) * cos(inc_sat) * sin(raan_sat)),
        -v_circ * (sin(arg_lat0) * sin(raan_sat) - cos(arg_lat0) * cos(inc_sat) * cos(raan_sat)),
        v_circ * cos(arg_lat0) * sin(inc_sat)
    };
    
    lg_vec3_t true_pos_f = lg_vec3_from_d(true_sat);
    lg_vec3_t true_vel_f = lg_vec3_from_d(true_vel_d);
    
    /* Gauss method should recover position within ~5% and velocity within ~20% */
    float pos_tol = a * 0.05f;
    float vel_tol = v_circ * 0.2f;
    ASSERT_NEAR(recovered_pos.x, true_pos_f.x, pos_tol);
    ASSERT_NEAR(recovered_pos.y, true_pos_f.y, pos_tol);
    ASSERT_NEAR(recovered_pos.z, true_pos_f.z, pos_tol);
    ASSERT_NEAR(recovered_vel.x, true_vel_f.x, vel_tol);
    ASSERT_NEAR(recovered_vel.y, true_vel_f.y, vel_tol);
    ASSERT_NEAR(recovered_vel.z, true_vel_f.z, vel_tol);
    
    printf("PASS: gauss_orbit_determination\n");
    return 0;
}

/*==========================================================================
 * Main
 *==========================================================================*/

typedef int (*test_func_t)(void);

typedef struct {
    const char* name;
    test_func_t func;
} test_t;

int main(void) {
    printf("=== Lagrange Gravity Tests ===\n\n");
    
    test_t tests[] = {
        {"j2_perturbation", test_j2_perturbation},
        {"j3_j4_perturbation", test_j3_j4_perturbation},
        {"atmospheric_drag", test_atmospheric_drag},
        {"solar_radiation_pressure", test_solar_radiation_pressure},
        {"third_body", test_third_body},
        {"kepler_solver", test_kepler_solver},
        {"mean_motion_anomaly", test_mean_motion_anomaly},
        {"flight_path_angle", test_flight_path_angle},
        {"hill_sphere", test_hill_sphere},
        {"synodic_period", test_synodic_period},
        {"gauss_orbit_determination", test_gauss_orbit_determination},
    };
    
    int num_tests = sizeof(tests) / sizeof(tests[0]);
    int passed = 0;
    int failed = 0;
    
    for (int i = 0; i < num_tests; i++) {
        int result = tests[i].func();
        if (result == 0) {
            passed++;
        } else {
            failed++;
            printf("FAILED: %s\n", tests[i].name);
        }
    }
    
    printf("\n=== Results ===\n");
    printf("Passed: %d/%d\n", passed, num_tests);
    printf("Failed: %d/%d\n", failed, num_tests);
    
    return failed > 0 ? 1 : 0;
}
