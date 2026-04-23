/*
 * Lagrange Tests - Attitude Control Gaps
 */

#include <stdio.h>
#include <math.h>
#include <stdbool.h>

#define LAGRANGE_IMPLEMENTATION
#include "lagrange.h"

#define EPSILON 1e-3f

#define ASSERT_NEAR(a, b, eps) do { \
    float _a = (a), _b = (b); \
    if (fabsf(_a - _b) > (eps)) { \
        printf("FAIL: %s:%d: Expected %f, got %f (diff: %e)\n", \
               __FILE__, __LINE__, _b, _a, fabsf(_a - _b)); \
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

/*============================================================================
 * Reaction Wheel Tests
 *===========================================================================*/

int test_rw_update(void) {
    lg_vec3_t axis = {0.0f, 0.0f, 1.0f};
    lg_reaction_wheel_t rw = lg_rw_init(axis, 0.5f, 6000.0f);
    
    ASSERT_NEAR(rw.omega, 0.0f, EPSILON);
    ASSERT_NEAR(rw.h, 0.0f, EPSILON);
    
    /* Apply torque command */
    rw.tau_command = 0.1f;
    lg_rw_update(&rw, 1.0f);
    
    /* omega = tau/J * dt = 0.1/0.5 * 1 = 0.2 rad/s */
    ASSERT_NEAR(rw.omega, 0.2f, EPSILON);
    ASSERT_NEAR(rw.h, 0.1f, EPSILON);
    
    /* Saturation */
    rw.tau_command = 1000.0f;
    lg_rw_update(&rw, 1000.0f);
    ASSERT_TRUE(fabsf(rw.omega) <= rw.max_omega * 1.001f);
    ASSERT_TRUE(rw.momentum_saturation <= 1.0f);
    
    printf("PASS: rw_update\n");
    return 0;
}

int test_rw_array_momentum(void) {
    lg_rw_array_t arr = lg_rw_array_triad();
    ASSERT_TRUE(arr.n_wheels == 3);
    
    /* Command some torques and update */
    arr.wheels[0].tau_command = 0.1f;
    arr.wheels[1].tau_command = -0.05f;
    arr.wheels[2].tau_command = 0.1f; /* clamped to max_torque=0.1f */
    lg_rw_array_update(&arr, 1.0f);
    
    lg_vec3_t h = lg_rw_array_momentum(&arr);
    /* h = J * omega = J * (tau/J * dt) = tau * dt */
    ASSERT_NEAR(h.x, 0.1f, EPSILON);
    ASSERT_NEAR(h.y, -0.05f, EPSILON);
    ASSERT_NEAR(h.z, 0.1f, EPSILON);
    
    /* Body torque should oppose motor torque */
    lg_vec3_t tau_body = lg_rw_body_torque(&arr.wheels[0]);
    ASSERT_NEAR(tau_body.x, -0.1f, EPSILON);
    
    printf("PASS: rw_array_momentum\n");
    return 0;
}

/*============================================================================
 * CMG Tests
 *===========================================================================*/

int test_cmg_update(void) {
    lg_vec3_t spin = {1.0f, 0.0f, 0.0f};
    lg_vec3_t gimbal = {0.0f, 0.0f, 1.0f};
    lg_cmg_t cmg = lg_cmg_init(spin, gimbal, 10.0f, 57.2958f); /* 57.3 deg/s = 1 rad/s */
    
    cmg.gimbal_rate_cmd = 1.0f; /* rad/s */
    lg_cmg_update(&cmg, 0.5f);
    
    ASSERT_NEAR(cmg.gimbal_angle, 0.5f, EPSILON);
    ASSERT_NEAR(cmg.gimbal_rate, 1.0f, EPSILON);
    
    /* Spin axis should have rotated about z by 0.5 rad */
    /* New spin axis ~ (cos(0.5), sin(0.5), 0) */
    ASSERT_NEAR(cmg.spin_axis.x, cosf(0.5f), EPSILON);
    ASSERT_NEAR(cmg.spin_axis.y, sinf(0.5f), EPSILON);
    ASSERT_NEAR(cmg.spin_axis.z, 0.0f, EPSILON);
    
    printf("PASS: cmg_update\n");
    return 0;
}

int test_cmg_torque(void) {
    lg_vec3_t spin = {1.0f, 0.0f, 0.0f};
    lg_vec3_t gimbal = {0.0f, 0.0f, 1.0f};
    lg_cmg_t cmg = lg_cmg_init(spin, gimbal, 10.0f, 57.2958f); /* 1 rad/s max gimbal rate */
    
    /* Set gimbal rate directly so spin axis stays exactly aligned with x */
    cmg.gimbal_rate = 1.0f;
    
    lg_vec3_t torque = lg_cmg_output_torque(&cmg);
    /* tau = gimbal_rate * (gimbal_axis x h_vec) */
    /* gimbal = (0,0,1), h_vec = (10,0,0) */
    /* cross((0,0,1), (10,0,0)) = (0, 10, 0) */
    ASSERT_NEAR(torque.x, 0.0f, 0.1f);
    ASSERT_NEAR(torque.y, 10.0f, 0.5f);
    ASSERT_NEAR(torque.z, 0.0f, 0.1f);
    
    /* Array torque */
    lg_cmg_array_t arr = lg_cmg_array_iss_like();
    lg_cmg_compute_jacobian(&arr);
    ASSERT_TRUE(arr.n_cmgs == 4);
    
    lg_vec3_t arr_torque = lg_cmg_array_torque(&arr);
    /* With zero gimbal rates, total torque should be zero */
    ASSERT_NEAR(arr_torque.x, 0.0f, EPSILON);
    ASSERT_NEAR(arr_torque.y, 0.0f, EPSILON);
    ASSERT_NEAR(arr_torque.z, 0.0f, EPSILON);
    
    printf("PASS: cmg_torque\n");
    return 0;
}

/*============================================================================
 * Gyroscopic Coupling
 *===========================================================================*/

int test_gyroscopic_coupling(void) {
    lg_vec3_t omega = {0.0f, 0.0f, 1.0f};
    lg_vec3_t h_wheels = {1.0f, 0.0f, 0.0f};
    
    lg_vec3_t tau = lg_gyroscopic_coupling(&omega, &h_wheels);
    /* cross(omega, h) = (0, 1, 0) */
    ASSERT_VEC3_NEAR(tau, 0.0f, 1.0f, 0.0f, EPSILON);
    
    /* Euler gyroscopic with zero external torque */
    lg_body_t body = lg_body(1.0f);
    body.inertia = lg_vec3(1.0f, 2.0f, 3.0f);
    body.inv_inertia = lg_vec3(1.0f, 0.5f, 0.333333f);
    body.angular_velocity = omega;
    
    lg_vec3_t alpha = lg_euler_gyroscopic(&body, &h_wheels, &(lg_vec3_t){0,0,0});
    /* I*omega = (0,0,3), h_total = (1,0,3), gyro = cross((0,0,1), (1,0,3)) = (0,1,0) */
    /* net_tau = -gyro = (0,-1,0), alpha = inv(I)*net_tau = (0, -0.5, 0) */
    ASSERT_NEAR(alpha.x, 0.0f, EPSILON);
    ASSERT_NEAR(alpha.y, -0.5f, EPSILON);
    ASSERT_NEAR(alpha.z, 0.0f, EPSILON);
    
    printf("PASS: gyroscopic_coupling\n");
    return 0;
}

/*============================================================================
 * Precession
 *===========================================================================*/

int test_precession(void) {
    lg_precession_state_t p = lg_precession_init(1.0f, 2.0f, 1.0f, 0.1f);
    
    /* Nutation angle causes wobble; z-component of omega should be preserved */
    float wz0 = p.omega_body.z;
    lg_precession_step(&p, 1.0f);
    ASSERT_NEAR(p.omega_body.z, wz0, EPSILON);
    
    /* Attitude quaternion should remain normalized */
    float qnorm = sqrtf(p.attitude.x*p.attitude.x + p.attitude.y*p.attitude.y +
                        p.attitude.z*p.attitude.z + p.attitude.w*p.attitude.w);
    ASSERT_NEAR(qnorm, 1.0f, EPSILON);
    
    printf("PASS: precession\n");
    return 0;
}

/*============================================================================
 * Main
 *===========================================================================*/

typedef int (*test_func_t)(void);

typedef struct {
    const char* name;
    test_func_t func;
} test_t;

int main(void) {
    printf("=== Lagrange Attitude Tests ===\n\n");
    
    test_t tests[] = {
        {"rw_update", test_rw_update},
        {"rw_array_momentum", test_rw_array_momentum},
        {"cmg_update", test_cmg_update},
        {"cmg_torque", test_cmg_torque},
        {"gyroscopic_coupling", test_gyroscopic_coupling},
        {"precession", test_precession},
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
