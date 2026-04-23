/*
 * Lagrange Tests - Koopman Operator Module
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
 * Observable Tests
 *===========================================================================*/

int test_koopman_observable(void) {
    lg_vec3_t pos = {1.0f, 2.0f, 3.0f};
    lg_vec3_t vel = {4.0f, 5.0f, 6.0f};
    float omega = 1.0f;
    
    /* Mode 0 is linear in pos->x */
    float complex psi0 = lg_koop_observable_dct(&pos, &vel, omega, 0);
    ASSERT_NEAR(crealf(psi0), 1.0f, EPSILON);
    ASSERT_NEAR(cimagf(psi0), 0.0f, EPSILON);
    
    /* Mode 1 is linear in pos->y */
    float complex psi1 = lg_koop_observable_dct(&pos, &vel, omega, 1);
    ASSERT_NEAR(crealf(psi1), 2.0f, EPSILON);
    ASSERT_NEAR(cimagf(psi1), 0.0f, EPSILON);
    
    /* Mode 3 is linear in vel->x */
    float complex psi3 = lg_koop_observable_dct(&pos, &vel, omega, 3);
    ASSERT_NEAR(crealf(psi3), 4.0f, EPSILON);
    ASSERT_NEAR(cimagf(psi3), 0.0f, EPSILON);
    
    /* Mode 6+ is spectral (complex exponential) */
    float complex psi6 = lg_koop_observable_dct(&pos, &vel, omega, 6);
    ASSERT_NEAR(cabsf(psi6), 1.0f, EPSILON);
    float phase6 = cargf(psi6);
    while (phase6 > 1.0f + M_PI) phase6 -= 2.0f * M_PI;
    while (phase6 < 1.0f - M_PI) phase6 += 2.0f * M_PI;
    ASSERT_NEAR(phase6, 1.0f, EPSILON); /* phase = omega * pos->x */
    
    printf("PASS: koopman_observable\n");
    return 0;
}

/*============================================================================
 * Fit / Predict Tests
 *===========================================================================*/

int test_koopman_fit_predict(void) {
    /* Three independent harmonic oscillators (rich 6D state)
     * Exact discrete map with dt = 0.1 for each pair:
     * [q ; vq]_{n+1} = [cos(w*dt), sin(w*dt); -sin(w*dt), cos(w*dt)] * [q ; vq]_n
     */
    float dt = 0.1f;
    float wx = 1.0f, wy = 1.5f, wz = 0.7f;
    float cx = cosf(wx * dt), sx = sinf(wx * dt);
    float cy = cosf(wy * dt), sy = sinf(wy * dt);
    float cz = cosf(wz * dt), sz = sinf(wz * dt);
    
    int n_snapshots = 150;
    int n_states = 6;
    
    float* X = (float*)calloc(n_states * n_snapshots, sizeof(float));
    float* Y = (float*)calloc(n_states * n_snapshots, sizeof(float));
    
    /* Generate training data from one 6D trajectory */
    float x = 1.0f, vx = 0.3f;
    float y = 0.5f, vy = 0.2f;
    float z = 0.3f, vz = 0.1f;
    for (int n = 0; n < n_snapshots; n++) {
        X[0 + n*n_states] = x;
        X[1 + n*n_states] = y;
        X[2 + n*n_states] = z;
        X[3 + n*n_states] = vx;
        X[4 + n*n_states] = vy;
        X[5 + n*n_states] = vz;
        
        float x_next = cx * x + sx * vx;
        float vx_next = -sx * x + cx * vx;
        float y_next = cy * y + sy * vy;
        float vy_next = -sy * y + cy * vy;
        float z_next = cz * z + sz * vz;
        float vz_next = -sz * z + cz * vz;
        
        Y[0 + n*n_states] = x_next;
        Y[1 + n*n_states] = y_next;
        Y[2 + n*n_states] = z_next;
        Y[3 + n*n_states] = vx_next;
        Y[4 + n*n_states] = vy_next;
        Y[5 + n*n_states] = vz_next;
        
        x = x_next + 0.03f * sinf(n * 0.7f);
        vx = vx_next + 0.03f * cosf(n * 0.7f);
        y = y_next + 0.02f * sinf(n * 0.5f + 1.0f);
        vy = vy_next + 0.02f * cosf(n * 0.5f + 1.0f);
        z = z_next + 0.01f * sinf(n * 0.3f + 2.0f);
        vz = vz_next + 0.01f * cosf(n * 0.3f + 2.0f);
    }
    
    lg_koop_data_t data = {
        .X = X,
        .Y = Y,
        .n_snapshots = n_snapshots,
        .n_states = n_states
    };
    
    lg_koopman_t km = {0};
    lg_koopman_init_cheap_compatible(&km, 6, 0.5f);
    
    lg_koopman_fit(&km, &data);
    
    /* Verify K and P were allocated and have finite values */
    ASSERT_TRUE(km.K != NULL);
    ASSERT_TRUE(km.P != NULL);
    
    bool k_finite = true, p_finite = true;
    for (int i = 0; i < km.n_observables * km.n_observables; i++) {
        if (isnan(crealf(km.K[i])) || isnan(cimagf(km.K[i]))) k_finite = false;
    }
    for (int i = 0; i < km.n_states * km.n_observables; i++) {
        if (isnan(crealf(km.P[i])) || isnan(cimagf(km.P[i]))) p_finite = false;
    }
    ASSERT_TRUE(k_finite);
    ASSERT_TRUE(p_finite);
    
    /* Predict one step from a test state */
    lg_vec3_t pos_in = {1.0f, 0.5f, 0.3f};
    lg_vec3_t vel_in = {0.3f, 0.2f, 0.1f};
    lg_vec3_t pos_out, vel_out;
    
    lg_koopman_predict(&km, &pos_in, &vel_in, &pos_out, &vel_out, 1);
    
    float x_expected = cx * 1.0f + sx * 0.3f;
    float vx_expected = -sx * 1.0f + cx * 0.3f;
    float y_expected = cy * 0.5f + sy * 0.2f;
    float vy_expected = -sy * 0.5f + cy * 0.2f;
    
    /* With 6 linear modes and sufficient data, EDMD should recover linear dynamics */
    ASSERT_NEAR(pos_out.x, x_expected, 0.05f);
    ASSERT_NEAR(vel_out.x, vx_expected, 0.05f);
    ASSERT_NEAR(pos_out.y, y_expected, 0.05f);
    ASSERT_NEAR(vel_out.y, vy_expected, 0.05f);
    ASSERT_NEAR(pos_out.z, cz * 0.3f + sz * 0.1f, 0.05f);
    
    /* Multi-step prediction should not diverge wildly */
    lg_koopman_predict(&km, &pos_in, &vel_in, &pos_out, &vel_out, 10);
    ASSERT_TRUE(fabsf(pos_out.x) < 5.0f); /* Should stay bounded for oscillator */
    ASSERT_TRUE(fabsf(vel_out.x) < 5.0f);
    
    free(X);
    free(Y);
    free(km.frequencies);
    free(km.scales);
    free(km.K);
    free(km.P);
    free(km.Psi);
    
    printf("PASS: koopman_fit_predict\n");
    return 0;
}

/*============================================================================
 * CHEAP Compatibility Test
 *===========================================================================*/

int test_koopman_cheap_init(void) {
    lg_koopman_t km = {0};
    lg_koopman_init_cheap_compatible(&km, 16, 0.3f);
    
    ASSERT_TRUE(km.type == LG_KOOP_DCT);
    ASSERT_TRUE(km.n_observables == 16);
    ASSERT_TRUE(km.n_states == 6);
    ASSERT_NEAR(km.hurst, 0.3f, EPSILON);
    
    /* Frequencies should follow DCT-II grid */
    ASSERT_NEAR(km.frequencies[0], M_PI * 0.5f / 16.0f, EPSILON);
    ASSERT_NEAR(km.frequencies[15], M_PI * 15.5f / 16.0f, EPSILON);
    
    free(km.frequencies);
    free(km.scales);
    
    printf("PASS: koopman_cheap_init\n");
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
    printf("=== Lagrange Koopman Tests ===\n\n");
    
    test_t tests[] = {
        {"koopman_observable", test_koopman_observable},
        {"koopman_cheap_init", test_koopman_cheap_init},
        {"koopman_fit_predict", test_koopman_fit_predict},
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
