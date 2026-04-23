/*
 * Lagrange Tests - Math Module
 */

#include <stdio.h>
#include <math.h>
#include <stdbool.h>

#define LAGRANGE_IMPLEMENTATION
#include "lagrange.h"

#define EPSILON 1e-5f

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
 * Vector Tests
 *===========================================================================*/

int test_vec3_creation(void) {
    lg_vec3_t v = lg_vec3(1.0f, 2.0f, 3.0f);
    ASSERT_VEC3_NEAR(v, 1.0f, 2.0f, 3.0f, EPSILON);
    
    lg_vec3_t zero = lg_vec3_zero();
    ASSERT_VEC3_NEAR(zero, 0.0f, 0.0f, 0.0f, EPSILON);
    
    lg_vec3_t one = lg_vec3_one();
    ASSERT_VEC3_NEAR(one, 1.0f, 1.0f, 1.0f, EPSILON);
    
    printf("PASS: vec3_creation\n");
    return 0;
}

int test_vec3_arithmetic(void) {
    lg_vec3_t a = lg_vec3(1.0f, 2.0f, 3.0f);
    lg_vec3_t b = lg_vec3(4.0f, 5.0f, 6.0f);
    
    /* Addition */
    lg_vec3_t sum = lg_vec3_add(a, b);
    ASSERT_VEC3_NEAR(sum, 5.0f, 7.0f, 9.0f, EPSILON);
    
    /* Subtraction */
    lg_vec3_t diff = lg_vec3_sub(b, a);
    ASSERT_VEC3_NEAR(diff, 3.0f, 3.0f, 3.0f, EPSILON);
    
    /* Scaling */
    lg_vec3_t scaled = lg_vec3_scale(a, 2.0f);
    ASSERT_VEC3_NEAR(scaled, 2.0f, 4.0f, 6.0f, EPSILON);
    
    /* Negation */
    lg_vec3_t neg = lg_vec3_neg(a);
    ASSERT_VEC3_NEAR(neg, -1.0f, -2.0f, -3.0f, EPSILON);
    
    printf("PASS: vec3_arithmetic\n");
    return 0;
}

int test_vec3_products(void) {
    lg_vec3_t a = lg_vec3(1.0f, 0.0f, 0.0f);
    lg_vec3_t b = lg_vec3(0.0f, 1.0f, 0.0f);
    
    /* Dot product */
    float dot = lg_vec3_dot(a, b);
    ASSERT_NEAR(dot, 0.0f, EPSILON);
    
    dot = lg_vec3_dot(a, a);
    ASSERT_NEAR(dot, 1.0f, EPSILON);
    
    /* Cross product */
    lg_vec3_t cross = lg_vec3_cross(a, b);
    ASSERT_VEC3_NEAR(cross, 0.0f, 0.0f, 1.0f, EPSILON);
    
    cross = lg_vec3_cross(b, a);
    ASSERT_VEC3_NEAR(cross, 0.0f, 0.0f, -1.0f, EPSILON);
    
    printf("PASS: vec3_products\n");
    return 0;
}

int test_vec3_length(void) {
    lg_vec3_t v = lg_vec3(3.0f, 4.0f, 0.0f);
    
    ASSERT_NEAR(lg_vec3_len(v), 5.0f, EPSILON);
    ASSERT_NEAR(lg_vec3_len_sq(v), 25.0f, EPSILON);
    
    /* Normalization */
    lg_vec3_t n = lg_vec3_norm(v);
    ASSERT_NEAR(lg_vec3_len(n), 1.0f, EPSILON);
    ASSERT_VEC3_NEAR(n, 0.6f, 0.8f, 0.0f, EPSILON);
    
    printf("PASS: vec3_length\n");
    return 0;
}

int test_vec3_interpolation(void) {
    lg_vec3_t a = lg_vec3(0.0f, 0.0f, 0.0f);
    lg_vec3_t b = lg_vec3(10.0f, 10.0f, 10.0f);
    
    lg_vec3_t mid = lg_vec3_lerp(a, b, 0.5f);
    ASSERT_VEC3_NEAR(mid, 5.0f, 5.0f, 5.0f, EPSILON);
    
    lg_vec3_t start = lg_vec3_lerp(a, b, 0.0f);
    ASSERT_VEC3_NEAR(start, 0.0f, 0.0f, 0.0f, EPSILON);
    
    lg_vec3_t end = lg_vec3_lerp(a, b, 1.0f);
    ASSERT_VEC3_NEAR(end, 10.0f, 10.0f, 10.0f, EPSILON);
    
    printf("PASS: vec3_interpolation\n");
    return 0;
}

/*============================================================================
 * Quaternion Tests
 *===========================================================================*/

int test_quat_identity(void) {
    lg_quat_t q = lg_quat_identity();
    ASSERT_NEAR(q.x, 0.0f, EPSILON);
    ASSERT_NEAR(q.y, 0.0f, EPSILON);
    ASSERT_NEAR(q.z, 0.0f, EPSILON);
    ASSERT_NEAR(q.w, 1.0f, EPSILON);
    
    printf("PASS: quat_identity\n");
    return 0;
}

int test_quat_rotation(void) {
    /* 90 degree rotation around Z axis */
    lg_quat_t q = lg_quat_from_axis_angle(lg_vec3(0.0f, 0.0f, 1.0f), (float)LG_PI / 2.0f);
    
    lg_vec3_t x_axis = lg_vec3(1.0f, 0.0f, 0.0f);
    lg_vec3_t rotated = lg_quat_rotate(q, x_axis);
    
    /* Should be (0, 1, 0) */
    ASSERT_VEC3_NEAR(rotated, 0.0f, 1.0f, 0.0f, 0.01f);
    
    printf("PASS: quat_rotation\n");
    return 0;
}

int test_quat_from_to(void) {
    lg_vec3_t from = lg_vec3(1.0f, 0.0f, 0.0f);
    lg_vec3_t to = lg_vec3(0.0f, 1.0f, 0.0f);
    
    lg_quat_t q = lg_quat_from_to(from, to);
    lg_vec3_t result = lg_quat_rotate(q, from);
    
    ASSERT_VEC3_NEAR(result, 0.0f, 1.0f, 0.0f, 0.01f);
    
    printf("PASS: quat_from_to\n");
    return 0;
}

/*============================================================================
 * Double Precision Tests
 *===========================================================================*/

int test_vec3d_precision(void) {
    lg_vec3d_t a = lg_vec3d(1e10, 2e10, 3e10);
    lg_vec3d_t b = lg_vec3d(4e10, 5e10, 6e10);
    
    lg_vec3d_t sum = lg_vec3d_add(a, b);
    ASSERT_NEAR((float)sum.x, 5e10f, 1e5f);
    ASSERT_NEAR((float)sum.y, 7e10f, 1e5f);
    ASSERT_NEAR((float)sum.z, 9e10f, 1e5f);
    
    /* Conversion */
    lg_vec3_t f = lg_vec3_from_d(a);
    ASSERT_NEAR(f.x, 1e10f, 1e3f);
    
    printf("PASS: vec3d_precision\n");
    return 0;
}

/*============================================================================
 * Matrix Tests
 *===========================================================================*/

int test_mat3_identity(void) {
    lg_mat3_t m = lg_mat3_identity();
    
    /* Diagonal should be 1 */
    ASSERT_NEAR(m.m[0], 1.0f, EPSILON);
    ASSERT_NEAR(m.m[4], 1.0f, EPSILON);
    ASSERT_NEAR(m.m[8], 1.0f, EPSILON);
    
    /* Off-diagonal should be 0 */
    ASSERT_NEAR(m.m[1], 0.0f, EPSILON);
    ASSERT_NEAR(m.m[2], 0.0f, EPSILON);
    ASSERT_NEAR(m.m[3], 0.0f, EPSILON);
    ASSERT_NEAR(m.m[5], 0.0f, EPSILON);
    ASSERT_NEAR(m.m[6], 0.0f, EPSILON);
    ASSERT_NEAR(m.m[7], 0.0f, EPSILON);
    
    printf("PASS: mat3_identity\n");
    return 0;
}

int test_mat3_determinant(void) {
    /* Identity det = 1 */
    lg_mat3_t i = lg_mat3_identity();
    ASSERT_NEAR(lg_mat3_determinant(i), 1.0f, EPSILON);
    
    /* Simple diagonal matrix det = product of diagonal */
    lg_mat3_t m = {{2.0f, 0.0f, 0.0f, 0.0f, 3.0f, 0.0f, 0.0f, 0.0f, 4.0f}};
    ASSERT_NEAR(lg_mat3_determinant(m), 24.0f, EPSILON);
    
    printf("PASS: mat3_determinant\n");
    return 0;
}

int test_mat3_inverse(void) {
    /* Test that M * M^-1 = I */
    lg_mat3_t m = lg_mat3_from_cols(
        lg_vec3(2.0f, 1.0f, 0.0f),
        lg_vec3(0.0f, 1.0f, 0.0f),
        lg_vec3(0.0f, 0.0f, 1.0f)
    );
    
    lg_mat3_t inv = lg_mat3_inverse(m);
    lg_mat3_t prod = lg_mat3_mul(m, inv);
    
    /* Should be close to identity */
    ASSERT_NEAR(prod.m[0], 1.0f, EPSILON);
    ASSERT_NEAR(prod.m[4], 1.0f, EPSILON);
    ASSERT_NEAR(prod.m[8], 1.0f, EPSILON);
    ASSERT_NEAR(prod.m[1], 0.0f, EPSILON);
    ASSERT_NEAR(prod.m[2], 0.0f, EPSILON);
    ASSERT_NEAR(prod.m[3], 0.0f, EPSILON);
    ASSERT_NEAR(prod.m[5], 0.0f, EPSILON);
    ASSERT_NEAR(prod.m[6], 0.0f, EPSILON);
    ASSERT_NEAR(prod.m[7], 0.0f, EPSILON);
    
    printf("PASS: mat3_inverse\n");
    return 0;
}

int test_mat3_mul_vec3(void) {
    /* Test matrix-vector multiplication */
    lg_mat3_t m = lg_mat3_from_cols(
        lg_vec3(1.0f, 0.0f, 0.0f),
        lg_vec3(0.0f, 2.0f, 0.0f),
        lg_vec3(0.0f, 0.0f, 3.0f)
    );
    
    lg_vec3_t v = lg_vec3(1.0f, 1.0f, 1.0f);
    lg_vec3_t result = lg_mat3_mul_vec3(m, v);
    
    /* Should scale each component */
    ASSERT_VEC3_NEAR(result, 1.0f, 2.0f, 3.0f, EPSILON);
    
    printf("PASS: mat3_mul_vec3\n");
    return 0;
}

int test_quat_slerp_edge_cases(void) {
    lg_quat_t a = lg_quat_from_axis_angle(lg_vec3(0, 0, 1), 0.0f);
    lg_quat_t b = lg_quat_from_axis_angle(lg_vec3(0, 0, 1), LG_PI / 2.0f);
    
    /* t=0 should give a */
    lg_quat_t r0 = lg_quat_slerp(a, b, 0.0f);
    ASSERT_NEAR(r0.w, a.w, EPSILON);
    ASSERT_NEAR(r0.z, a.z, EPSILON);
    
    /* t=1 should give b */
    lg_quat_t r1 = lg_quat_slerp(a, b, 1.0f);
    ASSERT_NEAR(r1.w, b.w, EPSILON);
    ASSERT_NEAR(r1.z, b.z, EPSILON);
    
    /* t=0.5 should be 45 degrees */
    lg_quat_t r05 = lg_quat_slerp(a, b, 0.5f);
    ASSERT_NEAR(r05.w, cosf(LG_PI / 8.0f), EPSILON);
    ASSERT_NEAR(r05.z, sinf(LG_PI / 8.0f), EPSILON);
    
    /* Opposite sign: slerp should take shorter path */
    lg_quat_t b_neg = lg_quat_scale(b, -1.0f);
    lg_quat_t r_neg = lg_quat_slerp(a, b_neg, 0.5f);
    ASSERT_NEAR(fabsf(r_neg.w), fabsf(r05.w), EPSILON);
    ASSERT_NEAR(fabsf(r_neg.z), fabsf(r05.z), EPSILON);
    
    /* Very close quaternions should use lerp fallback (no NaN) */
    lg_quat_t c = lg_quat_norm(lg_quat(1e-4f, 0.0f, 0.0f, 1.0f));
    lg_quat_t d = lg_quat_norm(lg_quat(-1e-4f, 0.0f, 0.0f, 1.0f));
    lg_quat_t r_close = lg_quat_slerp(c, d, 0.5f);
    ASSERT_TRUE(!isnan(r_close.w));
    ASSERT_TRUE(!isnan(r_close.x));
    
    printf("PASS: quat_slerp_edge_cases\n");
    return 0;
}

int test_simd_batch_operations(void) {
    /* Test fast rsqrt approximation */
    float x = 4.0f;
    float inv_sqrt = lg_rsqrt_fast(x);
    ASSERT_NEAR(inv_sqrt, 0.5f, 1e-3f);
    
    /* Test fast normalize */
    lg_vec3_t v = lg_vec3(3.0f, 0.0f, 4.0f);
    lg_vec3_t n = lg_vec3_norm_fast(v);
    ASSERT_NEAR(n.x, 0.6f, 1e-3f);
    ASSERT_NEAR(n.y, 0.0f, 1e-3f);
    ASSERT_NEAR(n.z, 0.8f, 1e-3f);
    
    /* Test batch normalize via generic dispatcher */
    float xs[8] = {3.0f, 1.0f, 0.0f, 2.0f, 0.0f, 0.0f, 1.0f, 0.0f};
    float ys[8] = {4.0f, 0.0f, 1.0f, 0.0f, 2.0f, 0.0f, 0.0f, 1.0f};
    float zs[8] = {0.0f, 0.0f, 0.0f, 1.0f, 0.0f, 3.0f, 0.0f, 0.0f};
    lg_vec3_batch_t batch = {.x = xs, .y = ys, .z = zs, .count = 8, .capacity = 8};
    lg_vec3_batch_normalize(&batch);
    
    for (int i = 0; i < 8; i++) {
        float len_sq = xs[i]*xs[i] + ys[i]*ys[i] + zs[i]*zs[i];
        if (len_sq > 1e-12f) {
            ASSERT_NEAR(sqrtf(len_sq), 1.0f, 1e-3f);
        }
    }
    
    printf("PASS: simd_batch_operations\n");
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
    printf("=== Lagrange Math Tests ===\n\n");
    
    test_t tests[] = {
        {"vec3_creation", test_vec3_creation},
        {"vec3_arithmetic", test_vec3_arithmetic},
        {"vec3_products", test_vec3_products},
        {"vec3_length", test_vec3_length},
        {"vec3_interpolation", test_vec3_interpolation},
        {"quat_identity", test_quat_identity},
        {"quat_rotation", test_quat_rotation},
        {"quat_from_to", test_quat_from_to},
        {"vec3d_precision", test_vec3d_precision},
        {"mat3_identity", test_mat3_identity},
        {"mat3_determinant", test_mat3_determinant},
        {"mat3_inverse", test_mat3_inverse},
        {"mat3_mul_vec3", test_mat3_mul_vec3},
        {"quat_slerp_edge_cases", test_quat_slerp_edge_cases},
        {"simd_batch_operations", test_simd_batch_operations},
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
