/*
 * Lagrange Tests - Wavelet Module
 */

#include <stdio.h>
#include <math.h>
#include <stdbool.h>
#include <string.h>

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

#define ASSERT_TRUE(cond) do { \
    if (!(cond)) { \
        printf("FAIL: %s:%d: Assertion failed: %s\n", __FILE__, __LINE__, #cond); \
        return 1; \
    } \
} while(0)

static float signal_max_diff(const float* a, const float* b, int n) {
    float max_diff = 0.0f;
    for (int i = 0; i < n; i++) {
        float diff = fabsf(a[i] - b[i]);
        if (diff > max_diff) max_diff = diff;
    }
    return max_diff;
}

/*============================================================================
 * Round-trip Tests
 *===========================================================================*/

int test_wavelet_haar_roundtrip(void) {
    float original[8] = {1.0f, 2.0f, 3.0f, 4.0f, 5.0f, 6.0f, 7.0f, 8.0f};
    float data[8];
    memcpy(data, original, sizeof(original));
    
    lg_wavelet_t w = {.type = LG_WAVELET_HAAR, .n = 8, .levels = 3};
    lg_wavelet_transform(&w, data);
    lg_wavelet_transform_inverse(&w, data);
    
    ASSERT_NEAR(signal_max_diff(original, data, 8), 0.0f, EPSILON);
    
    printf("PASS: wavelet_haar_roundtrip\n");
    return 0;
}

int test_wavelet_d4_roundtrip(void) {
    float original[8] = {1.0f, 0.5f, -0.2f, 0.3f, 0.8f, -0.1f, 0.0f, 0.4f};
    float data[8];
    memcpy(data, original, sizeof(original));
    
    lg_wavelet_t w = {.type = LG_WAVELET_D4, .n = 8, .levels = 2};
    lg_wavelet_transform(&w, data);
    lg_wavelet_transform_inverse(&w, data);
    
    ASSERT_NEAR(signal_max_diff(original, data, 8), 0.0f, 1e-2f);
    
    printf("PASS: wavelet_d4_roundtrip\n");
    return 0;
}

int test_wavelet_cdf97_roundtrip(void) {
    float original[16] = {
        1.0f, 0.8f, 0.6f, 0.4f, 0.2f, 0.0f, -0.2f, -0.4f,
        -0.2f, 0.0f, 0.2f, 0.4f, 0.6f, 0.8f, 1.0f, 1.2f
    };
    float data[16];
    memcpy(data, original, sizeof(original));
    
    lg_wavelet_t w = {.type = LG_WAVELET_CDF97, .n = 16, .levels = 3};
    lg_wavelet_transform(&w, data);
    lg_wavelet_transform_inverse(&w, data);
    
    ASSERT_NEAR(signal_max_diff(original, data, 16), 0.0f, 1e-2f);
    
    printf("PASS: wavelet_cdf97_roundtrip\n");
    return 0;
}

int test_wavelet_dm53_roundtrip(void) {
    float original[8] = {2.0f, 1.5f, 1.0f, 0.5f, 0.0f, 0.5f, 1.0f, 1.5f};
    float data[8];
    memcpy(data, original, sizeof(original));
    
    lg_wavelet_t w = {.type = LG_WAVELET_DM53, .n = 8, .levels = 2};
    lg_wavelet_transform(&w, data);
    lg_wavelet_transform_inverse(&w, data);
    
    ASSERT_NEAR(signal_max_diff(original, data, 8), 0.0f, 1e-2f);
    
    printf("PASS: wavelet_dm53_roundtrip\n");
    return 0;
}

/*============================================================================
 * Orbital Decomposition Test
 *===========================================================================*/

int test_wavelet_orbit_decompose(void) {
    lg_vec3_t positions[8] = {
        {3.0f, 4.0f, 0.0f},  /* magnitude = 5 */
        {0.0f, 1.0f, 0.0f},  /* magnitude = 1 */
        {1.0f, 2.0f, 2.0f},  /* magnitude = 3 */
        {0.0f, 0.0f, 0.0f},  /* magnitude = 0 */
        {2.0f, 0.0f, 0.0f},  /* magnitude = 2 */
        {0.0f, 3.0f, 4.0f},  /* magnitude = 5 */
        {1.0f, 1.0f, 1.0f},  /* magnitude = sqrt(3) */
        {6.0f, 8.0f, 0.0f},  /* magnitude = 10 */
    };
    
    float coeffs[8] = {0};
    lg_wavelet_orbit_t orbit = {
        .coeffs = coeffs,
        .n = 8,
        .n_levels = 3
    };
    
    lg_wavelet_orbit_decompose(&orbit, positions, 8, LG_WAVELET_HAAR);
    
    /* With Haar, the first coefficient after 2 levels should be the average magnitude */
    float avg = (5.0f + 1.0f + 3.0f + 0.0f + 2.0f + 5.0f + 1.73205f + 10.0f) / 8.0f;
    ASSERT_NEAR(coeffs[0], avg, EPSILON);
    
    printf("PASS: wavelet_orbit_decompose\n");
    return 0;
}

/*============================================================================
 * Compression Test
 *===========================================================================*/

int test_wavelet_compress_decompress(void) {
    float original[8] = {1.0f, 1.1f, 1.0f, 1.1f, 1.0f, 1.1f, 1.0f, 1.1f};
    
    lg_wavelet_compressed_t compressed = {0};
    lg_wavelet_compress(original, 8, 0.2f, &compressed, LG_WAVELET_HAAR);
    
    /* A nearly constant signal should have few significant detail coefficients */
    ASSERT_TRUE(compressed.n_nonzero < 8);
    
    float reconstructed[8];
    lg_wavelet_decompress(&compressed, reconstructed, 8);
    
    /* Reconstruction should be close to original for smooth data */
    ASSERT_NEAR(signal_max_diff(original, reconstructed, 8), 0.0f, 0.3f);
    
    free(compressed.indices);
    free(compressed.values);
    
    printf("PASS: wavelet_compress_decompress\n");
    return 0;
}

/*============================================================================
 * Maneuver Detection Test
 *===========================================================================*/

int test_wavelet_detect_maneuver(void) {
    float coeffs[16] = {0};
    /* Inject a large detail coefficient at the finest scale */
    coeffs[8] = 0.0f; /* half = 8, finest detail starts at index 8 */
    coeffs[12] = 5.0f; /* above threshold */
    
    lg_wavelet_orbit_t orbit = {
        .coeffs = coeffs,
        .n = 16,
        .n_levels = 4
    };
    
    int idx;
    bool detected = lg_wavelet_detect_maneuver(&orbit, 1.0f, &idx);
    ASSERT_TRUE(detected);
    ASSERT_TRUE(idx == 4); /* 12 - 8 = 4 */
    
    printf("PASS: wavelet_detect_maneuver\n");
    return 0;
}

/*============================================================================
 * Morlet CWT Test
 *===========================================================================*/

int test_wavelet_morlet(void) {
    float signal[32];
    for (int i = 0; i < 32; i++) {
        signal[i] = sinf(2.0f * M_PI * 0.1f * i); /* 10% of Nyquist */
    }
    
    float out_real[32] = {0};
    float out_imag[32] = {0};
    
    lg_wavelet_morlet_transform(signal, 32, out_real, out_imag, 1, 1.0f, 0.1f);
    
    /* CWT of a sinusoid should have significant energy */
    float energy = 0.0f;
    for (int i = 0; i < 32; i++) {
        energy += out_real[i] * out_real[i] + out_imag[i] * out_imag[i];
    }
    ASSERT_TRUE(energy > 1.0f);
    
    printf("PASS: wavelet_morlet\n");
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
    printf("=== Lagrange Wavelet Tests ===\n\n");
    
    test_t tests[] = {
        {"wavelet_haar_roundtrip", test_wavelet_haar_roundtrip},
        {"wavelet_d4_roundtrip", test_wavelet_d4_roundtrip},
        {"wavelet_cdf97_roundtrip", test_wavelet_cdf97_roundtrip},
        {"wavelet_dm53_roundtrip", test_wavelet_dm53_roundtrip},
        {"wavelet_orbit_decompose", test_wavelet_orbit_decompose},
        {"wavelet_compress_decompress", test_wavelet_compress_decompress},
        {"wavelet_detect_maneuver", test_wavelet_detect_maneuver},
        {"wavelet_morlet", test_wavelet_morlet},
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
