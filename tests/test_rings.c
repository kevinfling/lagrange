/*
 * Lagrange Tests - Planetary Ring System Module
 */

#include <stdio.h>
#include <math.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

#define LAGRANGE_IMPLEMENTATION
#include "lagrange.h"
#include "lagrange/rings.h"

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

/* Helper: create a minimal ring system without lg_node_t dependency */
static lg_ring_system_t* make_test_ring_system(void) {
    lg_ring_system_t* rings = (lg_ring_system_t*)calloc(1, sizeof(lg_ring_system_t));
    strcpy(rings->name, "Test");
    rings->planet_mass = 5.683e26f;
    rings->planet_radius = 5.8232e7f;
    rings->j2 = 0.016f;
    rings->j4 = -0.001f;
    rings->self_gravity_enabled = true;
    rings->zones_capacity = 16;
    rings->zones = (lg_ring_zone_t*)calloc(rings->zones_capacity, sizeof(lg_ring_zone_t));
    return rings;
}

/*============================================================================
 * Saturn Preset
 *===========================================================================*/

int test_ring_preset_saturn(void) {
    lg_ring_system_t* rings = make_test_ring_system();
    lg_ring_preset_saturn_main(rings);
    
    /* Should have D, C, B, Cassini, A, Encke, Keeler, F, G, E = 10 zones */
    ASSERT_TRUE(rings->n_zones >= 9);
    
    /* B ring should be densest */
    int b_ring_idx = -1;
    float max_tau = 0.0f;
    for (int i = 0; i < rings->n_zones; i++) {
        if (rings->zones[i].tau_normal > max_tau) {
            max_tau = rings->zones[i].tau_normal;
            b_ring_idx = i;
        }
    }
    ASSERT_TRUE(b_ring_idx >= 0);
    ASSERT_NEAR(max_tau, 1.5f, 0.1f);
    
    /* Check Encke Gap exists (~133.4-133.8e6 m) */
    bool has_encke = false;
    for (int i = 0; i < rings->n_zones; i++) {
        if (rings->zones[i].r_inner > 133000e3f && rings->zones[i].r_inner < 134000e3f &&
            rings->zones[i].type == LG_ZONE_GAP) {
            has_encke = true;
        }
    }
    ASSERT_TRUE(has_encke);
    
    /* Check G ring exists */
    bool has_g = false;
    for (int i = 0; i < rings->n_zones; i++) {
        if (rings->zones[i].r_inner > 160000e3f && rings->zones[i].r_inner < 170000e3f) {
            has_g = true;
        }
    }
    ASSERT_TRUE(has_g);
    
    free(rings->zones);
    free(rings);
    
    printf("PASS: ring_preset_saturn\n");
    return 0;
}

/*============================================================================
 * Particle Population
 *===========================================================================*/

int test_ring_zone_populate(void) {
    lg_ring_system_t* rings = make_test_ring_system();
    lg_ring_zone_t* zone = lg_ring_system_add_zone(
        rings, 90000e3f, 100000e3f, LG_ZONE_CONTINUOUS,
        LG_RING_ICE_WATER, &LG_DIST_SATURN_ICE);
    zone->tau_normal = 0.5f;
    zone->velocity_dispersion = 0.01f;
    zone->scale_height = 10.0f;
    
    int n_particles = 100;
    lg_ring_zone_populate(zone, rings, n_particles);
    
    ASSERT_TRUE(rings->particles.count == n_particles);
    ASSERT_TRUE(rings->particles.base.count == n_particles);
    ASSERT_TRUE(rings->particles.capacity >= n_particles);
    
    /* Check particles are within zone bounds */
    for (int i = 0; i < n_particles; i++) {
        float x = rings->particles.base.x[i];
        float y = rings->particles.base.y[i];
        float r = sqrtf(x*x + y*y);
        ASSERT_TRUE(r >= zone->r_inner * 0.9f && r <= zone->r_outer * 1.1f);
        ASSERT_TRUE(rings->particles.radius[i] > 0.0f);
        ASSERT_TRUE(rings->particles.mass[i] > 0.0f);
        ASSERT_TRUE(rings->particles.temperature[i] > 0.0f);
    }
    
    /* Clean up particle arrays */
    free(rings->particles.base.x);
    free(rings->particles.base.y);
    free(rings->particles.base.z);
    free(rings->particles.base.vx);
    free(rings->particles.base.vy);
    free(rings->particles.base.vz);
    free(rings->particles.radius);
    free(rings->particles.mass);
    free(rings->particles.composition);
    free(rings->particles.internal_density);
    free(rings->particles.temperature);
    free(rings->particles.charge);
    free(rings->particles.collision_count);
    free(rings->particles.last_collision_time);
    for (int d = 0; d < 3; d++) free(rings->particles.spin[d]);
    
    free(rings->zones);
    free(rings);
    
    printf("PASS: ring_zone_populate\n");
    return 0;
}

/*============================================================================
 * Physics: Toomre Q and Gap Width
 *===========================================================================*/

int test_ring_toomre_q(void) {
    float r = 100000e3f;
    float Omega = sqrtf(5.683e26f * 6.67430e-11f / (r*r*r));
    float sigma = 0.01f;
    float surface_density = 100.0f;
    
    float Q = lg_ring_toomre_Q(r, Omega, sigma, surface_density);
    ASSERT_TRUE(Q > 0.0f);
    ASSERT_TRUE(isfinite(Q));
    
    /* Higher surface density → lower Q (more unstable) */
    float Q2 = lg_ring_toomre_Q(r, Omega, sigma, surface_density * 10.0f);
    ASSERT_TRUE(Q2 < Q);
    
    printf("PASS: ring_toomre_q\n");
    return 0;
}

int test_ring_gap_width(void) {
    /* q = 1e-7, h = 0.01, alpha = 0.01 */
    float width = lg_ring_gap_width(1e-7f, 0.01f, 0.01f);
    ASSERT_TRUE(width > 0.0f);
    ASSERT_TRUE(isfinite(width));
    
    printf("PASS: ring_gap_width\n");
    return 0;
}

/*============================================================================
 * Collisional Evolution
 *===========================================================================*/

int test_ring_collisional_evolution(void) {
    lg_ring_system_t* rings = make_test_ring_system();
    lg_ring_zone_t* zone = lg_ring_system_add_zone(
        rings, 90000e3f, 100000e3f, LG_ZONE_CONTINUOUS,
        LG_RING_ICE_WATER, &LG_DIST_SATURN_ICE);
    zone->tau_normal = 0.5f;
    zone->velocity_dispersion = 0.01f;
    
    lg_ring_zone_populate(zone, rings, 50);
    
    /* Set initial temperatures to 0 */
    for (int i = 0; i < rings->particles.count; i++) {
        rings->particles.temperature[i] = 0.0f;
    }
    
    /* Run collisional evolution */
    rings->age = 0.0;
    lg_ring_collisional_evolution(rings, 1000.0f);
    
    /* Temperatures should have equilibrated upward */
    float T_sum = 0.0f;
    for (int i = 0; i < rings->particles.count; i++) {
        T_sum += rings->particles.temperature[i];
    }
    ASSERT_TRUE(T_sum > 0.0f);
    
    /* Collision counts should have increased */
    float count_sum = 0.0f;
    for (int i = 0; i < rings->particles.count; i++) {
        count_sum += rings->particles.collision_count[i];
    }
    ASSERT_TRUE(count_sum > 0.0f);
    
    /* Cleanup */
    free(rings->particles.base.x);
    free(rings->particles.base.y);
    free(rings->particles.base.z);
    free(rings->particles.base.vx);
    free(rings->particles.base.vy);
    free(rings->particles.base.vz);
    free(rings->particles.radius);
    free(rings->particles.mass);
    free(rings->particles.composition);
    free(rings->particles.internal_density);
    free(rings->particles.temperature);
    free(rings->particles.charge);
    free(rings->particles.collision_count);
    free(rings->particles.last_collision_time);
    for (int d = 0; d < 3; d++) free(rings->particles.spin[d]);
    free(rings->zones);
    free(rings);
    
    printf("PASS: ring_collisional_evolution\n");
    return 0;
}

/*============================================================================
 * Self-Gravity Wakes
 *===========================================================================*/

int test_ring_self_gravity_wakes(void) {
    lg_ring_system_t* rings = make_test_ring_system();
    rings->self_gravity_enabled = true;
    
    /* Create a dense zone that should be unstable (Q < 2) */
    lg_ring_zone_t* zone = lg_ring_system_add_zone(
        rings, 90000e3f, 100000e3f, LG_ZONE_CONTINUOUS,
        LG_RING_ICE_WATER, &LG_DIST_SATURN_ICE);
    zone->tau_normal = 2.0f;
    zone->velocity_dispersion = 0.001f;  /* Very cold → unstable */
    zone->surface_density = 10000.0f;
    zone->scale_height = 1.0f;
    
    float tau_before = zone->tau_normal;
    float sigma_before = zone->velocity_dispersion;
    
    lg_ring_self_gravity_wakes(rings, 100.0f);
    
    /* Unstable zone should form wakes: tau increases, dispersion decreases */
    ASSERT_TRUE(zone->tau_normal >= tau_before);
    ASSERT_TRUE(zone->velocity_dispersion <= sigma_before);
    
    free(rings->zones);
    free(rings);
    
    printf("PASS: ring_self_gravity_wakes\n");
    return 0;
}

/*============================================================================
 * Occultation Profile
 *===========================================================================*/

int test_ring_occultation(void) {
    lg_ring_system_t* rings = make_test_ring_system();
    lg_ring_preset_saturn_main(rings);
    
    float tau[256];
    lg_ring_compute_occultation(rings, tau, 256);
    
    /* Check that profile has variation (gaps should show as dips) */
    float tau_max = 0.0f, tau_min = 1e10f;
    for (int i = 0; i < 256; i++) {
        if (tau[i] > tau_max) tau_max = tau[i];
        if (tau[i] < tau_min) tau_min = tau[i];
    }
    ASSERT_TRUE(tau_max > tau_min);
    ASSERT_TRUE(tau_max > 0.5f);  /* B ring should be dense */
    
    free(rings->zones);
    free(rings);
    
    printf("PASS: ring_occultation\n");
    return 0;
}

/*============================================================================
 * System Update
 *===========================================================================*/

int test_ring_system_update(void) {
    lg_ring_system_t* rings = make_test_ring_system();
    lg_ring_zone_t* zone = lg_ring_system_add_zone(
        rings, 90000e3f, 100000e3f, LG_ZONE_CONTINUOUS,
        LG_RING_ICE_WATER, &LG_DIST_SATURN_ICE);
    zone->tau_normal = 0.5f;
    zone->velocity_dispersion = 0.01f;
    zone->scale_height = 10.0f;
    zone->viscosity_alpha = 0.01f;
    
    lg_ring_zone_populate(zone, rings, 20);
    
    /* Record initial positions */
    float x0 = rings->particles.base.x[0];
    float y0 = rings->particles.base.y[0];
    
    /* Run update */
    rings->age = 0.0;
    lg_ring_system_update(rings, 0.0, 10.0f);
    
    /* Particles should have moved */
    float dx = rings->particles.base.x[0] - x0;
    float dy = rings->particles.base.y[0] - y0;
    ASSERT_TRUE(dx*dx + dy*dy > 1e-6f);
    
    /* Age should have advanced */
    ASSERT_NEAR((float)rings->age, 10.0f, 0.1f);
    
    /* Cleanup */
    free(rings->particles.base.x);
    free(rings->particles.base.y);
    free(rings->particles.base.z);
    free(rings->particles.base.vx);
    free(rings->particles.base.vy);
    free(rings->particles.base.vz);
    free(rings->particles.radius);
    free(rings->particles.mass);
    free(rings->particles.composition);
    free(rings->particles.internal_density);
    free(rings->particles.temperature);
    free(rings->particles.charge);
    free(rings->particles.collision_count);
    free(rings->particles.last_collision_time);
    for (int d = 0; d < 3; d++) free(rings->particles.spin[d]);
    free(rings->zones);
    free(rings);
    
    printf("PASS: ring_system_update\n");
    return 0;
}

/*============================================================================
 * LOD Update
 *===========================================================================*/

int test_ring_lod_update(void) {
    lg_ring_system_t* rings = make_test_ring_system();
    lg_ring_preset_saturn_main(rings);
    
    /* Very close: full particles */
    lg_ring_lod_update(rings, 1e6f, 1e-6f);
    ASSERT_TRUE(rings->n_particles_target == 100000);
    
    /* Medium distance: aggregates */
    lg_ring_lod_update(rings, 1e12f, 1e-6f);
    ASSERT_TRUE(rings->n_particles_target == 10000);
    
    /* Far: texture */
    lg_ring_lod_update(rings, 1e13f, 1e-6f);
    ASSERT_TRUE(rings->n_particles_target == 1000);
    
    /* Very far: billboard */
    lg_ring_lod_update(rings, 1e14f, 1e-6f);
    ASSERT_TRUE(rings->n_particles_target == 0);
    
    free(rings->zones);
    free(rings);
    
    printf("PASS: ring_lod_update\n");
    return 0;
}

/*============================================================================
 * Size Sampling
 *===========================================================================*/

int test_ring_sample_size(void) {
    lg_ring_size_distribution_t dist = LG_DIST_SATURN_ICE;
    
    float a_sum = 0.0f;
    int n = 1000;
    for (int i = 0; i < n; i++) {
        float a = lg_ring_sample_size(&dist);
        ASSERT_TRUE(a >= dist.a_min && a <= dist.a_max);
        a_sum += a;
    }
    
    /* Average should be somewhere in the middle (not all at min or max) */
    float a_avg = a_sum / n;
    ASSERT_TRUE(a_avg > dist.a_min * 2.0f);
    ASSERT_TRUE(a_avg < dist.a_max * 0.5f);
    
    printf("PASS: ring_sample_size\n");
    return 0;
}

/*============================================================================
 * Density Wavelength
 *===========================================================================*/

int test_ring_density_wavelength(void) {
    float r = 100000e3f;
    float Omega = sqrtf(5.683e26f * 6.67430e-11f / (r*r*r));
    float sigma = 0.01f;
    float kappa = Omega;  /* Nearly circular */
    
    float lambda = lg_ring_density_wavelength(r, 2, Omega, sigma, kappa);
    ASSERT_TRUE(lambda > 0.0f);
    ASSERT_TRUE(isfinite(lambda));
    
    /* Evanescent case should return large value */
    float lambda_evan = lg_ring_density_wavelength(r, 100, Omega, sigma, 0.9f * Omega);
    ASSERT_TRUE(lambda_evan > 1e20f);
    
    printf("PASS: ring_density_wavelength\n");
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
    printf("=== Lagrange Rings Tests ===\n\n");
    
    test_t tests[] = {
        {"ring_preset_saturn", test_ring_preset_saturn},
        {"ring_zone_populate", test_ring_zone_populate},
        {"ring_toomre_q", test_ring_toomre_q},
        {"ring_gap_width", test_ring_gap_width},
        {"ring_collisional_evolution", test_ring_collisional_evolution},
        {"ring_self_gravity_wakes", test_ring_self_gravity_wakes},
        {"ring_occultation", test_ring_occultation},
        {"ring_system_update", test_ring_system_update},
        {"ring_lod_update", test_ring_lod_update},
        {"ring_sample_size", test_ring_sample_size},
        {"ring_density_wavelength", test_ring_density_wavelength},
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
