/*
 * Lagrange Tests - Orbital Mechanics Gaps
 */

#include <stdio.h>
#include <math.h>
#include <stdbool.h>

#define LAGRANGE_IMPLEMENTATION
#include "lagrange.h"

#define EPSILON 1e-3f
#define EPSILON_LOOSE 1e-2f

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
 * Orbital Elements Round-Trip
 *===========================================================================*/

int test_orbital_elements_roundtrip(void) {
    /* Circular orbit at 7000 km */
    double mu = LG_MU_EARTH;
    double a = 7e6;
    double v_circ = sqrt(mu / a);
    
    lg_vec3_t pos = { (float)a, 0.0f, 0.0f };
    lg_vec3_t vel = { 0.0f, (float)v_circ, 0.0f };
    
    lg_orbital_elements_t elem = lg_state_to_elements(pos, vel, mu);
    
    ASSERT_NEAR((float)elem.semi_major_axis, (float)a, 100.0f);
    ASSERT_NEAR((float)elem.eccentricity, 0.0f, 1e-4f);
    ASSERT_NEAR((float)elem.inclination, 0.0f, 1e-4f);
    
    /* Round-trip back to state vectors */
    lg_vec3_t pos2, vel2;
    lg_elements_to_state(&elem, mu, &pos2, &vel2);
    
    ASSERT_NEAR(pos2.x, pos.x, 1.0f);
    ASSERT_NEAR(pos2.y, pos.y, 1.0f);
    ASSERT_NEAR(vel2.x, vel.x, 0.1f);
    ASSERT_NEAR(vel2.y, vel.y, 0.1f);
    
    printf("PASS: orbital_elements_roundtrip\n");
    return 0;
}

/*============================================================================
 * Hohmann Transfer
 *===========================================================================*/

int test_hohmann_transfer(void) {
    float mu = 3.986004418e14f; /* Earth m^3/s^2 */
    float r1 = 6.6e6f;  /* ~230 km LEO */
    float r2 = 4.216e7f; /* GEO radius */
    
    lg_transfer_solution_t sol = lg_transfer_hohmann(mu, r1, r2);
    
    /* Sanity: total dv should be a few km/s */
    ASSERT_TRUE(sol.dv_total > 3000.0f);
    ASSERT_TRUE(sol.dv_total < 5000.0f);
    
    /* TOF should be about 5 hours for LEO->GEO */
    ASSERT_TRUE(sol.tof > 15000.0f);
    ASSERT_TRUE(sol.tof < 25000.0f);
    
    /* Efficiency should be 1.0 (optimal 2-impulse) */
    ASSERT_NEAR(sol.efficiency, 1.0f, EPSILON);
    
    /* Grid Hohmann baseline */
    float dv_grid[4] = {0};
    lg_transfer_grid_t grid = {
        .r1 = r1, .r2 = r2,
        .departure_times = NULL,
        .tofs = NULL,
        .dv_grid = dv_grid,
        .n_deps = 2, .n_tofs = 2
    };
    lg_transfer_grid_hohmann(&grid, mu);
    ASSERT_NEAR(dv_grid[0], sol.dv_total, EPSILON);
    ASSERT_NEAR(dv_grid[3], sol.dv_total, EPSILON);
    
    printf("PASS: hohmann_transfer\n");
    return 0;
}

/*============================================================================
 * Bi-elliptic Transfer
 *===========================================================================*/

int test_bielliptic_transfer(void) {
    float mu = 3.986004418e14f;
    float r1 = 6.6e6f;
    float r2 = 4.216e7f;
    float rb = 20.0f * r2; /* Very high intermediate */
    
    lg_transfer_solution_t hohm = lg_transfer_hohmann(mu, r1, r2);
    (void)hohm;
    lg_transfer_solution_t biell = lg_transfer_bielliptic(mu, r1, r2, rb);
    
    /* Bi-elliptic should have 3 burns */
    ASSERT_TRUE(biell.n_burns == 3);
    
    /* For LEO->GEO ratio ~6.4, Hohmann usually wins, so bi-elliptic dv > Hohmann */
    /* But function should return finite positive dv */
    ASSERT_TRUE(biell.dv_total > 0.0f);
    ASSERT_TRUE(biell.dv_total < 1e19f);
    
    /* Test the bielliptic_wins helper with an extreme ratio */
    bool wins = lg_transfer_bielliptic_wins(r1, r2, rb);
    ASSERT_TRUE(!wins); /* ratio < 11.94, should not win */
    
    /* Now test with extreme ratio > 11.94 */
    float r_big = r1 * 20.0f;
    lg_transfer_solution_t hohm2 = lg_transfer_hohmann(mu, r1, r_big);
    (void)hohm2;
    lg_transfer_solution_t biell2 = lg_transfer_bielliptic(mu, r1, r_big, r_big * 5.0f);
    ASSERT_TRUE(biell2.dv_total > 0.0f);
    /* bielliptic_wins requires rb_limit > 5*max(r1,r2); use 6x to satisfy strict inequality */
    ASSERT_TRUE(lg_transfer_bielliptic_wins(r1, r_big, r_big * 6.0f));
    
    printf("PASS: bielliptic_transfer\n");
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
    printf("=== Lagrange Orbital Tests ===\n\n");
    
    test_t tests[] = {
        {"orbital_elements_roundtrip", test_orbital_elements_roundtrip},
        {"hohmann_transfer", test_hohmann_transfer},
        {"bielliptic_transfer", test_bielliptic_transfer},
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
