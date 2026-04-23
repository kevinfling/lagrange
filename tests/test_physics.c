/*
 * Lagrange Tests - Physics Module
 */

#include <stdio.h>
#include <math.h>
#include <stdbool.h>

#define LAGRANGE_IMPLEMENTATION
#include "lagrange.h"

#define EPSILON 1e-4f

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
 * Body Tests
 *===========================================================================*/

int test_body_mass(void) {
    lg_body_t b = lg_body(10.0f);
    ASSERT_NEAR(b.mass, 10.0f, EPSILON);
    ASSERT_NEAR(b.inv_mass, 0.1f, EPSILON);
    
    lg_body_set_mass(&b, 5.0f);
    ASSERT_NEAR(b.mass, 5.0f, EPSILON);
    ASSERT_NEAR(b.inv_mass, 0.2f, EPSILON);
    
    /* Static body */
    lg_body_t s = lg_body_static();
    ASSERT_NEAR(s.mass, 0.0f, EPSILON);
    ASSERT_NEAR(s.inv_mass, 0.0f, EPSILON);
    
    printf("PASS: body_mass\n");
    return 0;
}

int test_body_forces(void) {
    lg_body_t b = lg_body(1.0f);
    
    lg_vec3_t force = lg_vec3(10.0f, 0.0f, 0.0f);
    lg_body_apply_force(&b, force);
    
    ASSERT_VEC3_NEAR(b.force, 10.0f, 0.0f, 0.0f, EPSILON);
    
    /* Apply another force */
    lg_body_apply_force(&b, lg_vec3(5.0f, 5.0f, 0.0f));
    ASSERT_VEC3_NEAR(b.force, 15.0f, 5.0f, 0.0f, EPSILON);
    
    /* Clear forces */
    lg_body_clear_forces(&b);
    ASSERT_VEC3_NEAR(b.force, 0.0f, 0.0f, 0.0f, EPSILON);
    
    printf("PASS: body_forces\n");
    return 0;
}

int test_body_impulse(void) {
    lg_body_t b = lg_body(2.0f);
    
    /* Apply impulse */
    lg_body_apply_impulse(&b, lg_vec3(10.0f, 0.0f, 0.0f));
    
    /* v = impulse / mass = 10 / 2 = 5 */
    ASSERT_VEC3_NEAR(b.velocity, 5.0f, 0.0f, 0.0f, EPSILON);
    
    printf("PASS: body_impulse\n");
    return 0;
}

int test_body_kinetic_energy(void) {
    /* KE = 0.5 * m * v² */
    lg_body_t b = lg_body(2.0f);
    b.velocity = lg_vec3(3.0f, 0.0f, 0.0f);
    
    float ke = lg_body_kinetic_energy(&b);
    /* KE = 0.5 * 2 * 9 = 9 */
    ASSERT_NEAR(ke, 9.0f, EPSILON);
    
    printf("PASS: body_kinetic_energy\n");
    return 0;
}

int test_body_sleep(void) {
    lg_body_t b = lg_body(1.0f);
    ASSERT_TRUE(!b.is_sleeping);
    
    b.allow_sleep = true;
    b.velocity = lg_vec3(0.001f, 0.0f, 0.0f);
    
    bool should_sleep = lg_body_should_sleep(&b, 0.01f);
    ASSERT_TRUE(should_sleep);
    
    lg_body_sleep(&b);
    ASSERT_TRUE(b.is_sleeping);
    
    lg_body_wake(&b);
    ASSERT_TRUE(!b.is_sleeping);
    
    printf("PASS: body_sleep\n");
    return 0;
}

/*============================================================================
 * Transform Tests
 *===========================================================================*/

int test_transform_identity(void) {
    lg_transform_t t = lg_transform();
    
    ASSERT_VEC3_NEAR(t.position, 0.0f, 0.0f, 0.0f, EPSILON);
    ASSERT_NEAR(t.rotation.w, 1.0f, EPSILON);
    ASSERT_VEC3_NEAR(t.scale, 1.0f, 1.0f, 1.0f, EPSILON);
    
    printf("PASS: transform_identity\n");
    return 0;
}

int test_transform_point(void) {
    lg_transform_t t = lg_transform_at(lg_vec3(10.0f, 0.0f, 0.0f));
    
    lg_vec3_t local = lg_vec3(1.0f, 0.0f, 0.0f);
    lg_vec3_t world = lg_transform_point(t, local);
    
    ASSERT_VEC3_NEAR(world, 11.0f, 0.0f, 0.0f, EPSILON);
    
    printf("PASS: transform_point\n");
    return 0;
}

/*============================================================================
 * Collider Tests
 *===========================================================================*/

int test_collider_volume(void) {
    /* Sphere volume = 4/3 * pi * r³ */
    lg_collider_t sphere = lg_collider_sphere(1.0f);
    float vol = lg_collider_volume(&sphere);
    ASSERT_NEAR(vol, 4.18879f, 0.001f); /* 4/3 * pi */
    
    /* Box volume = 8 * hx * hy * hz */
    lg_collider_t box = lg_collider_box(1.0f, 2.0f, 3.0f);
    vol = lg_collider_volume(&box);
    ASSERT_NEAR(vol, 48.0f, EPSILON); /* 8 * 1 * 2 * 3 */
    
    printf("PASS: collider_volume\n");
    return 0;
}

int test_collider_inertia(void) {
    /* Sphere inertia = 0.4 * m * r² */
    lg_collider_t sphere = lg_collider_sphere(1.0f);
    lg_vec3_t i = lg_collider_inertia(&sphere, 10.0f);
    ASSERT_NEAR(i.x, 4.0f, EPSILON); /* 0.4 * 10 * 1 */
    ASSERT_NEAR(i.y, 4.0f, EPSILON);
    ASSERT_NEAR(i.z, 4.0f, EPSILON);
    
    printf("PASS: collider_inertia\n");
    return 0;
}

/*============================================================================
 * World Tests
 *===========================================================================*/

int test_world_create_destroy(void) {
    lg_world_t* world = lg_world_create(NULL);
    ASSERT_TRUE(world != NULL);
    
    ASSERT_NEAR(world->config.gravity[1], -9.81f, EPSILON);
    
    lg_world_destroy(world);
    
    printf("PASS: world_create_destroy\n");
    return 0;
}

int test_entity_create_destroy(void) {
    lg_world_t* world = lg_world_create(NULL);
    
    lg_entity_t e1 = lg_entity_create(world);
    ASSERT_TRUE(e1 != LG_ENTITY_INVALID);
    
    lg_entity_t e2 = lg_entity_create(world);
    ASSERT_TRUE(e2 != LG_ENTITY_INVALID);
    ASSERT_TRUE(e1 != e2);
    
    ASSERT_TRUE(lg_entity_valid(world, e1));
    ASSERT_TRUE(lg_entity_valid(world, e2));
    
    lg_entity_destroy(world, e1);
    ASSERT_TRUE(!lg_entity_valid(world, e1));
    ASSERT_TRUE(lg_entity_valid(world, e2));
    
    lg_world_destroy(world);
    
    printf("PASS: entity_create_destroy\n");
    return 0;
}

int test_world_components(void) {
    lg_world_t* world = lg_world_create(NULL);
    
    lg_entity_t e = lg_entity_create(world);
    
    /* Set position */
    lg_set_position(world, e, lg_vec3(1.0f, 2.0f, 3.0f));
    
    lg_vec3_t pos = lg_get_position(world, e);
    ASSERT_VEC3_NEAR(pos, 1.0f, 2.0f, 3.0f, EPSILON);
    
    /* Set velocity */
    lg_set_velocity(world, e, lg_vec3(10.0f, 0.0f, 0.0f));
    
    lg_vec3_t vel = lg_get_velocity(world, e);
    ASSERT_VEC3_NEAR(vel, 10.0f, 0.0f, 0.0f, EPSILON);
    
    /* Set mass */
    lg_set_mass(world, e, 5.0f);
    ASSERT_NEAR(lg_get_mass(world, e), 5.0f, EPSILON);
    
    lg_world_destroy(world);
    
    printf("PASS: world_components\n");
    return 0;
}

int test_entity_recycling(void) {
    lg_world_t* world = lg_world_create(NULL);
    
    lg_entity_t e1 = lg_entity_create(world);
    lg_entity_t e2 = lg_entity_create(world);
    ASSERT_TRUE(e1 != LG_ENTITY_INVALID);
    ASSERT_TRUE(e2 != LG_ENTITY_INVALID);
    
    lg_entity_destroy(world, e1);
    
    /* Create new entity - should reuse e1's ID */
    lg_entity_t e3 = lg_entity_create(world);
    ASSERT_TRUE(e3 == e1);
    ASSERT_TRUE(!lg_entity_valid(world, e2) || e3 != e2); /* e3 should be e1, not e2 */
    ASSERT_TRUE(lg_entity_valid(world, e3));
    
    lg_world_destroy(world);
    
    printf("PASS: entity_recycling\n");
    return 0;
}

static int g_collision_count = 0;
static lg_entity_t g_collision_a = LG_ENTITY_INVALID;
static lg_entity_t g_collision_b = LG_ENTITY_INVALID;

static void test_collision_handler(lg_world_t* world, lg_entity_t a, lg_entity_t b, void* user_data) {
    (void)world;
    (void)user_data;
    g_collision_count++;
    g_collision_a = a;
    g_collision_b = b;
}

int test_collision_callback(void) {
    lg_world_t* world = lg_world_create(NULL);
    world->config.collision_iterations = 1;
    
    lg_world_set_collision_callback(world, test_collision_handler, NULL);
    
    /* Create two spheres that overlap */
    lg_entity_t e1 = lg_entity_create(world);
    lg_entity_t e2 = lg_entity_create(world);
    
    lg_set_position(world, e1, lg_vec3(0.0f, 0.0f, 0.0f));
    lg_set_position(world, e2, lg_vec3(0.5f, 0.0f, 0.0f));
    
    lg_collider_t c1 = lg_collider_sphere(0.5f);
    lg_collider_t c2 = lg_collider_sphere(0.5f);
    lg_set_collider(world, e1, &c1);
    lg_set_collider(world, e2, &c2);
    
    g_collision_count = 0;
    g_collision_a = LG_ENTITY_INVALID;
    g_collision_b = LG_ENTITY_INVALID;
    
    lg_world_step(world);
    
    ASSERT_TRUE(g_collision_count >= 1);
    ASSERT_TRUE((g_collision_a == e1 && g_collision_b == e2) ||
                (g_collision_a == e2 && g_collision_b == e1));
    
    lg_world_destroy(world);
    
    printf("PASS: collision_callback\n");
    return 0;
}

int test_transform_interpolation(void) {
    lg_world_t* world = lg_world_create(NULL);
    world->config.time_step = 1.0f; /* Large step for clear differences */
    world->use_gravity = false;
    
    lg_entity_t e = lg_entity_create(world);
    lg_set_position(world, e, lg_vec3(0.0f, 0.0f, 0.0f));
    
    lg_body_t body = lg_body(1.0f);
    body.velocity = lg_vec3(10.0f, 0.0f, 0.0f);
    body.linear_damping = 0.0f; /* Disable damping for exact test */
    lg_set_body(world, e, &body);
    
    /* Step once: position should move to (10, 0, 0) */
    lg_world_step(world);
    
    lg_vec3_t current = lg_get_position(world, e);
    ASSERT_VEC3_NEAR(current, 10.0f, 0.0f, 0.0f, EPSILON);
    
    /* Set accumulator to 0.5 for alpha = 0.5 */
    world->accumulator = 0.5f;
    
    lg_transform_t interp = lg_get_interpolated_transform(world, e);
    ASSERT_VEC3_NEAR(interp.position, 5.0f, 0.0f, 0.0f, EPSILON);
    
    lg_world_destroy(world);
    
    printf("PASS: transform_interpolation\n");
    return 0;
}

/*============================================================================
 * Integration Tests
 *===========================================================================*/

int test_integration_euler(void) {
    lg_world_t* world = lg_world_create(NULL);
    world->integrator = LG_INTEGRATOR_SEMI_IMPLICIT_EULER;
    
    lg_entity_t e = lg_entity_create(world);
    lg_set_position(world, e, lg_vec3(0.0f, 10.0f, 0.0f));
    
    /* Set up a dynamic body */
    lg_body_t body = lg_body(1.0f);
    lg_set_body(world, e, &body);
    
    /* Verify body is dynamic */
    lg_body_t* b = lg_get_body(world, e);
    ASSERT_TRUE(b != NULL);
    ASSERT_TRUE(b->type == LG_BODY_DYNAMIC);
    
    /* Simulate falling under gravity for 1 second */
    for (int i = 0; i < 60; i++) {
        lg_world_step(world);
    }
    
    lg_vec3_t pos = lg_get_position(world, e);
    
    /* Should have fallen significantly (rough check)
     * With g=-9.81, after 1s: y = 10 - 0.5*9.81*1² ≈ 5.1
     */
    ASSERT_TRUE(pos.y < 6.0f);
    
    lg_world_destroy(world);
    
    printf("PASS: integration_euler\n");
    return 0;
}

/*============================================================================
 * Gravity Tests
 *===========================================================================*/

int test_orbital_velocity(void) {
    /* Earth orbit at 7000 km from center */
    double r = 7e6;
    double v = lg_orbital_velocity_circular(LG_MU_EARTH, r);
    
    /* Should be around 7.5 km/s */
    ASSERT_NEAR((float)v, 7500.0f, 100.0f);
    
    /* Verify: v² = mu / r */
    double v_sq_calc = LG_MU_EARTH / r;
    ASSERT_NEAR((float)(v * v), (float)v_sq_calc, 1000.0f);
    
    printf("PASS: orbital_velocity\n");
    return 0;
}

int test_gravity_accel(void) {
    /* Surface gravity of Earth */
    lg_vec3_t surface = lg_vec3(LG_RADIUS_EARTH, 0.0f, 0.0f);
    lg_vec3_t center = lg_vec3_zero();
    
    lg_vec3_t accel = lg_gravity_accel(surface, center, LG_MU_EARTH);
    
    /* Should be approximately 9.8 m/s² towards center */
    float g = lg_vec3_len(accel);
    ASSERT_NEAR(g, 9.8f, 0.1f);
    
    printf("PASS: gravity_accel\n");
    return 0;
}

/*============================================================================
 * Integrator Tests (Beyond Euler)
 *===========================================================================*/

int test_integrator_verlet(void) {
    lg_world_t* world = lg_world_create(NULL);
    world->config.time_step = 1.0f / 60.0f;
    world->use_gravity = false;
    world->integrator = LG_INTEGRATOR_VELOCITY_VERLET;
    
    lg_entity_t e = lg_entity_create(world);
    lg_set_position(world, e, lg_vec3(1.0f, 0.0f, 0.0f));
    
    lg_body_t body = lg_body(1.0f);
    body.linear_damping = 0.0f;
    lg_set_body(world, e, &body);
    
    /* Spring-like restoring force F = -kx */
    lg_body_t* b = lg_get_body(world, e);
    b->force = lg_vec3(-1.0f, 0.0f, 0.0f);
    
    /* Reset verlet state for clean first step */
    lg_verlet_state_t* state = &world->storage.verlet_states[lg_storage_find(&world->storage, e)];
    state->first_step = true;
    
    float total_energy = 0.0f;
    float initial_energy = 0.0f;
    for (int i = 0; i < 120; i++) {
        (void)b; /* b is not needed here */
        /* Actually we need to apply force each step via body */
        lg_body_apply_force(b, lg_vec3(-1.0f * lg_get_position(world, e).x, 0.0f, 0.0f));
        lg_world_step(world);
        lg_body_clear_forces(b);
        
        float ke = 0.5f * b->mass * lg_vec3_len_sq(b->velocity);
        float pe = 0.5f * lg_get_position(world, e).x * lg_get_position(world, e).x;
        total_energy = ke + pe;
        if (i == 0) initial_energy = total_energy;
    }
    
    /* Verlet should conserve energy well */
    ASSERT_NEAR(total_energy, initial_energy, 0.05f * initial_energy);
    
    lg_world_destroy(world);
    printf("PASS: integrator_verlet\n");
    return 0;
}

int test_integrator_rk4(void) {
    lg_body_t body = lg_body(1.0f);
    body.linear_damping = 0.0f;
    lg_transform_t xform = lg_transform();
    xform.position = lg_vec3(1.0f, 0.0f, 0.0f);
    
    float dt = 1.0f / 60.0f;
    float initial_energy = 0.5f; /* PE = 0.5 * 1^2 */
    
    for (int i = 0; i < 120; i++) {
        body.force = lg_vec3(-xform.position.x, 0.0f, 0.0f);
        lg_integrate_rk4(&body, &xform, dt);
    }
    
    float ke = 0.5f * lg_vec3_len_sq(body.velocity);
    float pe = 0.5f * xform.position.x * xform.position.x;
    float total_energy = ke + pe;
    
    /* RK4 should conserve energy very well */
    ASSERT_NEAR(total_energy, initial_energy, 0.01f * initial_energy);
    
    printf("PASS: integrator_rk4\n");
    return 0;
}

/*============================================================================
 * Sleeping System
 *===========================================================================*/

int test_sleeping_system(void) {
    lg_world_t* world = lg_world_create(NULL);
    world->config.sleep_threshold = 0.001f;
    world->config.time_step = 0.1f;
    world->use_gravity = false;
    
    lg_entity_t e = lg_entity_create(world);
    lg_body_t body = lg_body(1.0f);
    body.linear_damping = 0.0f;
    body.allow_sleep = true;
    lg_set_body(world, e, &body);
    lg_set_position(world, e, lg_vec3(0.0f, 0.0f, 0.0f));
    
    lg_body_t* b = lg_get_body(world, e);
    b->velocity = lg_vec3(0.01f, 0.0f, 0.0f);
    
    /* Simulate until body should sleep */
    bool slept = false;
    for (int i = 0; i < 20; i++) {
        if (lg_integrator_check_sleep(b, 0.001f, 0.1f)) {
            b->is_sleeping = true;
            slept = true;
        }
        lg_world_step(world);
    }
    
    ASSERT_TRUE(slept);
    ASSERT_TRUE(b->is_sleeping);
    
    /* Wake on impulse */
    lg_body_wake(b);
    ASSERT_TRUE(!b->is_sleeping);
    
    lg_world_destroy(world);
    printf("PASS: sleeping_system\n");
    return 0;
}

/*============================================================================
 * Variable Timestep
 *===========================================================================*/

int test_variable_timestep(void) {
    lg_world_t* world = lg_world_create(NULL);
    world->config.time_step = 0.1f;
    world->use_gravity = false;
    
    lg_entity_t e = lg_entity_create(world);
    lg_set_position(world, e, lg_vec3(0.0f, 0.0f, 0.0f));
    lg_body_t body = lg_body(1.0f);
    body.linear_damping = 0.0f;
    lg_set_body(world, e, &body);
    lg_set_velocity(world, e, lg_vec3(1.0f, 0.0f, 0.0f));
    
    /* Call update with dt = 0.25s -> should produce 2 steps (0.2s) and accumulator = 0.05s */
    /* (Use 0.25 because lg_world_update clamps accumulator to 0.25 max) */
    lg_world_update(world, 0.25f);
    
    /* Should have stepped 2 times */
    ASSERT_TRUE(world->step_count == 2);
    ASSERT_NEAR(world->accumulator, 0.05f, 1e-4f);
    
    /* Position should be v * 2 * dt = 1 * 0.2 = 0.2 */
    lg_vec3_t pos = lg_get_position(world, e);
    ASSERT_NEAR(pos.x, 0.2f, 1e-3f);
    
    lg_world_destroy(world);
    printf("PASS: variable_timestep\n");
    return 0;
}

/*============================================================================
 * Momentum Conservation in Collisions
 *===========================================================================*/

int test_momentum_conservation_collision(void) {
    lg_world_t* world = lg_world_create(NULL);
    world->config.time_step = 1.0f / 60.0f;
    world->config.collision_iterations = 1;
    world->use_gravity = false;
    
    /* Two equal masses, elastic head-on collision */
    lg_entity_t e1 = lg_entity_create(world);
    lg_entity_t e2 = lg_entity_create(world);
    
    lg_set_position(world, e1, lg_vec3(-1.0f, 0.0f, 0.0f));
    lg_set_position(world, e2, lg_vec3(1.0f, 0.0f, 0.0f));
    
    lg_body_t b1 = lg_body(1.0f);
    b1.velocity = lg_vec3(5.0f, 0.0f, 0.0f);
    b1.linear_damping = 0.0f;
    lg_set_body(world, e1, &b1);
    
    lg_body_t b2 = lg_body(1.0f);
    b2.velocity = lg_vec3(-5.0f, 0.0f, 0.0f);
    b2.linear_damping = 0.0f;
    lg_set_body(world, e2, &b2);
    
    lg_collider_t c1 = lg_collider_sphere(1.5f);
    lg_collider_t c2 = lg_collider_sphere(1.5f);
    lg_material_t m = LG_MATERIAL_DEFAULT;
    m.restitution = 1.0f; /* Elastic */
    m.friction = 0.0f;
    
    lg_set_collider(world, e1, &c1);
    lg_set_collider(world, e2, &c2);
    lg_set_material(world, e1, &m);
    lg_set_material(world, e2, &m);
    
    float initial_momentum = 5.0f + (-5.0f); /* Should be 0 */
    
    /* Step until they collide and separate */
    for (int i = 0; i < 10; i++) {
        lg_world_step(world);
    }
    
    lg_vec3_t v1 = lg_get_velocity(world, e1);
    lg_vec3_t v2 = lg_get_velocity(world, e2);
    float final_momentum = v1.x + v2.x;
    
    /* Momentum should be conserved */
    ASSERT_NEAR(final_momentum, initial_momentum, 0.5f);
    
    lg_world_destroy(world);
    printf("PASS: momentum_conservation_collision\n");
    return 0;
}

/*============================================================================
 * Stack Stability
 *===========================================================================*/

int test_stack_stability(void) {
    lg_world_t* world = lg_world_create(NULL);
    world->config.time_step = 1.0f / 60.0f;
    world->config.gravity[1] = -9.81f;
    world->use_gravity = true;
    
    /* Create a stack of 5 boxes */
    lg_collider_t box = lg_collider_box(0.5f, 0.5f, 0.5f);
    lg_material_t mat = LG_MATERIAL_DEFAULT;
    mat.restitution = 0.0f; /* Inelastic to help settling */
    mat.friction = 0.5f;
    
    for (int i = 0; i < 5; i++) {
        lg_entity_t e = lg_entity_create(world);
        lg_set_position(world, e, lg_vec3(0.0f, 0.5f + i * 1.0f, 0.0f));
        lg_body_t body = lg_body(1.0f);
        lg_set_body(world, e, &body);
        lg_set_collider(world, e, &box);
        lg_set_material(world, e, &mat);
    }
    
    /* Add a ground plane */
    lg_entity_t ground = lg_entity_create(world);
    lg_set_position(world, ground, lg_vec3(0.0f, 0.0f, 0.0f));
    lg_collider_t plane = lg_collider_plane(lg_vec3(0.0f, 1.0f, 0.0f), 0.0f);
    lg_body_t ground_body = lg_body_static();
    lg_set_body(world, ground, &ground_body);
    lg_set_collider(world, ground, &plane);
    
    /* Simulate for 2 seconds */
    for (int i = 0; i < 120; i++) {
        lg_world_step(world);
    }
    
    /* Stack should settle without blowing up (all y positions should be reasonable) */
    for (size_t i = 1; i < world->storage.count; i++) { /* skip ground */
        lg_vec3_t pos = world->storage.transforms[i].position;
        ASSERT_TRUE(pos.y >= 0.0f);
        ASSERT_TRUE(pos.y < 10.0f); /* If stack explodes, y will be huge */
    }
    
    lg_world_destroy(world);
    printf("PASS: stack_stability\n");
    return 0;
}

/*============================================================================
 * Collision Detection Tests
 *===========================================================================*/

int test_collision_sphere_box(void) {
    lg_contact_t contact;
    bool hit = lg_collide_sphere_box(
        lg_vec3(1.0f, 0.0f, 0.0f), 0.6f,
        lg_vec3(0.0f, 0.0f, 0.0f), lg_vec3(1.0f, 1.0f, 1.0f),
        &contact
    );
    ASSERT_TRUE(hit);
    ASSERT_NEAR(contact.penetration, 0.6f, EPSILON);
    ASSERT_NEAR(contact.normal.x, 1.0f, EPSILON);

    /* No collision */
    hit = lg_collide_sphere_box(
        lg_vec3(2.0f, 0.0f, 0.0f), 0.5f,
        lg_vec3(0.0f, 0.0f, 0.0f), lg_vec3(1.0f, 1.0f, 1.0f),
        &contact
    );
    ASSERT_TRUE(!hit);

    printf("PASS: collision_sphere_box\n");
    return 0;
}

int test_collision_box_box(void) {
    lg_contact_t contact;
    bool hit = lg_collide_box_box(
        lg_vec3(0.0f, 0.0f, 0.0f), lg_vec3(1.0f, 1.0f, 1.0f),
        lg_vec3(1.5f, 0.0f, 0.0f), lg_vec3(1.0f, 1.0f, 1.0f),
        &contact
    );
    ASSERT_TRUE(hit);
    ASSERT_NEAR(contact.penetration, 0.5f, EPSILON);
    ASSERT_NEAR(contact.normal.x, -1.0f, EPSILON);

    /* No collision */
    hit = lg_collide_box_box(
        lg_vec3(0.0f, 0.0f, 0.0f), lg_vec3(1.0f, 1.0f, 1.0f),
        lg_vec3(3.0f, 0.0f, 0.0f), lg_vec3(1.0f, 1.0f, 1.0f),
        &contact
    );
    ASSERT_TRUE(!hit);

    printf("PASS: collision_box_box\n");
    return 0;
}

int test_collision_sphere_capsule(void) {
    lg_contact_t contact;
    bool hit = lg_collide_sphere_capsule(
        lg_vec3(0.0f, 1.8f, 0.0f), 0.5f,
        lg_vec3(0.0f, 0.0f, 0.0f), lg_vec3(0.0f, 1.0f, 0.0f), 1.0f, 0.5f,
        &contact
    );
    ASSERT_TRUE(hit);
    ASSERT_NEAR(contact.normal.y, 1.0f, EPSILON);
    ASSERT_NEAR(contact.penetration, 0.2f, EPSILON);

    /* No collision */
    hit = lg_collide_sphere_capsule(
        lg_vec3(0.0f, 3.0f, 0.0f), 0.5f,
        lg_vec3(0.0f, 0.0f, 0.0f), lg_vec3(0.0f, 1.0f, 0.0f), 1.0f, 0.5f,
        &contact
    );
    ASSERT_TRUE(!hit);

    printf("PASS: collision_sphere_capsule\n");
    return 0;
}

int test_collision_capsule_capsule(void) {
    lg_contact_t contact;
    bool hit = lg_collide_capsule_capsule(
        lg_vec3(0.0f, 0.0f, 0.0f), lg_vec3(0.0f, 1.0f, 0.0f), 1.0f, 0.5f,
        lg_vec3(0.0f, 2.5f, 0.0f), lg_vec3(0.0f, 1.0f, 0.0f), 1.0f, 0.5f,
        &contact
    );
    ASSERT_TRUE(hit);
    ASSERT_NEAR(contact.normal.y, -1.0f, EPSILON);
    ASSERT_NEAR(contact.penetration, 0.5f, EPSILON);

    /* No collision */
    hit = lg_collide_capsule_capsule(
        lg_vec3(0.0f, 0.0f, 0.0f), lg_vec3(0.0f, 1.0f, 0.0f), 1.0f, 0.5f,
        lg_vec3(0.0f, 4.0f, 0.0f), lg_vec3(0.0f, 1.0f, 0.0f), 1.0f, 0.5f,
        &contact
    );
    ASSERT_TRUE(!hit);

    printf("PASS: collision_capsule_capsule\n");
    return 0;
}

/*============================================================================
 * Attitude Control Tests
 *===========================================================================*/

int test_attitude_control_rw_torque(void) {
    lg_rw_array_t arr = lg_rw_array_triad();
    lg_vec3_t tau_cmd = lg_vec3(1.0f, 0.0f, 0.0f);
    lg_rw_array_command_torque(&arr, &tau_cmd);

    ASSERT_NEAR(arr.wheels[0].tau_command, 1.0f, EPSILON);
    ASSERT_NEAR(arr.wheels[1].tau_command, 0.0f, EPSILON);
    ASSERT_NEAR(arr.wheels[2].tau_command, 0.0f, EPSILON);
    ASSERT_TRUE(!arr.singularity_warn);

    printf("PASS: attitude_control_rw_torque\n");
    return 0;
}

int test_attitude_control_cmg_steering(void) {
    lg_cmg_array_t arr = lg_cmg_array_iss_like();
    lg_vec3_t tau_cmd = lg_vec3(0.1f, 0.0f, 0.0f);
    lg_cmg_steering_pinv(&arr, &tau_cmd, 0.01f);

    for (int i = 0; i < arr.n_cmgs; i++) {
        ASSERT_TRUE(!isnan(arr.cmgs[i].gimbal_rate_cmd));
        ASSERT_TRUE(!isinf(arr.cmgs[i].gimbal_rate_cmd));
    }

    printf("PASS: attitude_control_cmg_steering\n");
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
    printf("=== Lagrange Physics Tests ===\n\n");
    
    test_t tests[] = {
        /* Body tests */
        {"body_mass", test_body_mass},
        {"body_forces", test_body_forces},
        {"body_impulse", test_body_impulse},
        {"body_kinetic_energy", test_body_kinetic_energy},
        {"body_sleep", test_body_sleep},
        
        /* Transform tests */
        {"transform_identity", test_transform_identity},
        {"transform_point", test_transform_point},
        
        /* Collider tests */
        {"collider_volume", test_collider_volume},
        {"collider_inertia", test_collider_inertia},
        
        /* World tests */
        {"world_create_destroy", test_world_create_destroy},
        {"entity_create_destroy", test_entity_create_destroy},
        {"world_components", test_world_components},
        {"entity_recycling", test_entity_recycling},
        {"collision_callback", test_collision_callback},
        {"transform_interpolation", test_transform_interpolation},
        
        /* Integration tests */
        {"integration_euler", test_integration_euler},
        {"integrator_verlet", test_integrator_verlet},
        {"integrator_rk4", test_integrator_rk4},
        {"sleeping_system", test_sleeping_system},
        {"variable_timestep", test_variable_timestep},
        {"momentum_conservation_collision", test_momentum_conservation_collision},
        {"stack_stability", test_stack_stability},

        /* Collision detection tests */
        {"collision_sphere_box", test_collision_sphere_box},
        {"collision_box_box", test_collision_box_box},
        {"collision_sphere_capsule", test_collision_sphere_capsule},
        {"collision_capsule_capsule", test_collision_capsule_capsule},

        /* Attitude control tests */
        {"attitude_control_rw_torque", test_attitude_control_rw_torque},
        {"attitude_control_cmg_steering", test_attitude_control_cmg_steering},

        /* Gravity tests */
        {"orbital_velocity", test_orbital_velocity},
        {"gravity_accel", test_gravity_accel},
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
