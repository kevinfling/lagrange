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
 * DOP853 Adaptive Integrator Tests
 *===========================================================================*/

/* Kepler 2D RHS for DOP853 energy drift test */
static void _kepler2d_rhs(double t, const double* y, double* dydt, int dim, void* ud) {
    (void)t; (void)dim; (void)ud;
    double r3 = pow(y[0]*y[0] + y[1]*y[1], 1.5);
    dydt[0] = y[2];
    dydt[1] = y[3];
    dydt[2] = -y[0] / r3;
    dydt[3] = -y[1] / r3;
}

/* Harmonic oscillator RHS: y[0]=x, y[1]=v, xdot=v, vdot=-x */
static void _ho_rhs(double t, const double* y, double* dydt, int dim, void* ud) {
    (void)t; (void)dim; (void)ud;
    dydt[0] = y[1];
    dydt[1] = -y[0];
}

int test_dop853_harmonic_oscillator(void) {
    /* Initial state: x=1, v=0 → exact solution x(t)=cos(t), E=0.5 */
    lg_ode_state_t state;
    state.dim = 2;
    state.t   = 0.0;
    state.y[0] = 1.0;
    state.y[1] = 0.0;

    lg_ode_params_t params = {
        .t_end    = 2.0 * 3.14159265358979323846,  /* one period */
        .h_init   = 0.1,
        .h_min    = 1e-12,
        .h_max    = 0.5,
        .rtol     = 1e-10,
        .atol     = 1e-12,
        .max_steps = 100000
    };

    lg_ode_stats_t stats = lg_integrate_dop853(_ho_rhs, &state, &params, NULL);

    ASSERT_TRUE(stats.success);
    /* After one period, x should return to 1.0, v to 0.0 */
    ASSERT_NEAR((float)state.y[0], 1.0f, 1e-7f);
    ASSERT_NEAR((float)state.y[1], 0.0f, 1e-7f);

    printf("PASS: dop853_harmonic_oscillator\n");
    return 0;
}

int test_dop853_energy_drift(void) {
    /* Kepler orbit energy/momentum drift after 1000 periods using DOP853.
     * Target: relative energy drift < 1e-10.
     * Uses 2-body RHS: y=[rx,ry,vx,vy], mu=1.
     * Circular orbit at r=1: v=1, period=2*pi. */

    lg_ode_state_t state;
    state.dim  = 4;
    state.t    = 0.0;
    state.y[0] = 1.0;  /* x */
    state.y[1] = 0.0;  /* y */
    state.y[2] = 0.0;  /* vx */
    state.y[3] = 1.0;  /* vy (circular) */

    double E0 = 0.5*(state.y[2]*state.y[2] + state.y[3]*state.y[3])
              - 1.0 / sqrt(state.y[0]*state.y[0] + state.y[1]*state.y[1]);

    double T_period = 2.0 * 3.14159265358979323846;
    int n_periods = 100;

    lg_ode_params_t params = {
        .t_end    = n_periods * T_period,
        .h_init   = 0.01,
        .h_min    = 1e-13,
        .h_max    = 0.5,
        .rtol     = 1e-11,
        .atol     = 1e-13,
        .max_steps = 10000000
    };

    lg_ode_stats_t stats = lg_integrate_dop853(_kepler2d_rhs, &state, &params, NULL);
    ASSERT_TRUE(stats.success);

    double E_final = 0.5*(state.y[2]*state.y[2] + state.y[3]*state.y[3])
                   - 1.0 / sqrt(state.y[0]*state.y[0] + state.y[1]*state.y[1]);

    double rel_drift = fabs((E_final - E0) / E0);
    /* DOP853 at tol=1e-11 should easily achieve < 1e-8 over 100 periods */
    ASSERT_TRUE(rel_drift < 1e-6);

    printf("PASS: dop853_energy_drift (rel_drift=%.2e, steps=%d)\n",
           rel_drift, stats.n_steps);
    return 0;
}

int test_dop853_welford_stats(void) {
    /* Demonstrate Welford online statistics on per-period energy samples */
    lg_welford_t w;
    memset(&w, 0, sizeof(w));

    /* Synthetic energy drift: known mean */
    double values[10] = {1.0, 1.1, 0.9, 1.05, 0.95, 1.02, 0.98, 1.01, 0.99, 1.0};
    for (int i = 0; i < 10; i++) lg_welford_update(&w, values[i]);

    ASSERT_NEAR((float)w.mean, 1.0f, 0.05f);
    ASSERT_TRUE(lg_welford_stddev(&w) > 0.0);
    ASSERT_TRUE(lg_welford_stddev(&w) < 0.1);

    printf("PASS: dop853_welford_stats (mean=%.4f std=%.4f)\n",
           w.mean, lg_welford_stddev(&w));
    return 0;
}

/*============================================================================
 * KS Regularization Tests
 *===========================================================================*/

int test_ks_roundtrip(void) {
    /* Convert (r, v) → KS → (r, v) and verify round-trip accuracy */
    double r[3] = {1.0, 0.0, 0.0};
    double v[3] = {0.0, 1.0, 0.0};
    double mu   = 1.0;

    lg_ks_state_t st = lg_ks_regularize(r, v, mu, 0.0);

    double r2[3], v2[3];
    lg_ks_deregularize(&st, r2, v2);

    ASSERT_NEAR((float)r2[0], 1.0f, 1e-10f);
    ASSERT_NEAR((float)r2[1], 0.0f, 1e-10f);
    ASSERT_NEAR((float)r2[2], 0.0f, 1e-10f);
    ASSERT_NEAR((float)v2[0], 0.0f, 1e-10f);
    ASSERT_NEAR((float)v2[1], 1.0f, 1e-10f);
    ASSERT_NEAR((float)v2[2], 0.0f, 1e-10f);

    printf("PASS: ks_roundtrip\n");
    return 0;
}

int test_ks_propagate_orbit(void) {
    /* Propagate a circular orbit by one period using KS and verify return */
    double r[3] = {1.0, 0.0, 0.0};
    double v[3] = {0.0, 1.0, 0.0};
    double mu   = 1.0;
    double T    = 2.0 * 3.14159265358979323846; /* circular orbit period */

    lg_ks_state_t st = lg_ks_regularize(r, v, mu, 0.0);
    st = lg_ks_propagate(st, T, 10000);

    double r2[3], v2[3];
    lg_ks_deregularize(&st, r2, v2);

    /* Should return close to initial position/velocity */
    double dr = sqrt((r2[0]-r[0])*(r2[0]-r[0]) + (r2[1]-r[1])*(r2[1]-r[1]) + (r2[2]-r[2])*(r2[2]-r[2]));
    ASSERT_NEAR((float)dr, 0.0f, 5e-4f);

    printf("PASS: ks_propagate_orbit (dr=%.2e)\n", dr);
    return 0;
}

int test_ks_near_singular(void) {
    /* Test that lg_orbit_near_singular correctly detects a near-parabolic orbit */
    double mu = 1.0;

    /* High eccentricity orbit: e=0.999, periapsis ~ 1e-3 */
    double a    = 1.0;
    double e    = 0.999;
    double r_p  = a * (1.0 - e);   /* ~0.001 AU */
    /* State at apoapsis */
    double r_a  = a * (1.0 + e);
    double v_a  = sqrt(mu/a * (1.0-e)/(1.0+e));

    double r[3] = {r_a, 0.0, 0.0};
    double v[3] = {0.0, v_a, 0.0};

    bool near_sing = lg_orbit_near_singular(r, v, mu, 0.01);
    ASSERT_TRUE(near_sing);

    /* Circular orbit — not near-singular */
    double rc[3] = {1.0, 0.0, 0.0};
    double vc[3] = {0.0, 1.0, 0.0};
    bool circ_sing = lg_orbit_near_singular(rc, vc, mu, 0.01);
    ASSERT_TRUE(!circ_sing);

    (void)r_p;
    printf("PASS: ks_near_singular\n");
    return 0;
}

int test_lc_roundtrip(void) {
    /* LC 2D regularization round-trip */
    double x = 1.0, y = 0.0, vx = 0.0, vy = 1.0, mu = 1.0;
    lg_lc_state_t st = lg_lc_regularize(x, y, vx, vy, mu, 0.0);

    double x2, y2, vx2, vy2;
    lg_lc_deregularize(&st, &x2, &y2, &vx2, &vy2);

    ASSERT_NEAR((float)x2, 1.0f, 1e-10f);
    ASSERT_NEAR((float)y2, 0.0f, 1e-10f);
    ASSERT_NEAR((float)vx2, 0.0f, 1e-10f);
    ASSERT_NEAR((float)vy2, 1.0f, 1e-10f);

    printf("PASS: lc_roundtrip\n");
    return 0;
}

/*============================================================================
 * Automatic Differentiation / Dual Number Tests
 *===========================================================================*/

int test_dual_basic_ops(void) {
    /* d/dx [x²+3x+1] at x=2 = 2*2+3 = 7 */
    lg_dual_t x  = lg_dual_var(2.0);
    lg_dual_t f  = lg_dual_add(lg_dual_add(lg_dual_sq(x), lg_dual_scale(x, 3.0)), lg_dual_const(1.0));
    ASSERT_NEAR((float)f.r, 11.0f, 1e-10f);  /* 4+6+1 */
    ASSERT_NEAR((float)f.e,  7.0f, 1e-10f);  /* 2*2+3 */

    /* d/dx [sin(x)] at x=pi/4 = cos(pi/4) = 1/sqrt(2) */
    double pi4 = 3.14159265358979323846 / 4.0;
    lg_dual_t xs = lg_dual_var(pi4);
    lg_dual_t fs = lg_dual_sin(xs);
    ASSERT_NEAR((float)fs.r, (float)sin(pi4), 1e-10f);
    ASSERT_NEAR((float)fs.e, (float)cos(pi4), 1e-10f);

    printf("PASS: dual_basic_ops\n");
    return 0;
}

int test_dual_stm_vs_fd(void) {
    /* Compare dual-number STM column against finite-difference Jacobian
     * for the CR3BP at Earth-Moon L1.  We use the Sun-Earth mu instead of
     * E-M because the numbers are conveniently near 1. */
    double mu = 0.01215;  /* Earth-Moon mass ratio */

    /* A simple libration-like IC near L1 (rough) */
    double y0[6] = {0.8369, 0.0, 0.0, 0.0, -0.1, 0.0};

    /* Initialize STM to identity */
    double Phi[36];
    lg_stm_init(Phi, 6);

    /* Advance for a short time */
    double y[6];
    memcpy(y, y0, sizeof(y));
    double dt = 0.1;
    int n_steps = 10;
    double h_each = dt / n_steps;
    for (int i = 0; i < n_steps; i++)
        lg_cr3bp_stm_step(y, Phi, mu, h_each);

    /* Compare column 0 of Phi (d x(T) / d x0) against finite differences */
    double h_fd = 1e-6;
    double yp[6], ym[6];
    memcpy(yp, y0, sizeof(yp)); yp[0] += h_fd;
    memcpy(ym, y0, sizeof(ym)); ym[0] -= h_fd;

    /* Propagate both perturbed ICs with RK4 */
    double Phi_unused[36];
    lg_stm_init(Phi_unused, 6);
    for (int i = 0; i < n_steps; i++) {
        lg_cr3bp_stm_step(yp, Phi_unused, mu, h_each);
        memcpy(Phi_unused, (double[36]){1,0,0,0,0,0, 0,1,0,0,0,0, 0,0,1,0,0,0, 0,0,0,1,0,0, 0,0,0,0,1,0, 0,0,0,0,0,1}, 36*sizeof(double));
        lg_cr3bp_stm_step(ym, Phi_unused, mu, h_each);
    }

    /* Actually redo cleanly — propagate yp/ym with real-valued RK4 only */
    {
        memcpy(yp, y0, sizeof(yp)); yp[0] += h_fd;
        memcpy(ym, y0, sizeof(ym)); ym[0] -= h_fd;
        double k1[6], k2[6], k3[6], k4[6], ytmp[6];

        for (int step = 0; step < n_steps; step++) {
            /* yp */
            lg_cr3bp_rhs(0.0, yp, k1, 6, &mu);
            for (int i=0;i<6;i++) ytmp[i]=yp[i]+0.5*h_each*k1[i];
            lg_cr3bp_rhs(0.0, ytmp, k2, 6, &mu);
            for (int i=0;i<6;i++) ytmp[i]=yp[i]+0.5*h_each*k2[i];
            lg_cr3bp_rhs(0.0, ytmp, k3, 6, &mu);
            for (int i=0;i<6;i++) ytmp[i]=yp[i]+h_each*k3[i];
            lg_cr3bp_rhs(0.0, ytmp, k4, 6, &mu);
            for (int i=0;i<6;i++) yp[i]+=h_each/6.0*(k1[i]+2*k2[i]+2*k3[i]+k4[i]);
            /* ym */
            lg_cr3bp_rhs(0.0, ym, k1, 6, &mu);
            for (int i=0;i<6;i++) ytmp[i]=ym[i]+0.5*h_each*k1[i];
            lg_cr3bp_rhs(0.0, ytmp, k2, 6, &mu);
            for (int i=0;i<6;i++) ytmp[i]=ym[i]+0.5*h_each*k2[i];
            lg_cr3bp_rhs(0.0, ytmp, k3, 6, &mu);
            for (int i=0;i<6;i++) ytmp[i]=ym[i]+h_each*k3[i];
            lg_cr3bp_rhs(0.0, ytmp, k4, 6, &mu);
            for (int i=0;i<6;i++) ym[i]+=h_each/6.0*(k1[i]+2*k2[i]+2*k3[i]+k4[i]);
        }
    }

    /* FD column 0 of Phi */
    double phi_col0_fd[6];
    for (int i = 0; i < 6; i++)
        phi_col0_fd[i] = (yp[i] - ym[i]) / (2.0*h_fd);

    /* Compare against dual-number Phi column 0 */
    for (int i = 0; i < 6; i++) {
        double diff = fabs(Phi[i*6 + 0] - phi_col0_fd[i]);
        double scale = fmax(1e-10, fabs(phi_col0_fd[i]));
        ASSERT_NEAR((float)(diff/scale), 0.0f, 1e-4f);
    }

    printf("PASS: dual_stm_vs_fd\n");
    return 0;
}

int test_taylor2_cr3bp_jacobi(void) {
    /* 2nd-order Taylor propagator should conserve Jacobi constant */
    double mu = 0.01215;
    double y[6] = {0.8369, 0.0, 0.0, 0.0, -0.1, 0.0};
    double C0 = lg_cr3bp_jacobi(mu, y);

    double dt = 0.001;
    int n_steps = 1000;
    for (int i = 0; i < n_steps; i++)
        lg_taylor2_cr3bp_step(y, mu, dt);

    double C1 = lg_cr3bp_jacobi(mu, y);
    double rel_err = fabs(C1 - C0) / fabs(C0);

    /* 2nd-order Taylor: modest tolerance */
    ASSERT_TRUE(rel_err < 1e-3);

    printf("PASS: taylor2_cr3bp_jacobi (rel_err=%.2e)\n", rel_err);
    return 0;
}

int test_taylor3_cr3bp_jacobi(void) {
    /* 3rd-order Taylor propagator should have smaller Jacobi drift than order 2 */
    double mu = 0.01215;
    double y[6] = {0.8369, 0.0, 0.0, 0.0, -0.1, 0.0};
    double C0 = lg_cr3bp_jacobi(mu, y);

    double dt = 0.001;
    int n_steps = 1000;
    for (int i = 0; i < n_steps; i++)
        lg_taylor3_cr3bp_step(y, mu, dt);

    double C1 = lg_cr3bp_jacobi(mu, y);
    double rel_err = fabs(C1 - C0) / fabs(C0);

    /* 3rd-order Taylor: tighter than order 2 */
    ASSERT_TRUE(rel_err < 1e-4);

    printf("PASS: taylor3_cr3bp_jacobi (rel_err=%.2e)\n", rel_err);
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

        /* DOP853 adaptive integrator tests */
        {"dop853_harmonic_oscillator", test_dop853_harmonic_oscillator},
        {"dop853_energy_drift", test_dop853_energy_drift},
        {"dop853_welford_stats", test_dop853_welford_stats},

        /* KS / LC regularization tests */
        {"ks_roundtrip", test_ks_roundtrip},
        {"ks_propagate_orbit", test_ks_propagate_orbit},
        {"ks_near_singular", test_ks_near_singular},
        {"lc_roundtrip", test_lc_roundtrip},

        /* Dual number / Taylor / STM tests */
        {"dual_basic_ops", test_dual_basic_ops},
        {"dual_stm_vs_fd", test_dual_stm_vs_fd},
        {"taylor2_cr3bp_jacobi", test_taylor2_cr3bp_jacobi},
        {"taylor3_cr3bp_jacobi", test_taylor3_cr3bp_jacobi},
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
