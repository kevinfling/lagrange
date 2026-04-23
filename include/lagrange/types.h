/*
 * Lagrange Physics Library - Core Types
 * Header-only C library for physics simulation
 */

#ifndef LAGRANGE_TYPES_H
#define LAGRANGE_TYPES_H

#include <stdbool.h>
#include <stdint.h>
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Entity handle (opaque) */
typedef uint64_t lg_entity_t;
#define LG_ENTITY_INVALID ((lg_entity_t)0)

/*============================================================================
 * World Configuration
 *===========================================================================*/

typedef struct {
    /* Gravity (m/s²) */
    float gravity[3];
    
    /* Time step for fixed-step simulation (seconds) */
    float time_step;
    
    /* Solver iterations */
    int velocity_iterations;
    int position_iterations;
    
    /* Performance settings */
    bool enable_sleeping;
    float sleep_threshold;      /* Linear + angular velocity magnitude */
    float linear_damping;       /* Global linear damping */
    float angular_damping;      /* Global angular damping */
    
    /* Collision settings */
    bool enable_collision;
    int collision_iterations;
    float collision_margin;     /* Default collision margin */
    
    /* World bounds (for broad-phase) */
    bool use_bounds;
    float bounds_min[3];
    float bounds_max[3];
    
    /* Memory settings */
    size_t max_entities;
    size_t max_bodies;
    size_t max_colliders;
} lg_world_config_t;

/* Default configuration */
static const lg_world_config_t LG_WORLD_CONFIG_DEFAULT = {
    .gravity = {0.0f, -9.81f, 0.0f},
    .time_step = 1.0f / 60.0f,
    .velocity_iterations = 6,
    .position_iterations = 2,
    .enable_sleeping = true,
    .sleep_threshold = 0.01f,
    .linear_damping = 0.0f,
    .angular_damping = 0.0f,
    .enable_collision = true,
    .collision_iterations = 4,
    .collision_margin = 0.01f,
    .use_bounds = false,
    .bounds_min = {-1000.0f, -1000.0f, -1000.0f},
    .bounds_max = {1000.0f, 1000.0f, 1000.0f},
    .max_entities = 10000,
    .max_bodies = 5000,
    .max_colliders = 5000
};

/*============================================================================
 * Forward Declarations
 *===========================================================================*/

typedef struct lg_world_t lg_world_t;

/*============================================================================
 * Callback Types
 *===========================================================================*/

typedef void (*lg_collision_callback_t)(
    lg_world_t* world,
    lg_entity_t entity_a,
    lg_entity_t entity_b,
    void* user_data
);

typedef void (*lg_trigger_callback_t)(
    lg_world_t* world,
    lg_entity_t trigger,
    lg_entity_t other,
    void* user_data
);

#ifdef __cplusplus
}
#endif

#endif /* LAGRANGE_TYPES_H */
