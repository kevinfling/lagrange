/*
 * Lagrange Physics Library - Rigid Body Component
 * Mass, velocity, forces, and dynamics
 */

#ifndef LAGRANGE_BODY_H
#define LAGRANGE_BODY_H

#include "types.h"
#include "math.h"

#ifdef __cplusplus
extern "C" {
#endif

/*============================================================================
 * Body Type
 *===========================================================================*/

typedef enum {
    LG_BODY_DYNAMIC,    /* Affected by forces, can move */
    LG_BODY_KINEMATIC,  /* Not affected by forces, moved by code */
    LG_BODY_STATIC      /* Immovable, infinite mass */
} lg_body_type_t;

/*============================================================================
 * Rigid Body Component
 *===========================================================================*/

typedef struct {
    /* Mass properties */
    float mass;
    float inv_mass;         /* 1/mass, 0 for static */
    float radius;           /* For collision, visualization, SOI calc */
    
    /* Inertia tensor (diagonal for now) */
    lg_vec3_t inertia;      /* Principal moments */
    lg_vec3_t inv_inertia;  /* 1/inertia */
    
    /* Linear motion */
    lg_vec3_t velocity;
    lg_vec3_t force;
    
    /* Angular motion */
    lg_vec3_t angular_velocity;
    lg_vec3_t torque;
    
    /* Damping */
    float linear_damping;
    float angular_damping;
    
    /* Type and state */
    lg_body_type_t type;
    bool allow_sleep;
    bool is_sleeping;
    
    /* For integration */
    lg_vec3_t prev_position;
    float sleep_timer;
} lg_body_t;

/*============================================================================
 * Material Properties
 *===========================================================================*/

typedef struct {
    float restitution;      /* Bounciness, 0-1 */
    float friction;         /* Static and dynamic friction */
    float rolling_friction; /* For rolling objects */
} lg_material_t;

static const lg_material_t LG_MATERIAL_DEFAULT = {
    .restitution = 0.3f,
    .friction = 0.5f,
    .rolling_friction = 0.1f
};

/*============================================================================
 * Body Construction
 *===========================================================================*/

static inline lg_body_t lg_body(float mass) {
    lg_body_t b = {0};
    b.mass = mass;
    b.inv_mass = (mass > 0.0f) ? 1.0f / mass : 0.0f;
    b.inertia = lg_vec3_one();
    b.inv_inertia = (mass > 0.0f) ? lg_vec3_one() : lg_vec3_zero();
    b.type = (mass > 0.0f) ? LG_BODY_DYNAMIC : LG_BODY_STATIC;
    b.allow_sleep = true;
    b.linear_damping = 0.01f;
    b.angular_damping = 0.01f;
    return b;
}

static inline lg_body_t lg_body_static(void) {
    return lg_body(0.0f);
}

static inline lg_body_t lg_body_kinematic(void) {
    lg_body_t b = lg_body(1.0f);
    b.type = LG_BODY_KINEMATIC;
    b.inv_mass = 0.0f;
    return b;
}

/*============================================================================
 * Body Operations
 *===========================================================================*/

static inline void lg_body_set_mass(lg_body_t* b, float mass) {
    b->mass = mass;
    b->inv_mass = (mass > 0.0f && b->type == LG_BODY_DYNAMIC) ? 1.0f / mass : 0.0f;
}

static inline void lg_body_set_inertia(lg_body_t* b, lg_vec3_t inertia) {
    b->inertia = inertia;
    b->inv_inertia.x = (inertia.x > 0.0f) ? 1.0f / inertia.x : 0.0f;
    b->inv_inertia.y = (inertia.y > 0.0f) ? 1.0f / inertia.y : 0.0f;
    b->inv_inertia.z = (inertia.z > 0.0f) ? 1.0f / inertia.z : 0.0f;
}

static inline void lg_body_apply_force(lg_body_t* b, lg_vec3_t force) {
    if (b->type != LG_BODY_DYNAMIC) return;
    b->force = lg_vec3_add(b->force, force);
    b->is_sleeping = false;
    b->sleep_timer = 0.0f;
}

static inline void lg_body_apply_force_at(lg_body_t* b, lg_vec3_t force, lg_vec3_t point, lg_vec3_t center) {
    if (b->type != LG_BODY_DYNAMIC) return;
    b->force = lg_vec3_add(b->force, force);
    lg_vec3_t lever = lg_vec3_sub(point, center);
    b->torque = lg_vec3_add(b->torque, lg_vec3_cross(lever, force));
    b->is_sleeping = false;
    b->sleep_timer = 0.0f;
}

static inline void lg_body_apply_torque(lg_body_t* b, lg_vec3_t torque) {
    if (b->type != LG_BODY_DYNAMIC) return;
    b->torque = lg_vec3_add(b->torque, torque);
    b->is_sleeping = false;
    b->sleep_timer = 0.0f;
}

static inline void lg_body_apply_impulse(lg_body_t* b, lg_vec3_t impulse) {
    if (b->type != LG_BODY_DYNAMIC) return;
    b->velocity = lg_vec3_add(b->velocity, lg_vec3_scale(impulse, b->inv_mass));
    b->is_sleeping = false;
    b->sleep_timer = 0.0f;
}

static inline void lg_body_apply_impulse_at(lg_body_t* b, lg_vec3_t impulse, lg_vec3_t point, lg_vec3_t center) {
    if (b->type != LG_BODY_DYNAMIC) return;
    b->velocity = lg_vec3_add(b->velocity, lg_vec3_scale(impulse, b->inv_mass));
    lg_vec3_t lever = lg_vec3_sub(point, center);
    lg_vec3_t angular_impulse = lg_vec3_cross(lever, impulse);
    b->angular_velocity = lg_vec3_add(b->angular_velocity, lg_vec3_mul(angular_impulse, b->inv_inertia));
    b->is_sleeping = false;
    b->sleep_timer = 0.0f;
}

static inline void lg_body_apply_angular_impulse(lg_body_t* b, lg_vec3_t impulse) {
    if (b->type != LG_BODY_DYNAMIC) return;
    b->angular_velocity = lg_vec3_add(b->angular_velocity, lg_vec3_mul(impulse, b->inv_inertia));
    b->is_sleeping = false;
    b->sleep_timer = 0.0f;
}

static inline void lg_body_clear_forces(lg_body_t* b) {
    b->force = lg_vec3_zero();
    b->torque = lg_vec3_zero();
}

static inline void lg_body_wake(lg_body_t* b) {
    b->is_sleeping = false;
    b->sleep_timer = 0.0f;
}

static inline void lg_body_sleep(lg_body_t* b) {
    if (b->allow_sleep) {
        b->is_sleeping = true;
        b->velocity = lg_vec3_zero();
        b->angular_velocity = lg_vec3_zero();
    }
}

/*============================================================================
 * Kinetic Energy
 *===========================================================================*/

static inline float lg_body_kinetic_energy(const lg_body_t* b) {
    float linear = 0.5f * b->mass * lg_vec3_len_sq(b->velocity);
    lg_vec3_t w = b->angular_velocity;
    lg_vec3_t i = b->inertia;
    float angular = 0.5f * (i.x * w.x * w.x + i.y * w.y * w.y + i.z * w.z * w.z);
    return linear + angular;
}

static inline float lg_body_linear_speed(const lg_body_t* b) {
    return lg_vec3_len(b->velocity);
}

static inline float lg_body_angular_speed(const lg_body_t* b) {
    return lg_vec3_len(b->angular_velocity);
}

/*============================================================================
 * Sleep Check
 *===========================================================================*/

static inline bool lg_body_should_sleep(const lg_body_t* b, float threshold) {
    if (!b->allow_sleep || b->type != LG_BODY_DYNAMIC) return false;
    float linear = lg_vec3_len_sq(b->velocity);
    float angular = lg_vec3_len_sq(b->angular_velocity);
    return (linear < threshold && angular < threshold);
}

#ifdef __cplusplus
}
#endif

#endif /* LAGRANGE_BODY_H */
