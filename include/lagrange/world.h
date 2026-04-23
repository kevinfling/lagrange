/*
 * Lagrange Physics Library - World API
 * Core simulation context and entity management
 */

#ifndef LAGRANGE_WORLD_H
#define LAGRANGE_WORLD_H

#include "types.h"
#include "math.h"
#include "transform.h"
#include "body.h"
#include "collider.h"
#include "integrator.h"

#include <stdlib.h>
#include <string.h>
#include <assert.h>

#ifdef __cplusplus
extern "C" {
#endif

/*============================================================================
 * Internal Component Storage
 *===========================================================================*/

typedef struct {
    lg_entity_t* entities;
    lg_transform_t* transforms;
    lg_transform_t* prev_transforms;
    lg_body_t* bodies;
    lg_collider_t* colliders;
    lg_material_t* materials;
    lg_verlet_state_t* verlet_states;
    
    size_t count;
    size_t capacity;
} lg_storage_t;

static inline bool lg_storage_init(lg_storage_t* s, size_t capacity) {
    s->entities = (lg_entity_t*)calloc(capacity, sizeof(lg_entity_t));
    s->transforms = (lg_transform_t*)calloc(capacity, sizeof(lg_transform_t));
    s->prev_transforms = (lg_transform_t*)calloc(capacity, sizeof(lg_transform_t));
    s->bodies = (lg_body_t*)calloc(capacity, sizeof(lg_body_t));
    s->colliders = (lg_collider_t*)calloc(capacity, sizeof(lg_collider_t));
    s->materials = (lg_material_t*)calloc(capacity, sizeof(lg_material_t));
    s->verlet_states = (lg_verlet_state_t*)calloc(capacity, sizeof(lg_verlet_state_t));
    
    if (!s->entities || !s->transforms || !s->prev_transforms || !s->bodies || 
        !s->colliders || !s->materials || !s->verlet_states) {
        return false;
    }
    
    s->count = 0;
    s->capacity = capacity;
    return true;
}

static inline void lg_storage_free(lg_storage_t* s) {
    free(s->entities);
    free(s->transforms);
    free(s->prev_transforms);
    free(s->bodies);
    free(s->colliders);
    free(s->materials);
    free(s->verlet_states);
    s->count = 0;
    s->capacity = 0;
}

static inline int lg_storage_find(lg_storage_t* s, lg_entity_t entity) {
    for (size_t i = 0; i < s->count; i++) {
        if (s->entities[i] == entity) return (int)i;
    }
    return -1;
}

static inline int lg_storage_add(lg_storage_t* s, lg_entity_t entity) {
    if (s->count >= s->capacity) return -1;
    
    size_t idx = s->count++;
    s->entities[idx] = entity;
    s->transforms[idx] = LG_TRANSFORM_IDENTITY;
    s->bodies[idx] = lg_body(1.0f);
    s->colliders[idx] = lg_collider_sphere(0.5f);
    s->materials[idx] = LG_MATERIAL_DEFAULT;
    s->verlet_states[idx] = (lg_verlet_state_t){.first_step = true};
    
    return (int)idx;
}

static inline void lg_storage_remove(lg_storage_t* s, size_t idx) {
    if (idx >= s->count) return;
    
    size_t last = s->count - 1;
    if (idx != last) {
        s->entities[idx] = s->entities[last];
        s->transforms[idx] = s->transforms[last];
        s->prev_transforms[idx] = s->prev_transforms[last];
        s->bodies[idx] = s->bodies[last];
        s->colliders[idx] = s->colliders[last];
        s->materials[idx] = s->materials[last];
        s->verlet_states[idx] = s->verlet_states[last];
    }
    s->count--;
}

/*============================================================================
 * World Structure
 *===========================================================================*/

struct lg_world_t {
    lg_world_config_t config;
    lg_storage_t storage;
    
    /* Simulation state */
    float time;
    float accumulator;
    uint64_t step_count;
    
    /* Performance */
    double last_step_time;
    
    /* Settings */
    lg_integrator_type_t integrator;
    lg_vec3_t gravity;
    bool use_gravity;
    
    /* Entity ID generation and recycling */
    uint64_t next_id;
    lg_entity_t* free_list;
    size_t free_list_count;
    size_t free_list_capacity;
    
    /* Collision callbacks */
    lg_collision_callback_t collision_callback;
    void* collision_callback_user_data;
};

/*============================================================================
 * World Lifecycle
 *===========================================================================*/

static inline lg_world_t* lg_world_create(const lg_world_config_t* config) {
    lg_world_t* world = (lg_world_t*)calloc(1, sizeof(lg_world_t));
    if (!world) return NULL;
    
    if (config) {
        world->config = *config;
    } else {
        world->config = LG_WORLD_CONFIG_DEFAULT;
    }
    
    if (!lg_storage_init(&world->storage, world->config.max_entities)) {
        free(world);
        return NULL;
    }
    
    world->time = 0.0f;
    world->accumulator = 0.0f;
    world->step_count = 0;
    world->integrator = LG_INTEGRATOR_SEMI_IMPLICIT_EULER;
    world->gravity = lg_vec3(world->config.gravity[0], 
                              world->config.gravity[1], 
                              world->config.gravity[2]);
    world->use_gravity = true;
    world->next_id = 1;
    
    return world;
}

static inline void lg_world_destroy(lg_world_t* world) {
    if (!world) return;
    lg_storage_free(&world->storage);
    free(world->free_list);
    free(world);
}

static inline void lg_world_reset(lg_world_t* world) {
    if (!world) return;
    world->storage.count = 0;
    world->time = 0.0f;
    world->accumulator = 0.0f;
    world->step_count = 0;
    world->next_id = 1;
    world->free_list_count = 0;
    memset(world->storage.prev_transforms, 0, world->storage.capacity * sizeof(lg_transform_t));
}

/*============================================================================
 * Entity Management
 *===========================================================================*/

static inline lg_entity_t lg_entity_create(lg_world_t* world) {
    if (!world) return LG_ENTITY_INVALID;
    
    lg_entity_t entity;
    if (world->free_list_count > 0) {
        entity = world->free_list[--world->free_list_count];
    } else {
        entity = world->next_id++;
    }
    
    int idx = lg_storage_add(&world->storage, entity);
    if (idx < 0) {
        /* Push ID back onto free list if storage is full */
        if (world->free_list_count >= world->free_list_capacity) {
            size_t new_cap = world->free_list_capacity ? world->free_list_capacity * 2 : 16;
            lg_entity_t* new_list = (lg_entity_t*)realloc(world->free_list, new_cap * sizeof(lg_entity_t));
            if (!new_list) return LG_ENTITY_INVALID;
            world->free_list = new_list;
            world->free_list_capacity = new_cap;
        }
        world->free_list[world->free_list_count++] = entity;
        return LG_ENTITY_INVALID;
    }
    
    /* Initialize previous transform for interpolation */
    world->storage.prev_transforms[idx] = world->storage.transforms[idx];
    
    return entity;
}

static inline void lg_entity_destroy(lg_world_t* world, lg_entity_t entity) {
    if (!world || entity == LG_ENTITY_INVALID) return;
    
    int idx = lg_storage_find(&world->storage, entity);
    if (idx >= 0) {
        lg_storage_remove(&world->storage, (size_t)idx);
        
        if (world->free_list_count >= world->free_list_capacity) {
            size_t new_cap = world->free_list_capacity ? world->free_list_capacity * 2 : 16;
            lg_entity_t* new_list = (lg_entity_t*)realloc(world->free_list, new_cap * sizeof(lg_entity_t));
            if (!new_list) return;
            world->free_list = new_list;
            world->free_list_capacity = new_cap;
        }
        world->free_list[world->free_list_count++] = entity;
    }
}

static inline bool lg_entity_valid(lg_world_t* world, lg_entity_t entity) {
    if (!world || entity == LG_ENTITY_INVALID) return false;
    return lg_storage_find(&world->storage, entity) >= 0;
}

/*============================================================================
 * Component Access
 *===========================================================================*/

static inline lg_transform_t* lg_get_transform(lg_world_t* world, lg_entity_t entity) {
    int idx = lg_storage_find(&world->storage, entity);
    if (idx < 0) return NULL;
    return &world->storage.transforms[idx];
}

static inline lg_body_t* lg_get_body(lg_world_t* world, lg_entity_t entity) {
    int idx = lg_storage_find(&world->storage, entity);
    if (idx < 0) return NULL;
    return &world->storage.bodies[idx];
}

static inline lg_collider_t* lg_get_collider(lg_world_t* world, lg_entity_t entity) {
    int idx = lg_storage_find(&world->storage, entity);
    if (idx < 0) return NULL;
    return &world->storage.colliders[idx];
}

static inline lg_material_t* lg_get_material(lg_world_t* world, lg_entity_t entity) {
    int idx = lg_storage_find(&world->storage, entity);
    if (idx < 0) return NULL;
    return &world->storage.materials[idx];
}

/*============================================================================
 * Component Setters (Convenience)
 *===========================================================================*/

static inline void lg_set_transform(lg_world_t* world, lg_entity_t entity, const lg_transform_t* transform) {
    lg_transform_t* t = lg_get_transform(world, entity);
    if (t) *t = *transform;
}

static inline void lg_set_body(lg_world_t* world, lg_entity_t entity, const lg_body_t* body) {
    lg_body_t* b = lg_get_body(world, entity);
    if (b) {
        *b = *body;
        /* Update inverse mass based on type */
        if (b->type != LG_BODY_DYNAMIC) {
            b->inv_mass = 0.0f;
        } else {
            b->inv_mass = (b->mass > 0.0f) ? 1.0f / b->mass : 0.0f;
        }
    }
}

static inline void lg_set_collider(lg_world_t* world, lg_entity_t entity, const lg_collider_t* collider) {
    lg_collider_t* c = lg_get_collider(world, entity);
    if (c) *c = *collider;
}

static inline void lg_set_material(lg_world_t* world, lg_entity_t entity, const lg_material_t* material) {
    lg_material_t* m = lg_get_material(world, entity);
    if (m) *m = *material;
}

/*============================================================================
 * Position/Rotation/Velocity Shortcuts
 *===========================================================================*/

static inline lg_vec3_t lg_get_position(lg_world_t* world, lg_entity_t entity) {
    lg_transform_t* t = lg_get_transform(world, entity);
    return t ? t->position : lg_vec3_zero();
}

static inline void lg_set_position(lg_world_t* world, lg_entity_t entity, lg_vec3_t pos) {
    lg_transform_t* t = lg_get_transform(world, entity);
    if (t) t->position = pos;
}

static inline lg_quat_t lg_get_rotation(lg_world_t* world, lg_entity_t entity) {
    lg_transform_t* t = lg_get_transform(world, entity);
    return t ? t->rotation : lg_quat_identity();
}

static inline void lg_set_rotation(lg_world_t* world, lg_entity_t entity, lg_quat_t rot) {
    lg_transform_t* t = lg_get_transform(world, entity);
    if (t) t->rotation = rot;
}

static inline lg_vec3_t lg_get_velocity(lg_world_t* world, lg_entity_t entity) {
    lg_body_t* b = lg_get_body(world, entity);
    return b ? b->velocity : lg_vec3_zero();
}

static inline void lg_set_velocity(lg_world_t* world, lg_entity_t entity, lg_vec3_t vel) {
    lg_body_t* b = lg_get_body(world, entity);
    if (b) {
        b->velocity = vel;
        lg_body_wake(b);
    }
}

static inline float lg_get_mass(lg_world_t* world, lg_entity_t entity) {
    lg_body_t* b = lg_get_body(world, entity);
    return b ? b->mass : 0.0f;
}

static inline void lg_set_mass(lg_world_t* world, lg_entity_t entity, float mass) {
    lg_body_t* b = lg_get_body(world, entity);
    if (b) lg_body_set_mass(b, mass);
}

/*============================================================================
 * Force Application Shortcuts
 *===========================================================================*/

static inline void lg_apply_force(lg_world_t* world, lg_entity_t entity, lg_vec3_t force) {
    lg_body_t* b = lg_get_body(world, entity);
    if (b) lg_body_apply_force(b, force);
}

static inline void lg_apply_impulse(lg_world_t* world, lg_entity_t entity, lg_vec3_t impulse) {
    lg_body_t* b = lg_get_body(world, entity);
    if (b) lg_body_apply_impulse(b, impulse);
}

static inline void lg_apply_torque(lg_world_t* world, lg_entity_t entity, lg_vec3_t torque) {
    lg_body_t* b = lg_get_body(world, entity);
    if (b) lg_body_apply_torque(b, torque);
}

/*============================================================================
 * World Settings
 *===========================================================================*/

static inline void lg_world_set_gravity_v(lg_world_t* world, lg_vec3_t gravity) {
    if (world) world->gravity = gravity;
}

static inline lg_vec3_t lg_world_get_gravity(lg_world_t* world) {
    return world ? world->gravity : lg_vec3_zero();
}

static inline void lg_world_set_integrator(lg_world_t* world, lg_integrator_type_t type) {
    if (world) world->integrator = type;
}

static inline lg_integrator_type_t lg_world_get_integrator(lg_world_t* world) {
    return world ? world->integrator : LG_INTEGRATOR_SEMI_IMPLICIT_EULER;
}

/*============================================================================
 * Collision Callbacks
 *===========================================================================*/

static inline void lg_world_set_collision_callback(
    lg_world_t* world,
    lg_collision_callback_t callback,
    void* user_data
) {
    if (!world) return;
    world->collision_callback = callback;
    world->collision_callback_user_data = user_data;
}

/*============================================================================
 * Iteration
 *===========================================================================*/

typedef void (*lg_entity_callback_t)(lg_world_t* world, lg_entity_t entity, void* user_data);

static inline void lg_world_for_each(lg_world_t* world, lg_entity_callback_t callback, void* user_data) {
    if (!world || !callback) return;
    
    for (size_t i = 0; i < world->storage.count; i++) {
        callback(world, world->storage.entities[i], user_data);
    }
}

/* Get count of entities */
static inline size_t lg_world_entity_count(lg_world_t* world) {
    return world ? world->storage.count : 0;
}

/* Get entity at index (for iteration) */
static inline lg_entity_t lg_world_get_entity(lg_world_t* world, size_t index) {
    if (!world || index >= world->storage.count) return LG_ENTITY_INVALID;
    return world->storage.entities[index];
}

#ifdef __cplusplus
}
#endif

#endif /* LAGRANGE_WORLD_H */
