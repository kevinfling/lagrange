/**
 * lagrange_spatial_index.h - libspatial Integration
 *
 * Provides BVH-accelerated broad-phase collision detection and octree-based
 * spatial queries for LaGrange worlds, using libspatial as the underlying
 * spatial indexing engine.
 *
 * libspatial is configured for float coordinates to match LaGrange's
 * game-physics precision model.
 */

#ifndef LAGRANGE_SPATIAL_INDEX_H
#define LAGRANGE_SPATIAL_INDEX_H

/* Configure libspatial to match LaGrange coordinate types */
#define SPATIAL_NUMTYPE float

/* Suppress libspatial's cast-align warnings (third-party code) */
#if defined(__clang__) || defined(__GNUC__)
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wcast-align"
#endif

#include <libspatial/bvh.h>
#include <libspatial/octree.h>

#if defined(__clang__) || defined(__GNUC__)
    #pragma GCC diagnostic pop
#endif

#include "world.h"

#ifdef __cplusplus
extern "C" {
#endif

/*============================================================================
 * Collider World-Space AABB
 *===========================================================================*/

static inline void lg_collider_world_aabb(const lg_collider_t* c, const lg_transform_t* t,
                                          float out_min[3], float out_max[3]) {
    lg_vec3_t local_min, local_max;
    lg_collider_bounds(c, &local_min, &local_max);

    /* For infinite planes, clamp to a large but finite bound */
    if (c->type == LG_SHAPE_PLANE) {
        out_min[0] = -1e15f; out_min[1] = -1e15f; out_min[2] = -1e15f;
        out_max[0] =  1e15f; out_max[1] =  1e15f; out_max[2] =  1e15f;
        return;
    }

    /* Transform 8 corners and find extremes */
    out_min[0] = out_min[1] = out_min[2] = 1e30f;
    out_max[0] = out_max[1] = out_max[2] = -1e30f;

    for (int corner = 0; corner < 8; corner++) {
        lg_vec3_t local = lg_vec3(
            (corner & 1) ? local_max.x : local_min.x,
            (corner & 2) ? local_max.y : local_min.y,
            (corner & 4) ? local_max.z : local_min.z
        );
        lg_vec3_t world = lg_vec3_add(t->position, lg_quat_rotate(t->rotation, local));

        out_min[0] = fminf(out_min[0], world.x);
        out_min[1] = fminf(out_min[1], world.y);
        out_min[2] = fminf(out_min[2], world.z);
        out_max[0] = fmaxf(out_max[0], world.x);
        out_max[1] = fmaxf(out_max[1], world.y);
        out_max[2] = fmaxf(out_max[2], world.z);
    }
}

/*============================================================================
 * BVH-Accelerated Broad Phase
 *===========================================================================*/

typedef struct {
    lg_world_t* world;
    lg_contact_t* contacts;
    int max_contacts;
    int* out_count;
    size_t query_idx;
} lg_bvh_query_state_t;

static inline bool lg_bvh_broad_phase_callback(const float* min, const float* max,
                                                void* data, void* udata) {
    (void)min;
    (void)max;
    lg_bvh_query_state_t* state = (lg_bvh_query_state_t*)udata;
    size_t j = (size_t)(uintptr_t)data;

    /* Only process pairs where j > query_idx to avoid duplicates */
    if (j <= state->query_idx) return true;

    if (*state->out_count >= state->max_contacts) return false; /* stop */

    lg_contact_t contact = {0};
    if (lg_narrow_phase(&state->world->storage, state->query_idx, j, &contact)) {
        state->contacts[(*state->out_count)++] = contact;
    }
    return true;
}

/**
 * BVH-accelerated broad-phase collision detection.
 *
 * For small entity counts (< 32), falls back to brute-force O(N^2).
 * For larger counts, builds a BVH from collider AABBs and only tests
 * overlapping pairs, giving O(N log N) behavior.
 */
static inline void lg_broad_phase_bvh(lg_world_t* world, lg_contact_t* contacts,
                                       int max_contacts, int* out_count) {
    *out_count = 0;
    lg_storage_t* s = &world->storage;
    size_t n = s->count;

    if (n < 2) return;

    /* Brute force is faster for tiny worlds (no BVH build overhead) */
    if (n < 32) {
        lg_broad_phase(world, contacts, max_contacts, out_count);
        return;
    }

    /* Build BVH from all collider AABBs */
    spatial_bvh* bvh = spatial_bvh_new();
    if (!bvh) return;

    float min[3], max[3];
    for (size_t i = 0; i < n; i++) {
        lg_collider_t* c = &s->colliders[i];
        if (c->is_trigger) continue;

        lg_collider_world_aabb(c, &s->transforms[i], min, max);
        spatial_bvh_insert(bvh, min, max, (void*)(uintptr_t)i);
    }

    if (!spatial_bvh_build(bvh)) {
        spatial_bvh_free(bvh);
        return;
    }

    /* Query each collider against the BVH */
    lg_bvh_query_state_t state = {
        .world = world,
        .contacts = contacts,
        .max_contacts = max_contacts,
        .out_count = out_count,
        .query_idx = 0
    };

    for (size_t i = 0; i < n && *out_count < max_contacts; i++) {
        lg_collider_t* c = &s->colliders[i];
        if (c->is_trigger) continue;

        lg_collider_world_aabb(c, &s->transforms[i], min, max);
        state.query_idx = i;
        spatial_bvh_search(bvh, min, max, lg_bvh_broad_phase_callback, &state);
    }

    spatial_bvh_free(bvh);
}

/*============================================================================
 * Octree Spatial Queries for Worlds
 *===========================================================================*/

typedef struct {
    lg_world_t* world;
    lg_entity_callback_t callback;
    void* user_data;
} lg_octree_query_state_t;

static inline bool lg_octree_entity_callback(const float* min, const float* max,
                                              void* data, void* udata) {
    (void)min;
    (void)max;
    lg_octree_query_state_t* state = (lg_octree_query_state_t*)udata;
    size_t idx = (size_t)(uintptr_t)data;
    state->callback(state->world, state->world->storage.entities[idx], state->user_data);
    return true;
}

/**
 * Query all entities whose collider AABB overlaps the given AABB.
 * Uses an octree for O(log N) query time.
 */
static inline void lg_world_query_aabb(lg_world_t* world,
                                        float qmin[3], float qmax[3],
                                        lg_entity_callback_t callback,
                                        void* user_data) {
    if (!world || !callback) return;

    lg_storage_t* s = &world->storage;
    size_t n = s->count;
    if (n == 0) return;

    /* Small worlds: brute force */
    if (n < 64) {
        for (size_t i = 0; i < n; i++) {
            float emin[3], emax[3];
            lg_collider_world_aabb(&s->colliders[i], &s->transforms[i], emin, emax);
            /* AABB overlap test */
            if (emin[0] <= qmax[0] && emax[0] >= qmin[0] &&
                emin[1] <= qmax[1] && emax[1] >= qmin[1] &&
                emin[2] <= qmax[2] && emax[2] >= qmin[2]) {
                callback(world, s->entities[i], user_data);
            }
        }
        return;
    }

    spatial_octree* oct = spatial_octree_new();
    if (!oct) return;

    float min[3], max[3];
    for (size_t i = 0; i < n; i++) {
        lg_collider_world_aabb(&s->colliders[i], &s->transforms[i], min, max);
        spatial_octree_insert(oct, min, max, (void*)(uintptr_t)i);
    }

    lg_octree_query_state_t state = {
        .world = world,
        .callback = callback,
        .user_data = user_data
    };
    spatial_octree_search(oct, qmin, qmax, lg_octree_entity_callback, &state);
    spatial_octree_free(oct);
}

/**
 * Query all entities whose collider AABB overlaps the given sphere.
 */
static inline void lg_world_query_sphere(lg_world_t* world,
                                          lg_vec3_t center, float radius,
                                          lg_entity_callback_t callback,
                                          void* user_data) {
    float qmin[3] = { center.x - radius, center.y - radius, center.z - radius };
    float qmax[3] = { center.x + radius, center.y + radius, center.z + radius };
    lg_world_query_aabb(world, qmin, qmax, callback, user_data);
}

#ifdef __cplusplus
}
#endif

#endif /* LAGRANGE_SPATIAL_INDEX_H */
