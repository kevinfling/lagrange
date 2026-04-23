/*
 * Lagrange Physics Library - Transform Component
 * Position, rotation, and scale for entities
 */

#ifndef LAGRANGE_TRANSFORM_H
#define LAGRANGE_TRANSFORM_H

#include "types.h"
#include "math.h"

#ifdef __cplusplus
extern "C" {
#endif

/*============================================================================
 * Transform Component
 *===========================================================================*/

typedef struct {
    lg_vec3_t position;
    lg_quat_t rotation;
    lg_vec3_t scale;
} lg_transform_t;

static const lg_transform_t LG_TRANSFORM_IDENTITY = {
    .position = {0.0f, 0.0f, 0.0f},
    .rotation = {0.0f, 0.0f, 0.0f, 1.0f},
    .scale = {1.0f, 1.0f, 1.0f}
};

/* Construction */
static inline lg_transform_t lg_transform(void) {
    return LG_TRANSFORM_IDENTITY;
}

static inline lg_transform_t lg_transform_at(lg_vec3_t position) {
    lg_transform_t t = LG_TRANSFORM_IDENTITY;
    t.position = position;
    return t;
}

static inline lg_transform_t lg_transform_full(
    lg_vec3_t position,
    lg_quat_t rotation,
    lg_vec3_t scale
) {
    lg_transform_t t = {position, rotation, scale};
    return t;
}

/* Direction vectors */
static inline lg_vec3_t lg_transform_forward(lg_transform_t t) {
    return lg_quat_rotate(t.rotation, lg_vec3_forward());
}

static inline lg_vec3_t lg_transform_right(lg_transform_t t) {
    return lg_quat_rotate(t.rotation, lg_vec3_right());
}

static inline lg_vec3_t lg_transform_up(lg_transform_t t) {
    return lg_quat_rotate(t.rotation, lg_vec3_up());
}

/* Transform point from local to world space */
static inline lg_vec3_t lg_transform_point(lg_transform_t t, lg_vec3_t local) {
    lg_vec3_t scaled = lg_vec3_mul(local, t.scale);
    lg_vec3_t rotated = lg_quat_rotate(t.rotation, scaled);
    return lg_vec3_add(rotated, t.position);
}

/* Transform point from world to local space */
static inline lg_vec3_t lg_transform_inv_point(lg_transform_t t, lg_vec3_t world) {
    lg_vec3_t translated = lg_vec3_sub(world, t.position);
    lg_quat_t inv_rot = lg_quat_conj(t.rotation);
    lg_vec3_t rotated = lg_quat_rotate(inv_rot, translated);
    return lg_vec3_div(rotated, t.scale);
}

/* Transform direction from local to world space */
static inline lg_vec3_t lg_transform_direction(lg_transform_t t, lg_vec3_t local) {
    return lg_quat_rotate(t.rotation, local);
}

/* Convert to matrix */
static inline lg_mat4_t lg_transform_to_matrix(lg_transform_t t) {
    lg_mat4_t translation = lg_mat4_translation(t.position);
    lg_mat4_t rotation = lg_mat4_rotation(t.rotation);
    lg_mat4_t scale = lg_mat4_scale(t.scale);
    
    // Combine: T * R * S
    lg_mat4_t result;
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            float sum = 0.0f;
            for (int k = 0; k < 4; k++) {
                sum += translation.m[i + k * 4] * rotation.m[k + j * 4];
            }
            result.m[i + j * 4] = sum;
        }
    }
    
    lg_mat4_t final;
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            float sum = 0.0f;
            for (int k = 0; k < 4; k++) {
                sum += result.m[i + k * 4] * scale.m[k + j * 4];
            }
            final.m[i + j * 4] = sum;
        }
    }
    
    return final;
}

/* Interpolation */
static inline lg_transform_t lg_transform_lerp(lg_transform_t a, lg_transform_t b, float t) {
    return lg_transform_full(
        lg_vec3_lerp(a.position, b.position, t),
        lg_quat_slerp(a.rotation, b.rotation, t),
        lg_vec3_lerp(a.scale, b.scale, t)
    );
}

#ifdef __cplusplus
}
#endif

#endif /* LAGRANGE_TRANSFORM_H */
