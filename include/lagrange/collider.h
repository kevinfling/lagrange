/*
 * Lagrange Physics Library - Collider Component
 * Collision shapes and geometry
 */

#ifndef LAGRANGE_COLLIDER_H
#define LAGRANGE_COLLIDER_H

#include "types.h"
#include "math.h"

#ifdef __cplusplus
extern "C" {
#endif

/*============================================================================
 * Collider Shape Types
 *===========================================================================*/

typedef enum {
    LG_SHAPE_SPHERE,
    LG_SHAPE_BOX,
    LG_SHAPE_CAPSULE,
    LG_SHAPE_CYLINDER,
    LG_SHAPE_PLANE
} lg_shape_type_t;

/*============================================================================
 * Collider Component
 *===========================================================================*/

typedef struct {
    lg_shape_type_t type;
    bool is_trigger;          /* No collision response, just detection */
    float margin;             /* Collision margin for stability */
    
    union {
        /* Sphere: radius */
        struct { float radius; } sphere;
        
        /* Box: half extents (width/2, height/2, depth/2) */
        struct { lg_vec3_t half_extents; } box;
        
        /* Capsule: radius and half height */
        struct { float radius; float half_height; } capsule;
        
        /* Cylinder: radius and half height */
        struct { float radius; float half_height; } cylinder;
        
        /* Plane: normal and distance from origin */
        struct { lg_vec3_t normal; float distance; } plane;
    };
} lg_collider_t;

/*============================================================================
 * Collider Construction
 *===========================================================================*/

static inline lg_collider_t lg_collider_sphere(float radius) {
    lg_collider_t c = {0};
    c.type = LG_SHAPE_SPHERE;
    c.margin = 0.01f;
    c.sphere.radius = radius;
    return c;
}

static inline lg_collider_t lg_collider_box(float hx, float hy, float hz) {
    lg_collider_t c = {0};
    c.type = LG_SHAPE_BOX;
    c.margin = 0.01f;
    c.box.half_extents = lg_vec3(hx, hy, hz);
    return c;
}

static inline lg_collider_t lg_collider_box_v(lg_vec3_t half_extents) {
    lg_collider_t c = {0};
    c.type = LG_SHAPE_BOX;
    c.margin = 0.01f;
    c.box.half_extents = half_extents;
    return c;
}

static inline lg_collider_t lg_collider_capsule(float radius, float height) {
    lg_collider_t c = {0};
    c.type = LG_SHAPE_CAPSULE;
    c.margin = 0.01f;
    c.capsule.radius = radius;
    c.capsule.half_height = height * 0.5f;
    return c;
}

static inline lg_collider_t lg_collider_cylinder(float radius, float height) {
    lg_collider_t c = {0};
    c.type = LG_SHAPE_CYLINDER;
    c.margin = 0.01f;
    c.cylinder.radius = radius;
    c.cylinder.half_height = height * 0.5f;
    return c;
}

static inline lg_collider_t lg_collider_plane(lg_vec3_t normal, float distance) {
    lg_collider_t c = {0};
    c.type = LG_SHAPE_PLANE;
    c.margin = 0.0f;
    c.plane.normal = lg_vec3_norm(normal);
    c.plane.distance = distance;
    return c;
}

/*============================================================================
 * Collider Properties
 *===========================================================================*/

static inline float lg_collider_volume(const lg_collider_t* c) {
    switch (c->type) {
        case LG_SHAPE_SPHERE: {
            float r = c->sphere.radius;
            return (4.0f / 3.0f) * (float)M_PI * r * r * r;
        }
        case LG_SHAPE_BOX: {
            lg_vec3_t h = c->box.half_extents;
            return 8.0f * h.x * h.y * h.z;
        }
        case LG_SHAPE_CAPSULE: {
            float r = c->capsule.radius;
            float hh = c->capsule.half_height;
            float sphere_part = (4.0f / 3.0f) * (float)M_PI * r * r * r;
            float cylinder_part = (float)M_PI * r * r * (2.0f * hh);
            return sphere_part + cylinder_part;
        }
        case LG_SHAPE_CYLINDER: {
            float r = c->cylinder.radius;
            float hh = c->cylinder.half_height;
            return (float)M_PI * r * r * (2.0f * hh);
        }
        case LG_SHAPE_PLANE:
            return 0.0f; /* Infinite */
    }
    return 0.0f;
}

/* Compute inertia tensor for a solid shape with given mass */
static inline lg_vec3_t lg_collider_inertia(const lg_collider_t* c, float mass) {
    switch (c->type) {
        case LG_SHAPE_SPHERE: {
            float r = c->sphere.radius;
            float i = 0.4f * mass * r * r;
            return lg_vec3(i, i, i);
        }
        case LG_SHAPE_BOX: {
            float x = c->box.half_extents.x * 2.0f;
            float y = c->box.half_extents.y * 2.0f;
            float z = c->box.half_extents.z * 2.0f;
            float ix = (1.0f / 12.0f) * mass * (y * y + z * z);
            float iy = (1.0f / 12.0f) * mass * (x * x + z * z);
            float iz = (1.0f / 12.0f) * mass * (x * x + y * y);
            return lg_vec3(ix, iy, iz);
        }
        case LG_SHAPE_CAPSULE: {
            float r = c->capsule.radius;
            float hh = c->capsule.half_height;
            // Approximate as cylinder + point masses at ends
            float ixz = 0.5f * mass * r * r + (1.0f / 3.0f) * mass * hh * hh;
            float iy = mass * r * r;
            return lg_vec3(ixz, iy, ixz);
        }
        case LG_SHAPE_CYLINDER: {
            float r = c->cylinder.radius;
            float hh = c->cylinder.half_height;
            float ixz = (1.0f / 12.0f) * mass * (3.0f * r * r + 4.0f * hh * hh);
            float iy = 0.5f * mass * r * r;
            return lg_vec3(ixz, iy, ixz);
        }
        case LG_SHAPE_PLANE:
            return lg_vec3_zero(); /* Infinite mass */
    }
    return lg_vec3_zero();
}

/* Compute bounding box in local space */
static inline void lg_collider_bounds(const lg_collider_t* c, lg_vec3_t* out_min, lg_vec3_t* out_max) {
    switch (c->type) {
        case LG_SHAPE_SPHERE: {
            float r = c->sphere.radius + c->margin;
            *out_min = lg_vec3(-r, -r, -r);
            *out_max = lg_vec3(r, r, r);
            break;
        }
        case LG_SHAPE_BOX: {
            lg_vec3_t h = c->box.half_extents;
            float m = c->margin;
            *out_min = lg_vec3(-h.x - m, -h.y - m, -h.z - m);
            *out_max = lg_vec3(h.x + m, h.y + m, h.z + m);
            break;
        }
        case LG_SHAPE_CAPSULE: {
            float r = c->capsule.radius + c->margin;
            float hh = c->capsule.half_height;
            *out_min = lg_vec3(-r, -hh - r, -r);
            *out_max = lg_vec3(r, hh + r, r);
            break;
        }
        case LG_SHAPE_CYLINDER: {
            float r = c->cylinder.radius + c->margin;
            float hh = c->cylinder.half_height;
            *out_min = lg_vec3(-r, -hh, -r);
            *out_max = lg_vec3(r, hh, r);
            break;
        }
        case LG_SHAPE_PLANE: {
            *out_min = lg_vec3(-1e30f, -1e30f, -1e30f);
            *out_max = lg_vec3(1e30f, 1e30f, 1e30f);
            break;
        }
    }
}

#ifdef __cplusplus
}
#endif

#endif /* LAGRANGE_COLLIDER_H */
