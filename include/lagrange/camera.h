#pragma once
#include "types.h"
#include <stdbool.h>

typedef enum {
    LG_PROJ_PERSPECTIVE,
    LG_PROJ_ORTHO
} lg_proj_type_t;

typedef struct {
    lg_vec3d_t    position;     // world-space (double for orbital scale)
    lg_quat_t     rotation;     // orientation (no gimbal lock)
    double        fov_y;        // radians, ignored for ortho
    double        aspect;
    double        near_z;
    double        far_z;
    double        ortho_height; // for LG_PROJ_ORTHO
    lg_proj_type_t proj_type;

    lg_mat4_t     view;         // cached
    lg_mat4_t     proj;         // cached
    lg_mat4_t     view_proj;    // cached
    bool          dirty;
} _Alignas(64) lg_camera_t;

// Initialize with sensible defaults (orbital viewer)
static inline void lg_camera_init(lg_camera_t *restrict cam) {
    *cam = (lg_camera_t){
        .position = {0, 0, 10},
        .rotation = {0, 0, 0, 1},  // identity quat
        .fov_y = 1.0,               // ~57° vertical
        .aspect = 16.0/9.0,
        .near_z = 0.1,
        .far_z = 1e9,               // orbital scale
        .ortho_height = 20.0,
        .proj_type = LG_PROJ_PERSPECTIVE,
        .dirty = true
    };
    lg_camera_update(cam);  // compute initial matrices
}

// Recompute matrices only when dirty — cache hot path
static inline void lg_camera_update(lg_camera_t *restrict cam) {
    if (!cam->dirty) return;

    // view = inverse(transform) — you already have lg_mat4_from_trs or similar in transform.h
    lg_mat4_t trs = /* your existing transform compose */;
    cam->view = lg_mat4_inverse(&trs);  // or pre-multiply if you store it inverted

    if (cam->proj_type == LG_PROJ_PERSPECTIVE) {
        cam->proj = lg_mat4_perspective(cam->fov_y, cam->aspect, cam->near_z, cam->far_z);
    } else {
        cam->proj = lg_mat4_ortho(-cam->aspect*cam->ortho_height*0.5,
                                   cam->aspect*cam->ortho_height*0.5,
                                  -cam->ortho_height*0.5, cam->ortho_height*0.5,
                                   cam->near_z, cam->far_z);
    }

    cam->view_proj = lg_mat4_mul(&cam->view, &cam->proj);
    cam->dirty = false;
}

// Project world point → terminal (x,y) in [0, w)×[0, h) plus depth for z-order/color
static inline bool lg_camera_project(const lg_camera_t *restrict cam,
                                     lg_vec3d_t world,
                                     size_t term_w, size_t term_h,
                                     float *restrict screen_x,
                                     float *restrict screen_y,
                                     float *restrict depth) {
    lg_vec4_t clip = lg_mat4_transform_vec4(&cam->view_proj, (lg_vec4_t){world.x, world.y, world.z, 1.0});
    if (clip.w <= 0.0) return false;  // behind camera

    lg_vec3_t ndc = {clip.x/clip.w, clip.y/clip.w, clip.z/clip.w};
    *screen_x = (ndc.x * 0.5f + 0.5f) * (float)term_w;
    *screen_y = (1.0f - (ndc.y * 0.5f + 0.5f)) * (float)term_h;  // terminal Y-down
    *depth    = (float)ndc.z;
    return true;
}

// Convenience helpers you’ll thank me for at 3 a.m.
static inline void lg_camera_look_at(lg_camera_t *restrict cam, lg_vec3d_t target, lg_vec3d_t up) {
    // build rotation quat from forward/up — you already have the math
    lg_vec3d_t forward = lg_vec3d_normalize(lg_vec3d_sub(target, cam->position));
    // ... quat from two vectors (reuse your existing code)
    cam->dirty = true;
}

static inline void lg_camera_orbit(lg_camera_t *restrict cam, lg_vec3d_t center, double azimuth, double elevation) {
    // classic orbital camera for your patched-conic demos
    cam->position = /* spherical offset from center */;
    cam->rotation = /* look-at quat toward center */;
    cam->dirty = true;
}

