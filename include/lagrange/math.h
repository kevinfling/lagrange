/*
 * Lagrange Physics Library - Math Module
 * Header-only C library for 3D math operations
 * 
 * Usage:
 *   #include <lagrange/math.h>
 *   
 *   // In exactly one source file:
 *   #define LAGRANGE_IMPLEMENTATION
 *   #include <lagrange/math.h>
 */

#ifndef LAGRANGE_MATH_H
#define LAGRANGE_MATH_H

#include <math.h>
#include <stdbool.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/*============================================================================
 * 3D Vector (float)
 *===========================================================================*/

typedef struct {
    float x, y, z;
} lg_vec3_t;

typedef struct {
    float x, y, z, w;
} lg_vec4_t;

/* Constants */
static const lg_vec3_t LG_VEC3_ZERO = {0.0f, 0.0f, 0.0f};
static const lg_vec3_t LG_VEC3_ONE = {1.0f, 1.0f, 1.0f};
static const lg_vec3_t LG_VEC3_UP = {0.0f, 1.0f, 0.0f};
static const lg_vec3_t LG_VEC3_RIGHT = {1.0f, 0.0f, 0.0f};
static const lg_vec3_t LG_VEC3_FORWARD = {0.0f, 0.0f, 1.0f};

/* Construction */
static inline lg_vec3_t lg_vec3(float x, float y, float z) {
    lg_vec3_t v = {x, y, z};
    return v;
}

static inline lg_vec3_t lg_vec3_zero(void) { return LG_VEC3_ZERO; }
static inline lg_vec3_t lg_vec3_one(void) { return LG_VEC3_ONE; }
static inline lg_vec3_t lg_vec3_up(void) { return LG_VEC3_UP; }
static inline lg_vec3_t lg_vec3_right(void) { return LG_VEC3_RIGHT; }
static inline lg_vec3_t lg_vec3_forward(void) { return LG_VEC3_FORWARD; }

/* Basic operations */
static inline lg_vec3_t lg_vec3_add(lg_vec3_t a, lg_vec3_t b) {
    return lg_vec3(a.x + b.x, a.y + b.y, a.z + b.z);
}

static inline lg_vec3_t lg_vec3_sub(lg_vec3_t a, lg_vec3_t b) {
    return lg_vec3(a.x - b.x, a.y - b.y, a.z - b.z);
}

static inline lg_vec3_t lg_vec3_mul(lg_vec3_t a, lg_vec3_t b) {
    return lg_vec3(a.x * b.x, a.y * b.y, a.z * b.z);
}

static inline lg_vec3_t lg_vec3_div(lg_vec3_t a, lg_vec3_t b) {
    return lg_vec3(a.x / b.x, a.y / b.y, a.z / b.z);
}

static inline lg_vec3_t lg_vec3_scale(lg_vec3_t v, float s) {
    return lg_vec3(v.x * s, v.y * s, v.z * s);
}

static inline lg_vec3_t lg_vec3_neg(lg_vec3_t v) {
    return lg_vec3(-v.x, -v.y, -v.z);
}

/* Vector properties */
static inline float lg_vec3_dot(lg_vec3_t a, lg_vec3_t b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

static inline lg_vec3_t lg_vec3_cross(lg_vec3_t a, lg_vec3_t b) {
    return lg_vec3(
        a.y * b.z - a.z * b.y,
        a.z * b.x - a.x * b.z,
        a.x * b.y - a.y * b.x
    );
}

static inline float lg_vec3_len_sq(lg_vec3_t v) {
    return v.x * v.x + v.y * v.y + v.z * v.z;
}

static inline float lg_vec3_len(lg_vec3_t v) {
    return sqrtf(lg_vec3_len_sq(v));
}

static inline float lg_vec3_dist_sq(lg_vec3_t a, lg_vec3_t b) {
    return lg_vec3_len_sq(lg_vec3_sub(a, b));
}

static inline float lg_vec3_dist(lg_vec3_t a, lg_vec3_t b) {
    return sqrtf(lg_vec3_dist_sq(a, b));
}

static inline lg_vec3_t lg_vec3_norm(lg_vec3_t v) {
    float len = lg_vec3_len(v);
    if (len < 1e-6f) return LG_VEC3_ZERO;
    return lg_vec3_scale(v, 1.0f / len);
}

/* Interpolation */
static inline lg_vec3_t lg_vec3_lerp(lg_vec3_t a, lg_vec3_t b, float t) {
    return lg_vec3_add(a, lg_vec3_scale(lg_vec3_sub(b, a), t));
}

/* Projection/rejection */
static inline lg_vec3_t lg_vec3_proj(lg_vec3_t v, lg_vec3_t onto) {
    lg_vec3_t n = lg_vec3_norm(onto);
    return lg_vec3_scale(n, lg_vec3_dot(v, n));
}

static inline lg_vec3_t lg_vec3_rej(lg_vec3_t v, lg_vec3_t onto) {
    return lg_vec3_sub(v, lg_vec3_proj(v, onto));
}

static inline lg_vec3_t lg_vec3_reflect(lg_vec3_t v, lg_vec3_t normal) {
    lg_vec3_t n = lg_vec3_norm(normal);
    return lg_vec3_sub(v, lg_vec3_scale(n, 2.0f * lg_vec3_dot(v, n)));
}

/* Comparison */
static inline bool lg_vec3_eq(lg_vec3_t a, lg_vec3_t b, float eps) {
    return fabsf(a.x - b.x) < eps && fabsf(a.y - b.y) < eps && fabsf(a.z - b.z) < eps;
}

static inline bool lg_vec3_is_zero(lg_vec3_t v, float eps) {
    return fabsf(v.x) < eps && fabsf(v.y) < eps && fabsf(v.z) < eps;
}

/* Min/max */
static inline lg_vec3_t lg_vec3_min(lg_vec3_t a, lg_vec3_t b) {
    return lg_vec3(
        a.x < b.x ? a.x : b.x,
        a.y < b.y ? a.y : b.y,
        a.z < b.z ? a.z : b.z
    );
}

static inline lg_vec3_t lg_vec3_max(lg_vec3_t a, lg_vec3_t b) {
    return lg_vec3(
        a.x > b.x ? a.x : b.x,
        a.y > b.y ? a.y : b.y,
        a.z > b.z ? a.z : b.z
    );
}

static inline lg_vec3_t lg_vec3_clamp(lg_vec3_t v, lg_vec3_t min, lg_vec3_t max) {
    return lg_vec3(
        v.x < min.x ? min.x : (v.x > max.x ? max.x : v.x),
        v.y < min.y ? min.y : (v.y > max.y ? max.y : v.y),
        v.z < min.z ? min.z : (v.z > max.z ? max.z : v.z)
    );
}

/*============================================================================
 * 3D Vector (double) - for orbital mechanics precision
 *===========================================================================*/

typedef struct {
    double x, y, z;
} lg_vec3d_t;

static const lg_vec3d_t LG_VEC3D_ZERO = {0.0, 0.0, 0.0};

static inline lg_vec3d_t lg_vec3d(double x, double y, double z) {
    lg_vec3d_t v = {x, y, z};
    return v;
}

static inline lg_vec3d_t lg_vec3d_zero(void) { return LG_VEC3D_ZERO; }

static inline lg_vec3d_t lg_vec3d_add(lg_vec3d_t a, lg_vec3d_t b) {
    return lg_vec3d(a.x + b.x, a.y + b.y, a.z + b.z);
}

static inline lg_vec3d_t lg_vec3d_sub(lg_vec3d_t a, lg_vec3d_t b) {
    return lg_vec3d(a.x - b.x, a.y - b.y, a.z - b.z);
}

static inline lg_vec3d_t lg_vec3d_scale(lg_vec3d_t v, double s) {
    return lg_vec3d(v.x * s, v.y * s, v.z * s);
}

static inline double lg_vec3d_dot(lg_vec3d_t a, lg_vec3d_t b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

static inline lg_vec3d_t lg_vec3d_cross(lg_vec3d_t a, lg_vec3d_t b) {
    return lg_vec3d(
        a.y * b.z - a.z * b.y,
        a.z * b.x - a.x * b.z,
        a.x * b.y - a.y * b.x
    );
}

static inline double lg_vec3d_len_sq(lg_vec3d_t v) {
    return v.x * v.x + v.y * v.y + v.z * v.z;
}

static inline double lg_vec3d_len(lg_vec3d_t v) {
    return sqrt(lg_vec3d_len_sq(v));
}

static inline lg_vec3d_t lg_vec3d_norm(lg_vec3d_t v) {
    double len = lg_vec3d_len(v);
    if (len < 1e-12) return LG_VEC3D_ZERO;
    return lg_vec3d_scale(v, 1.0 / len);
}

/* Conversion between float and double */
static inline lg_vec3_t lg_vec3_from_d(lg_vec3d_t v) {
    return lg_vec3((float)v.x, (float)v.y, (float)v.z);
}

static inline lg_vec3d_t lg_vec3_to_d(lg_vec3_t v) {
    return lg_vec3d((double)v.x, (double)v.y, (double)v.z);
}

/*============================================================================
 * Quaternion
 *===========================================================================*/

typedef struct {
    float x, y, z, w;
} lg_quat_t;

static const lg_quat_t LG_QUAT_IDENTITY = {0.0f, 0.0f, 0.0f, 1.0f};

/* Construction */
static inline lg_quat_t lg_quat(float x, float y, float z, float w) {
    lg_quat_t q = {x, y, z, w};
    return q;
}

static inline lg_quat_t lg_quat_identity(void) { return LG_QUAT_IDENTITY; }

/* Basic operations */
static inline lg_quat_t lg_quat_add(lg_quat_t a, lg_quat_t b) {
    return lg_quat(a.x + b.x, a.y + b.y, a.z + b.z, a.w + b.w);
}

static inline lg_quat_t lg_quat_sub(lg_quat_t a, lg_quat_t b) {
    return lg_quat(a.x - b.x, a.y - b.y, a.z - b.z, a.w - b.w);
}

static inline lg_quat_t lg_quat_scale(lg_quat_t q, float s) {
    return lg_quat(q.x * s, q.y * s, q.z * s, q.w * s);
}

/* Quaternion multiplication (hamilton product) */
static inline lg_quat_t lg_quat_mul(lg_quat_t a, lg_quat_t b) {
    return lg_quat(
        a.w * b.x + a.x * b.w + a.y * b.z - a.z * b.y,
        a.w * b.y - a.x * b.z + a.y * b.w + a.z * b.x,
        a.w * b.z + a.x * b.y - a.y * b.x + a.z * b.w,
        a.w * b.w - a.x * b.x - a.y * b.y - a.z * b.z
    );
}

/* Conjugate and inverse */
static inline lg_quat_t lg_quat_conj(lg_quat_t q) {
    return lg_quat(-q.x, -q.y, -q.z, q.w);
}

static inline float lg_quat_len_sq(lg_quat_t q) {
    return q.x * q.x + q.y * q.y + q.z * q.z + q.w * q.w;
}

static inline float lg_quat_len(lg_quat_t q) {
    return sqrtf(lg_quat_len_sq(q));
}

static inline lg_quat_t lg_quat_inv(lg_quat_t q) {
    float len_sq = lg_quat_len_sq(q);
    if (len_sq < 1e-6f) return LG_QUAT_IDENTITY;
    return lg_quat_scale(lg_quat_conj(q), 1.0f / len_sq);
}

static inline lg_quat_t lg_quat_norm(lg_quat_t q) {
    float len = lg_quat_len(q);
    if (len < 1e-6f) return LG_QUAT_IDENTITY;
    return lg_quat_scale(q, 1.0f / len);
}

/* Rotate vector by quaternion */
static inline lg_vec3_t lg_quat_rotate(lg_quat_t q, lg_vec3_t v) {
    lg_quat_t p = lg_quat(v.x, v.y, v.z, 0.0f);
    lg_quat_t q_conj = lg_quat_conj(q);
    lg_quat_t result = lg_quat_mul(lg_quat_mul(q, p), q_conj);
    return lg_vec3(result.x, result.y, result.z);
}

/* Create quaternion from axis-angle */
static inline lg_quat_t lg_quat_from_axis_angle(lg_vec3_t axis, float angle) {
    float half_angle = angle * 0.5f;
    float s = sinf(half_angle);
    lg_vec3_t n = lg_vec3_norm(axis);
    return lg_quat_norm(lg_quat(n.x * s, n.y * s, n.z * s, cosf(half_angle)));
}

/* Create quaternion from two vectors (rotation from a to b) */
static inline lg_quat_t lg_quat_from_to(lg_vec3_t from, lg_vec3_t to) {
    lg_vec3_t a = lg_vec3_norm(from);
    lg_vec3_t b = lg_vec3_norm(to);
    float dot = lg_vec3_dot(a, b);
    
    if (dot > 0.999999f) {
        return LG_QUAT_IDENTITY;
    }
    
    if (dot < -0.999999f) {
        // 180 degree rotation, find orthogonal axis
        lg_vec3_t axis = lg_vec3_cross(a, lg_vec3(1.0f, 0.0f, 0.0f));
        if (lg_vec3_len_sq(axis) < 0.001f) {
            axis = lg_vec3_cross(a, lg_vec3(0.0f, 1.0f, 0.0f));
        }
        axis = lg_vec3_norm(axis);
        return lg_quat_from_axis_angle(axis, (float)M_PI);
    }
    
    lg_vec3_t cross = lg_vec3_cross(a, b);
    lg_quat_t q = lg_quat(cross.x, cross.y, cross.z, 1.0f + dot);
    return lg_quat_norm(q);
}

/* Spherical linear interpolation */
static inline lg_quat_t lg_quat_slerp(lg_quat_t a, lg_quat_t b, float t) {
    float dot = a.x * b.x + a.y * b.y + a.z * b.z + a.w * b.w;
    
    // If dot is negative, negate one quaternion to take shorter path
    if (dot < 0.0f) {
        b = lg_quat_scale(b, -1.0f);
        dot = -dot;
    }
    
    // If quaternions are very close, use lerp
    if (dot > 0.9995f) {
        lg_quat_t result = lg_quat_add(a, lg_quat_scale(lg_quat_sub(b, a), t));
        return lg_quat_norm(result);
    }
    
    float theta_0 = acosf(dot);
    float theta = theta_0 * t;
    float sin_theta = sinf(theta);
    float sin_theta_0 = sinf(theta_0);
    
    float s0 = cosf(theta) - dot * sin_theta / sin_theta_0;
    float s1 = sin_theta / sin_theta_0;
    
    return lg_quat_add(lg_quat_scale(a, s0), lg_quat_scale(b, s1));
}

/*============================================================================
 * 4x4 Matrix
 *===========================================================================*/

typedef struct {
    float m[16];
} lg_mat4_t;

static inline lg_mat4_t lg_mat4_identity(void) {
    lg_mat4_t m = {{
        1.0f, 0.0f, 0.0f, 0.0f,
        0.0f, 1.0f, 0.0f, 0.0f,
        0.0f, 0.0f, 1.0f, 0.0f,
        0.0f, 0.0f, 0.0f, 1.0f
    }};
    return m;
}

static inline lg_mat4_t lg_mat4_translation(lg_vec3_t t) {
    lg_mat4_t m = lg_mat4_identity();
    m.m[12] = t.x;
    m.m[13] = t.y;
    m.m[14] = t.z;
    return m;
}

static inline lg_mat4_t lg_mat4_scale(lg_vec3_t s) {
    lg_mat4_t m = lg_mat4_identity();
    m.m[0] = s.x;
    m.m[5] = s.y;
    m.m[10] = s.z;
    return m;
}

static inline lg_mat4_t lg_mat4_rotation(lg_quat_t q) {
    lg_quat_t n = lg_quat_norm(q);
    float xx = n.x * n.x;
    float yy = n.y * n.y;
    float zz = n.z * n.z;
    float xy = n.x * n.y;
    float xz = n.x * n.z;
    float yz = n.y * n.z;
    float wx = n.w * n.x;
    float wy = n.w * n.y;
    float wz = n.w * n.z;
    
    lg_mat4_t m = {{
        1.0f - 2.0f * (yy + zz), 2.0f * (xy - wz),     2.0f * (xz + wy),     0.0f,
        2.0f * (xy + wz),        1.0f - 2.0f * (xx + zz), 2.0f * (yz - wx),     0.0f,
        2.0f * (xz - wy),        2.0f * (yz + wx),     1.0f - 2.0f * (xx + yy), 0.0f,
        0.0f,                    0.0f,                 0.0f,                 1.0f
    }};
    return m;
}

static inline lg_vec3_t lg_mat4_mul_vec3(lg_mat4_t m, lg_vec3_t v) {
    lg_vec3_t result;
    result.x = m.m[0] * v.x + m.m[4] * v.y + m.m[8] * v.z + m.m[12];
    result.y = m.m[1] * v.x + m.m[5] * v.y + m.m[9] * v.z + m.m[13];
    result.z = m.m[2] * v.x + m.m[6] * v.y + m.m[10] * v.z + m.m[14];
    return result;
}

/*============================================================================
 * 3x3 Matrix (for attitude control, inertia tensors, etc.)
 *===========================================================================*/

typedef struct {
    float m[9];  /* Column-major: m[col*3 + row] */
} lg_mat3_t;

/* 3x3 matrix from column vectors */
static inline lg_mat3_t lg_mat3_from_cols(lg_vec3_t c0, lg_vec3_t c1, lg_vec3_t c2) {
    lg_mat3_t m;
    m.m[0] = c0.x; m.m[1] = c0.y; m.m[2] = c0.z;
    m.m[3] = c1.x; m.m[4] = c1.y; m.m[5] = c1.z;
    m.m[6] = c2.x; m.m[7] = c2.y; m.m[8] = c2.z;
    return m;
}

/* 3x3 identity matrix */
static inline lg_mat3_t lg_mat3_identity(void) {
    lg_mat3_t m = {{
        1.0f, 0.0f, 0.0f,
        0.0f, 1.0f, 0.0f,
        0.0f, 0.0f, 1.0f
    }};
    return m;
}

/* 3x3 matrix determinant */
static inline float lg_mat3_determinant(lg_mat3_t m) {
    return m.m[0] * (m.m[4] * m.m[8] - m.m[5] * m.m[7])
         - m.m[3] * (m.m[1] * m.m[8] - m.m[2] * m.m[7])
         + m.m[6] * (m.m[1] * m.m[5] - m.m[2] * m.m[4]);
}

/* 3x3 matrix inverse (returns zero matrix if singular) */
static inline lg_mat3_t lg_mat3_inverse(lg_mat3_t m) {
    float det = lg_mat3_determinant(m);
    if (fabsf(det) < 1e-10f) {
        return lg_mat3_identity();  /* Singular matrix - return identity as fallback */
    }
    
    float inv_det = 1.0f / det;
    lg_mat3_t inv;
    
    inv.m[0] = (m.m[4] * m.m[8] - m.m[5] * m.m[7]) * inv_det;
    inv.m[1] = (m.m[2] * m.m[7] - m.m[1] * m.m[8]) * inv_det;
    inv.m[2] = (m.m[1] * m.m[5] - m.m[2] * m.m[4]) * inv_det;
    inv.m[3] = (m.m[5] * m.m[6] - m.m[3] * m.m[8]) * inv_det;
    inv.m[4] = (m.m[0] * m.m[8] - m.m[2] * m.m[6]) * inv_det;
    inv.m[5] = (m.m[2] * m.m[3] - m.m[0] * m.m[5]) * inv_det;
    inv.m[6] = (m.m[3] * m.m[7] - m.m[4] * m.m[6]) * inv_det;
    inv.m[7] = (m.m[1] * m.m[6] - m.m[0] * m.m[7]) * inv_det;
    inv.m[8] = (m.m[0] * m.m[4] - m.m[1] * m.m[3]) * inv_det;
    
    return inv;
}

/* 3x3 matrix transpose */
static inline lg_mat3_t lg_mat3_transpose(lg_mat3_t m) {
    lg_mat3_t t;
    t.m[0] = m.m[0]; t.m[1] = m.m[3]; t.m[2] = m.m[6];
    t.m[3] = m.m[1]; t.m[4] = m.m[4]; t.m[5] = m.m[7];
    t.m[6] = m.m[2]; t.m[7] = m.m[5]; t.m[8] = m.m[8];
    return t;
}

/* 3x3 matrix-vector multiplication */
static inline lg_vec3_t lg_mat3_mul_vec3(lg_mat3_t m, lg_vec3_t v) {
    lg_vec3_t result;
    result.x = m.m[0] * v.x + m.m[3] * v.y + m.m[6] * v.z;
    result.y = m.m[1] * v.x + m.m[4] * v.y + m.m[7] * v.z;
    result.z = m.m[2] * v.x + m.m[5] * v.y + m.m[8] * v.z;
    return result;
}

/* 3x3 matrix multiplication */
static inline lg_mat3_t lg_mat3_mul(lg_mat3_t a, lg_mat3_t b) {
    lg_mat3_t result;
    for (int col = 0; col < 3; col++) {
        for (int row = 0; row < 3; row++) {
            result.m[col*3 + row] = 
                a.m[0*3 + row] * b.m[col*3 + 0] +
                a.m[1*3 + row] * b.m[col*3 + 1] +
                a.m[2*3 + row] * b.m[col*3 + 2];
        }
    }
    return result;
}

#ifdef __cplusplus
}
#endif

#endif /* LAGRANGE_MATH_H */
