/**
 * lagrange_math_simd.h - AVX2/NEON Acceleration for Lagrange Math
 * 
 * Drop this alongside math.h. It overrides hot paths with SIMD 
 * intrinsics while maintaining API compatibility.
 * 
 * For your Coffee Lake: AVX2 (256-bit), FMA3, 32-byte alignment.
 * For Jetson Orin Nano: NEON (128-bit), can process 4 floats.
 */

#ifndef LAGRANGE_MATH_SIMD_H
#define LAGRANGE_MATH_SIMD_H

#include "math.h"  /* Pull in scalar types */

#if defined(__x86_64__) && defined(__AVX2__)
    #include <immintrin.h>
    #define LG_SIMD_AVX2 1
    #define LG_SIMD_WIDTH 8  /* 8 floats */
#elif defined(__ARM_NEON) || defined(__aarch64__)
    #include <arm_neon.h>
    #define LG_SIMD_NEON 1
    #define LG_SIMD_WIDTH 4  /* 4 floats */
#endif

#ifdef __cplusplus
extern "C" {
#endif

/*============================================================================
 * SIMD-Aligned Storage (for batch constellation propagation)
 *===========================================================================*/

typedef struct {
    float* x;
    float* y;
    float* z;
    int count;      /* Must be multiple of LG_SIMD_WIDTH */
    int capacity;
} lg_vec3_batch_t;

/*============================================================================
 * Fast Approximations (Newton-Raphson refinement)
 *===========================================================================*/

/* Reciprocal square root: rsqrt(x) with 2 Newton iterations */
static inline float lg_rsqrt_fast(float x) {
    float y;
    
    #if defined(__AVX512VL__)
        /* AVX-512VL provides 14-bit accurate rsqrt */
        __m128 vx = _mm_set_ss(x);
        __m128 vr = _mm_rsqrt14_ss(vx, vx);
        y = _mm_cvtss_f32(vr);
    #elif defined(__SSE__)
        /* SSE provides ~12-bit accurate rsqrt */
        __m128 vx = _mm_set_ss(x);
        __m128 vr = _mm_rsqrt_ss(vx, vx);
        y = _mm_cvtss_f32(vr);
    #else
        /* Portable bit-level initial guess (magic number) */
        union { float f; int32_t i; } u = {x};
        u.i = 0x5f3759df - (u.i >> 1);  /* Quake III fast rsqrt */
        y = u.f;
    #endif
    
    /* Newton-Raphson: y = y * (1.5f - 0.5f * x * y * y) */
    float xh = 0.5f * x;
    y = y * (1.5f - xh * y * y);  /* Iteration 1 */
    y = y * (1.5f - xh * y * y);  /* Iteration 2 - enough for float precision */
    return y;
}

/* Fast normalize using rsqrt */
static inline lg_vec3_t lg_vec3_norm_fast(lg_vec3_t v) {
    float len_sq = lg_vec3_len_sq(v);
    float inv_len = lg_rsqrt_fast(len_sq);
    return lg_vec3_scale(v, inv_len);
}

/*============================================================================
 * AVX2 Vectorized Operations (8-way parallel)
 *===========================================================================*/

#ifdef LG_SIMD_AVX2

/* Dot product of 8 vectors simultaneously (SoA layout) */
static inline void lg_vec3_dot_batch_avx2(const float* ax, const float* ay, const float* az,
                                          const float* bx, const float* by, const float* bz,
                                          float* out, int n) {
    for (int i = 0; i < n; i += 8) {
        __m256 xa = _mm256_loadu_ps(ax + i);
        __m256 ya = _mm256_loadu_ps(ay + i);
        __m256 za = _mm256_loadu_ps(az + i);
        
        __m256 xb = _mm256_loadu_ps(bx + i);
        __m256 yb = _mm256_loadu_ps(by + i);
        __m256 zb = _mm256_loadu_ps(bz + i);
        
        /* FMA: x1*x2 + y1*y2 */
        __m256 dot = _mm256_mul_ps(xa, xb);
        dot = _mm256_fmadd_ps(ya, yb, dot);  /* FMA: dot + y1*y2 */
        dot = _mm256_fmadd_ps(za, zb, dot);  /* FMA: dot + z1*z2 */
        
        _mm256_storeu_ps(out + i, dot);
    }
}

/* Normalize 8 vectors at once using fast rsqrt */
static inline void lg_vec3_norm_batch_avx2(float* x, float* y, float* z, int n) {
    for (int i = 0; i < n; i += 8) {
        __m256 vx = _mm256_loadu_ps(x + i);
        __m256 vy = _mm256_loadu_ps(y + i);
        __m256 vz = _mm256_loadu_ps(z + i);
        
        /* len_sq = x^2 + y^2 + z^2 */
        __m256 len_sq = _mm256_mul_ps(vx, vx);
        len_sq = _mm256_fmadd_ps(vy, vy, len_sq);
        len_sq = _mm256_fmadd_ps(vz, vz, len_sq);
        
        /* rsqrt with 1 Newton iteration for vectorized case */
        __m256 rsqrt = _mm256_rsqrt_ps(len_sq);  /* Hardware estimate */
        /* Refine: rsqrt = rsqrt * (1.5 - 0.5 * len_sq * rsqrt * rsqrt) */
        __m256 half = _mm256_set1_ps(0.5f);
        __m256 three_half = _mm256_set1_ps(1.5f);
        __m256 rsqrt2 = _mm256_mul_ps(rsqrt, rsqrt);
        __m256 rsqrt3 = _mm256_mul_ps(rsqrt2, len_sq);
        __m256 rsqrt4 = _mm256_mul_ps(rsqrt3, half);
        __m256 rsqrt5 = _mm256_sub_ps(three_half, rsqrt4);
        rsqrt = _mm256_mul_ps(rsqrt, rsqrt5);
        
        _mm256_storeu_ps(x + i, _mm256_mul_ps(vx, rsqrt));
        _mm256_storeu_ps(y + i, _mm256_mul_ps(vy, rsqrt));
        _mm256_storeu_ps(z + i, _mm256_mul_ps(vz, rsqrt));
    }
}

/* Cross product: a × b for 8 pairs */
static inline void lg_vec3_cross_batch_avx2(const float* ax, const float* ay, const float* az,
                                            const float* bx, const float* by, const float* bz,
                                            float* cx, float* cy, float* cz, int n) {
    for (int i = 0; i < n; i += 8) {
        __m256 xa = _mm256_loadu_ps(ax + i);
        __m256 ya = _mm256_loadu_ps(ay + i);
        __m256 za = _mm256_loadu_ps(az + i);
        __m256 xb = _mm256_loadu_ps(bx + i);
        __m256 yb = _mm256_loadu_ps(by + i);
        __m256 zb = _mm256_loadu_ps(bz + i);
        
        /* c.x = a.y*b.z - a.z*b.y */
        __m256 cx = _mm256_fmsub_ps(ya, zb, _mm256_mul_ps(za, yb));
        /* c.y = a.z*b.x - a.x*b.z */
        __m256 cy = _mm256_fmsub_ps(za, xb, _mm256_mul_ps(xa, zb));
        /* c.z = a.x*b.y - a.y*b.x */
        __m256 cz = _mm256_fmsub_ps(xa, yb, _mm256_mul_ps(ya, xb));
        
        _mm256_storeu_ps(cx + i, cx);
        _mm256_storeu_ps(cy + i, cy);
        _mm256_storeu_ps(cz + i, cz);
    }
}

/* Quaternion multiplication (Hamilton product) - 4 quats at once */
static inline void lg_quat_mul_batch_avx2(const float* ax, const float* ay, const float* az, const float* aw,
                                          const float* bx, const float* by, const float* bz, const float* bw,
                                          float* cx, float* cy, float* cz, float* cw, int n) {
    for (int i = 0; i < n; i += 8) {
        __m256 xa = _mm256_loadu_ps(ax + i);
        __m256 ya = _mm256_loadu_ps(ay + i);
        __m256 za = _mm256_loadu_ps(az + i);
        __m256 wa = _mm256_loadu_ps(aw + i);
        
        __m256 xb = _mm256_loadu_ps(bx + i);
        __m256 yb = _mm256_loadu_ps(by + i);
        __m256 zb = _mm256_loadu_ps(bz + i);
        __m256 wb = _mm256_loadu_ps(bw + i);
        
        /* c.x = w1*x2 + x1*w2 + y1*z2 - z1*y2 */
        __m256 cx = _mm256_mul_ps(wa, xb);
        cx = _mm256_fmadd_ps(xa, wb, cx);
        cx = _mm256_fmadd_ps(ya, zb, cx);
        cx = _mm256_fnmadd_ps(za, yb, cx);  /* FMA with negation */
        
        /* c.y = w1*y2 - x1*z2 + y1*w2 + z1*x2 */
        __m256 cy = _mm256_mul_ps(wa, yb);
        cy = _mm256_fnmadd_ps(xa, zb, cy);
        cy = _mm256_fmadd_ps(ya, wb, cy);
        cy = _mm256_fmadd_ps(za, xb, cy);
        
        /* c.z = w1*z2 + x1*y2 - y1*x2 + z1*w2 */
        __m256 cz = _mm256_mul_ps(wa, zb);
        cz = _mm256_fmadd_ps(xa, yb, cz);
        cz = _mm256_fnmadd_ps(ya, xb, cz);
        cz = _mm256_fmadd_ps(za, wb, cz);
        
        /* c.w = w1*w2 - x1*x2 - y1*y2 - z1*z2 */
        __m256 cw = _mm256_mul_ps(wa, wb);
        cw = _mm256_fnmadd_ps(xa, xb, cw);
        cw = _mm256_fnmadd_ps(ya, yb, cw);
        cw = _mm256_fnmadd_ps(za, zb, cw);
        
        _mm256_storeu_ps(cx + i, cx);
        _mm256_storeu_ps(cy + i, cy);
        _mm256_storeu_ps(cz + i, cz);
        _mm256_storeu_ps(cw + i, cw);
    }
}

#endif /* LG_SIMD_AVX2 */

/*============================================================================
 * NEON Implementation (Jetson Orin Nano)
 *===========================================================================*/

#ifdef LG_SIMD_NEON

static inline void lg_vec3_dot_batch_neon(const float* ax, const float* ay, const float* az,
                                          const float* bx, const float* by, const float* bz,
                                          float* out, int n) {
    for (int i = 0; i < n; i += 4) {
        float32x4_t xa = vld1q_f32(ax + i);
        float32x4_t ya = vld1q_f32(ay + i);
        float32x4_t za = vld1q_f32(az + i);
        
        float32x4_t xb = vld1q_f32(bx + i);
        float32x4_t yb = vld1q_f32(by + i);
        float32x4_t zb = vld1q_f32(bz + i);
        
        float32x4_t dot = vmulq_f32(xa, xb);
        dot = vfmaq_f32(dot, ya, yb);  /* FMA: dot + y1*y2 */
        dot = vfmaq_f32(dot, za, zb);
        
        vst1q_f32(out + i, dot);
    }
}

/* Fast NEON rsqrt using VRSQRTE + 1 Newton iteration */
static inline void lg_vec3_norm_batch_neon(float* x, float* y, float* z, int n) {
    for (int i = 0; i < n; i += 4) {
        float32x4_t vx = vld1q_f32(x + i);
        float32x4_t vy = vld1q_f32(y + i);
        float32x4_t vz = vld1q_f32(z + i);
        
        float32x4_t len_sq = vmulq_f32(vx, vx);
        len_sq = vfmaq_f32(len_sq, vy, vy);
        len_sq = vfmaq_f32(len_sq, vz, vz);
        
        /* Hardware estimate */
        float32x4_t rsqrt = vrsqrteq_f32(len_sq);
        
        /* Newton: rsqrt = rsqrt * (3.0 - len_sq * rsqrt^2) * 0.5 */
        float32x4_t rsqrt2 = vmulq_f32(rsqrt, rsqrt);
        float32x4_t three = vmovq_n_f32(3.0f);
        float32x4_t half = vmovq_n_f32(0.5f);
        float32x4_t term = vfmsq_f32(three, len_sq, rsqrt2);  /* 3 - x*rsqrt^2 */
        rsqrt = vmulq_f32(vmulq_f32(rsqrt, term), half);
        
        vst1q_f32(x + i, vmulq_f32(vx, rsqrt));
        vst1q_f32(y + i, vmulq_f32(vy, rsqrt));
        vst1q_f32(z + i, vmulq_f32(vz, rsqrt));
    }
}

#endif /* LG_SIMD_NEON */

/*============================================================================
 * Generic Batch Dispatcher
 *===========================================================================*/

static inline void lg_vec3_batch_normalize(lg_vec3_batch_t* batch) {
    int n = batch->count;
    #if defined(LG_SIMD_AVX2)
        lg_vec3_norm_batch_avx2(batch->x, batch->y, batch->z, n);
    #elif defined(LG_SIMD_NEON)
        lg_vec3_norm_batch_neon(batch->x, batch->y, batch->z, n);
    #else
        /* Scalar fallback */
        for (int i = 0; i < n; i++) {
            float len_sq = batch->x[i]*batch->x[i] + batch->y[i]*batch->y[i] + batch->z[i]*batch->z[i];
            float inv_len = 1.0f / sqrtf(len_sq);
            batch->x[i] *= inv_len; batch->y[i] *= inv_len; batch->z[i] *= inv_len;
        }
    #endif
}

/*============================================================================
 * Constellation Gravity Integration (Batch Verlet)
 * 
 * Process 8 satellites simultaneously using SoA layout.
 * This is where you get 8x throughput on your Coffee Lake.
 *===========================================================================*/

typedef struct {
    lg_vec3_batch_t pos;    /* SoA positions */
    lg_vec3_batch_t vel;    /* SoA velocities */
    lg_vec3_batch_t acc;    /* SoA accelerations (for Verlet) */
    float* inv_mass;        /* Per-satellite mass */
    float* mu;              /* Gravitational parameters */
    int n;                  /* Must be multiple of 8 for AVX2 */
} lg_constellation_t;

/* Batch gravity acceleration: a_i = -mu * r_i / |r_i|^3 */
static inline void lg_constellation_gravity_batch(lg_constellation_t* cons) {
    int n = cons->n;
    
    #if defined(LG_SIMD_AVX2)
        __m256 zero = _mm256_setzero_ps();
        __m256 ones = _mm256_set1_ps(1.0f);
        
        for (int i = 0; i < n; i += 8) {
            __m256 px = _mm256_loadu_ps(cons->pos.x + i);
            __m256 py = _mm256_loadu_ps(cons->pos.y + i);
            __m256 pz = _mm256_loadu_ps(cons->pos.z + i);
            
            /* r^2 = x^2 + y^2 + z^2 */
            __m256 r2 = _mm256_mul_ps(px, px);
            r2 = _mm256_fmadd_ps(py, py, r2);
            r2 = _mm256_fmadd_ps(pz, pz, r2);
            
            /* r^3 for gravity law */
            __m256 r = _mm256_sqrt_ps(r2);
            __m256 r3 = _mm256_mul_ps(r, r2);  /* r^3 */
            
            /* Load mu values */
            __m256 mu = _mm256_loadu_ps(cons->mu + i);
            
            /* factor = -mu / r^3 */
            __m256 factor = _mm256_div_ps(mu, r3);
            factor = _mm256_xor_ps(factor, _mm256_set1_ps(-0.0f)); /* Negate */
            
            _mm256_storeu_ps(cons->acc.x + i, _mm256_mul_ps(px, factor));
            _mm256_storeu_ps(cons->acc.y + i, _mm256_mul_ps(py, factor));
            _mm256_storeu_ps(cons->acc.z + i, _mm256_mul_ps(pz, factor));
        }
    #else
        /* Scalar fallback */
        for (int i = 0; i < n; i++) {
            float px = cons->pos.x[i], py = cons->pos.y[i], pz = cons->pos.z[i];
            float r2 = px*px + py*py + pz*pz;
            float r = sqrtf(r2);
            float factor = -cons->mu[i] / (r * r2);
            cons->acc.x[i] = px * factor;
            cons->acc.y[i] = py * factor;
            cons->acc.z[i] = pz * factor;
        }
    #endif
}

#ifdef __cplusplus
}
#endif

#endif /* LAGRANGE_MATH_SIMD_H */

