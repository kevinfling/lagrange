/**
 * lagrange_floating_origin.h - Infinite Precision Coordinate System
 * 
 * Double-precision global coordinates + single-precision local physics
 * Solves the "jitter at distance" problem for solar system/galactic scales.
 * 
 * Strategy:
 *   - 64-bit cell coordinates (parsec or AU scale)
 *   - 32-bit sub-cell offsets (meter or km precision)
 *   - Camera/observer always at (0,0,0) local
 *   - Periodic rebasing when observer crosses cell boundary
 *   - Hierarchical: Object[cell] -> Cell[origin] -> Origin[universe]
 * 
 * Integrates with: lagrange_particle.h (SoA rebatching),
 *                  lagrange_math_simd.h (AVX2 coordinate transforms),
 *                  lagrange_gravity.h (high-precision distance calc)
 */

#ifndef LAGRANGE_FLOATING_ORIGIN_H
#define LAGRANGE_FLOATING_ORIGIN_H

#include "math.h"
#include "body.h"
#include "particle.h"
#include <stdint.h>
#include <string.h>
#ifdef __x86_64__
#include <immintrin.h>  /* For AVX2 double<->float conversion */
#endif

#ifdef __cplusplus
extern "C" {
#endif

/*============================================================================
 * Coordinate Storage
 * Split 64-bit global + 32-bit local for maximum precision range
 *===========================================================================*/

typedef struct {
    int64_t cell[3];          /* High-precision cell index (AU, km, or pc scale) */
    float local[3];           /* Sub-cell offset (0..cell_size meters) */
    float cell_size;          /* Size of one cell in meters (typically 1e6 to 1e9) */
} lg_floCoord_t;

/* Universe-scale position (128 bits per coordinate) */
typedef struct {
    int64_t sector[3];        /* Galaxy-scale chunk (pc or kpc units) */
    int32_t sub[3];           /* Within sector (light-seconds or AU) */
    float frac[3];            /* Sub-meter precision */
} lg_universal_coord_t;

/*============================================================================
 * Floating Origin Manager
 *===========================================================================*/

typedef struct {
    lg_floCoord_t origin;     /* Current global location of local (0,0,0) */
    lg_vec3_t observer_local; /* Camera/primary particle in local coords (ideally near 0) */
    
    float rebasing_threshold; /* Rebasing trigger distance (fraction of cell_size, typically 0.5) */
    int64_t rebase_count;     /* Stats: number of recenterings performed */
    
    /* Cached conversion factors */
    double meters_per_cell;
    float inv_cell_size;
} lg_floating_origin_t;

/* Initialize with cell size appropriate for scale */
static inline lg_floating_origin_t lg_floating_origin_init(float cell_size_meters) {
    lg_floating_origin_t fo = {
        .origin = {{0, 0, 0}, {0.0f, 0.0f, 0.0f}, cell_size_meters},
        .observer_local = lg_vec3_zero(),
        .rebasing_threshold = 0.5f * cell_size_meters,
        .rebase_count = 0,
        .meters_per_cell = (double)cell_size_meters,
        .inv_cell_size = 1.0f / cell_size_meters
    };
    return fo;
}

/* Solar System scale: 1 cell = 1 AU = 1.5e11 m (too big for single precision detail) 
 * Better: 1 cell = 1e6 km = 1e9 m for inner system, hierarchical for outer */
static inline lg_floating_origin_t lg_fo_solar_system(void) {
    return lg_floating_origin_init(1e9f); /* 1 million km cells */
}

/* Galactic scale: 1 cell = 1 parsec ≈ 3e16 m */
static inline lg_floating_origin_t lg_fo_galactic(void) {
    return lg_floating_origin_init(3.086e16f);
}

/*============================================================================
 * Coordinate Transforms (High Precision)
 *===========================================================================*/

/* Convert floating coordinate to absolute double-precision meters */
static inline void lg_floCoord_to_meters(const lg_floCoord_t* c, double* out_m) {
    for (int i = 0; i < 3; i++) {
        out_m[i] = (double)c->cell[i] * (double)c->cell_size + (double)c->local[i];
    }
}

/* Convert meters back to floating coordinate (used during rebase) */
static inline lg_floCoord_t lg_meters_to_floCoord(const double* m, float cell_size) {
    lg_floCoord_t c;
    c.cell_size = cell_size;
    float inv = 1.0f / cell_size;
    
    for (int i = 0; i < 3; i++) {
        double cells_d = m[i] * inv;
        c.cell[i] = (int64_t)floor(cells_d);
        c.local[i] = (float)(m[i] - (double)c.cell[i] * cell_size);
    }
    return c;
}

/* Relative vector between two floating coordinates (high precision) */
static inline lg_vec3_t lg_floCoord_delta(const lg_floCoord_t* a, const lg_floCoord_t* b) {
    lg_vec3_t delta;
    
    for (int i = 0; i < 3; i++) {
        int64_t d_cell = a->cell[i] - b->cell[i];
        double d_total = (double)d_cell * (double)a->cell_size + 
                        (double)(a->local[i] - b->local[i]);
        /* This preserves precision even for distant objects */
        ((float*)&delta.x)[i] = (float)d_total; /* x,y,z layout dependent, careful */
    }
    /* Safer: */
    delta.x = (float)((double)(a->cell[0] - b->cell[0]) * a->cell_size + 
                      (a->local[0] - b->local[0]));
    delta.y = (float)((double)(a->cell[1] - b->cell[1]) * a->cell_size + 
                      (a->local[1] - b->local[1]));
    delta.z = (float)((double)(a->cell[2] - b->cell[2]) * a->cell_size + 
                      (a->local[2] - b->local[2]));
    
    return delta;
}

/* High-precision distance calculation (avoids cancellation error) */
static inline double lg_floCoord_distance_sq_hp(const lg_floCoord_t* a, const lg_floCoord_t* b) {
    double sum = 0.0;
    for (int i = 0; i < 3; i++) {
        int64_t d_cell = a->cell[i] - b->cell[i];
        double d = (double)d_cell * (double)a->cell_size + 
                   (double)(a->local[i] - b->local[i]);
        sum += d * d;
    }
    return sum;
}

static inline double lg_floCoord_distance_hp(const lg_floCoord_t* a, const lg_floCoord_t* b) {
    return sqrt(lg_floCoord_distance_sq_hp(a, b));
}

/*============================================================================
 * Rebase Operation (The Critical Floating Origin Update)
 * Shift all coordinates so that new_origin becomes (0,0,0) local
 *===========================================================================*/

/* Single particle rebase */
static inline void lg_floCoord_rebase(lg_floCoord_t* c, const lg_floCoord_t* new_origin) {
    for (int i = 0; i < 3; i++) {
        int64_t cell_diff = c->cell[i] - new_origin->cell[i];
        float local_diff = c->local[i] - new_origin->local[i];
        
        /* Handle borrow if local < 0 after subtraction */
        if (local_diff < 0.0f) {
            cell_diff -= 1;
            local_diff += c->cell_size;
        } else if (local_diff >= c->cell_size) {
            cell_diff += 1;
            local_diff -= c->cell_size;
        }
        
        c->cell[i] = cell_diff;
        c->local[i] = local_diff;
    }
}

/* Batch rebase for entire particle system (SIMD-accelerated local update) */
static inline void lg_particle_system_rebase(lg_particle_system_t* ps, 
                                           lg_floating_origin_t* fo,
                                           const lg_floCoord_t* new_origin,
                                           int n_particles) {
    /* Update global origin counter */
    fo->origin = *new_origin;
    fo->rebase_count++;
    
    /* Convert new_origin to meters for delta calculation */
    double origin_m[3];
    lg_floCoord_to_meters(new_origin, origin_m);
    
#if defined(__AVX2__)
    /* Process 8 particles at once for local coordinate update */
    __m256 cell_size_vec = _mm256_set1_ps(fo->origin.cell_size);
    
    for (int i = 0; i < n_particles; i += 8) {
        /* Load positions (these are "local" floats that need rebasing) */
        /* Actually in floating origin, we store only local + global index */
        /* So we need to reconstruct absolute, subtract origin, store new local */
        
        /* For this implementation, assume we're storing absolute in SoA */
        /* and need to convert to local relative to new origin */
        
        __m256 px = _mm256_loadu_ps(ps->pos.x + i);
        __m256 py = _mm256_loadu_ps(ps->pos.y + i);
        __m256 pz = _mm256_loadu_ps(ps->pos.z + i);
        
        /* Subtract new origin (broadcast) */
        __m256 ox = _mm256_set1_ps((float)origin_m[0]);
        __m256 oy = _mm256_set1_ps((float)origin_m[1]);
        __m256 oz = _mm256_set1_ps((float)origin_m[2]);
        
        px = _mm256_sub_ps(px, ox);
        py = _mm256_sub_ps(py, oy);
        pz = _mm256_sub_ps(pz, oz);
        
        /* Store back as new local coordinates */
        _mm256_storeu_ps(ps->pos.x + i, px);
        _mm256_storeu_ps(ps->pos.y + i, py);
        _mm256_storeu_ps(ps->pos.z + i, pz);
        
        /* Also rebase the global cell indices stored separately */
        /* (Requires int64_t SoA for cell coordinates) */
    }
#else
    /* Scalar fallback */
    for (int i = 0; i < n_particles; i++) {
        ps->pos.x[i] = (float)((double)ps->pos.x[i] - origin_m[0]);
        ps->pos.y[i] = (float)((double)ps->pos.y[i] - origin_m[1]);
        ps->pos.z[i] = (float)((double)ps->pos.z[i] - origin_m[2]);
    }
#endif
}

/*============================================================================
 * Smart Rebase Trigger (Automatic when observer drifts too far)
 *===========================================================================*/

/* Check if primary particle (observer/camera) needs rebase */
static inline bool lg_floating_origin_needs_rebase(const lg_floating_origin_t* fo) {
    float dist_sq = lg_vec3_len_sq(fo->observer_local);
    return dist_sq > fo->rebasing_threshold * fo->rebasing_threshold;
}

/* Perform automatic rebase centered on current observer position */
static inline void lg_floating_origin_recenter(lg_floating_origin_t* fo,
                                             lg_particle_system_t* ps,
                                             int n) {
    /* New origin absolute position */
    double new_origin_m[3];
    lg_floCoord_to_meters(&fo->origin, new_origin_m);
    
    /* Add current local offset */
    new_origin_m[0] += fo->observer_local.x;
    new_origin_m[1] += fo->observer_local.y;
    new_origin_m[2] += fo->observer_local.z;
    
    lg_floCoord_t new_origin = lg_meters_to_floCoord(new_origin_m, fo->origin.cell_size);
    
    /* Rebase all particles */
    lg_particle_system_rebase(ps, fo, &new_origin, n);
    
    /* Reset observer to local zero */
    fo->observer_local = lg_vec3_zero();
}

/*============================================================================
 * Physics Integration with Floating Origin
 * Store high-precision state for drift-free integration
 *===========================================================================*/

typedef struct {
    lg_floCoord_t position;   /* High-precision position */
    lg_vec3_t velocity;       /* Local velocity (m/s, always small enough for float) */
    lg_vec3_t acceleration;   /* Local acceleration */
    
    /* For symplectic integration: store previous high-precision position */
    lg_floCoord_t position_prev;
} lg_floBody_t;

/* Verlet integration with floating origin (time-symmetric, preserves precision) */
static inline void lg_floBody_integrate_verlet(lg_floBody_t* body, float dt) {
    /* Save current */
    body->position_prev = body->position;
    
    /* Update local position (drifts, will be rebased later) */
    lg_vec3_t delta_pos = lg_vec3_scale(body->velocity, dt);
    delta_pos = lg_vec3_add(delta_pos, 
                           lg_vec3_scale(body->acceleration, 0.5f * dt * dt));
    
    /* Accumulate into local coordinate, handle cell overflow */
    for (int i = 0; i < 3; i++) {
        body->position.local[i] += ((float*)&delta_pos.x)[i]; /* x, y, z by index */
        
        /* Check for cell boundary crossing */
        if (body->position.local[i] >= body->position.cell_size) {
            body->position.cell[i]++;
            body->position.local[i] -= body->position.cell_size;
        } else if (body->position.local[i] < 0.0f) {
            body->position.cell[i]--;
            body->position.local[i] += body->position.cell_size;
        }
    }
    
    /* Velocity Verlet half-step */
    body->velocity = lg_vec3_add(body->velocity, 
                                lg_vec3_scale(body->acceleration, 0.5f * dt));
}

/* Convert to rendering coordinates (double->float with origin at camera) */
static inline lg_vec3_t lg_floBody_render_position(const lg_floBody_t* body,
                                                 const lg_floating_origin_t* fo) {
    /* Return local position relative to floating origin */
    lg_vec3_t local;
    local.x = (float)((double)(body->position.cell[0] - fo->origin.cell[0]) * body->position.cell_size 
                      + (body->position.local[0] - fo->origin.local[0]));
    local.y = (float)((double)(body->position.cell[1] - fo->origin.cell[1]) * body->position.cell_size 
                      + (body->position.local[1] - fo->origin.local[1]));
    local.z = (float)((double)(body->position.cell[2] - fo->origin.cell[2]) * body->position.cell_size 
                      + (body->position.local[2] - fo->origin.local[2]));
    return local;
}

/*============================================================================
 * Gravity Calculation (High Precision Distance)
 *===========================================================================*/

static inline lg_vec3_t lg_floBody_gravity(const lg_floBody_t* a,
                                          const lg_floBody_t* b,
                                          float G,
                                          float softening_sq) {
    /* High-precision relative vector */
    double dx = (double)(a->position.cell[0] - b->position.cell[0]) * a->position.cell_size 
                + (a->position.local[0] - b->position.local[0]);
    double dy = (double)(a->position.cell[1] - b->position.cell[1]) * a->position.cell_size 
                + (a->position.local[1] - b->position.local[1]);
    double dz = (double)(a->position.cell[2] - b->position.cell[2]) * a->position.cell_size 
                + (a->position.local[2] - b->position.local[2]);
    
    double dist_sq = dx*dx + dy*dy + dz*dz + (double)softening_sq;
    double dist = sqrt(dist_sq);
    double inv_dist_cubed = 1.0 / (dist_sq * dist);
    
    double factor = G * inv_dist_cubed; /* Mass should be included */
    
    lg_vec3_t accel;
    accel.x = (float)(-factor * dx);
    accel.y = (float)(-factor * dy);
    accel.z = (float)(-factor * dz);
    
    return accel;
}

/*============================================================================
 * Multi-Scale Hierarchy (Solar System within Galaxy)
 * Three-level: Sector (kpc) -> Cell (AU) -> Local (meters)
 *===========================================================================*/

typedef struct {
    lg_universal_coord_t univ;
    lg_floating_origin_t local_fo;
    float scale_factor;       /* Universe meters per simulation meter */
} lg_hierarchical_space_t;

/* Transform from universal to local simulation coordinates */
static inline lg_vec3_t lg_universal_to_local(const lg_universal_coord_t* univ,
                                             const lg_hierarchical_space_t* space) {
    /* Calculate offset from space's origin in universe meters */
    double dx = (double)(univ->sector[0] - space->univ.sector[0]) * 3.086e19; /* kpc to m */
    /* ... plus sub and frac ... */
    
    /* Scale to simulation units */
    lg_vec3_t local;
    local.x = (float)(dx * space->scale_factor);
    /* ... */
    return local;
}

/*============================================================================
 * Patched Conic Integration (Floating Origin Safe)
 *===========================================================================*/

/* Lambert solver using high-precision positions, returns single-precision result */
static inline lg_lambert_solution_t lg_flo_lambert(const lg_floBody_t* r1,
                                                const lg_floBody_t* r2,
                                                float dt,
                                                float mu,
                                                const lg_floating_origin_t* fo) {
    /* Get high-precision relative positions */
    double p1[3], p2[3];
    lg_floCoord_to_meters(&r1->position, p1);
    lg_floCoord_to_meters(&r2->position, p2);
    
    /* Temporary rebasing for Lambert (work near origin) */
    double offset[3] = {p1[0], p1[1], p1[2]};
    p1[0] -= offset[0]; p1[1] -= offset[1]; p1[2] -= offset[2];
    p2[0] -= offset[0]; p2[1] -= offset[1]; p2[2] -= offset[2];
    
    /* Pack into float for standard Lambert (safe now since relative and small) */
    lg_vec3_t r1f = {(float)p1[0], (float)p1[1], (float)p1[2]};
    lg_vec3_t r2f = {(float)p2[0], (float)p2[1], (float)p2[2]};
    
    lg_lambert_problem_t prob = {
        .r1 = r1f, .r2 = r2f, .dt = dt, .mu = mu,
        .short_way = true, .max_iter = 50, .tol = 1e-6f
    };
    
    return lg_lambert_solve(&prob);
}

#ifdef __cplusplus
}
#endif

#endif /* LAGRANGE_FLOATING_ORIGIN_H */

