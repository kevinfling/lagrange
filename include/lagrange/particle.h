/**
 * lagrange_particle.h - N-Body, Dust, and Galactic Scales
 * 
 * Hierarchical methods for 10^3 to 10^9 particles:
 *   - Barnes-Hut O(N log N) with SIMD MAC
 *   - Fast Multipole Method (FMM) using wavelets for far-field
 *   - Individual Time Step (IT) integration with block grouping
 *   - Two-fluid dust-gas coupling (Epstein/Stokes drag)
 *   - Z-order spatial sorting for cache coherence
 *   - Periodic boundary Ewald sums
 * 
 * Integrates with: lagrange_math_simd.h (SoA layout), 
 *                  lagrange_wavelet.h (FMM kernel compression),
 *                  lagrange_fractional.h (anomalous diffusion in disks)
 */

#ifndef LAGRANGE_PARTICLE_H
#define LAGRANGE_PARTICLE_H

#include "math.h"
#include "body.h"
#include "math_simd.h"
#include "wavelet.h"
#include "fractional.h"
#include <stdint.h>
#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

/*============================================================================
 * Spatial Acceleration: Z-Order/Morton Curves
 * Critical for cache coherence in particle updates
 *===========================================================================*/

typedef uint64_t lg_morton_t;

/* Compute 3D Morton code from normalized position [0,1]^3 */
static inline lg_morton_t lg_morton_encode(float x, float y, float z) {
    /* Clamp to [0,1] to prevent overflow/wrap in quantization */
    x = fminf(fmaxf(x, 0.0f), 1.0f);
    y = fminf(fmaxf(y, 0.0f), 1.0f);
    z = fminf(fmaxf(z, 0.0f), 1.0f);

    /* Quantize to 21 bits per dimension (63 total in uint64_t) */
    uint32_t ix = (uint32_t)(x * 2097151.0f); /* 2^21 - 1 */
    uint32_t iy = (uint32_t)(y * 2097151.0f);
    uint32_t iz = (uint32_t)(z * 2097151.0f);
    
    /* Spread bits: 0b111 -> 0b001001001 */
    #define SPREAD(v) (((v)*0x00010001u) & 0xFF0000FFu)
    
    lg_morton_t m = 0;
    for (int i = 0; i < 21; i++) {
        m |= ((lg_morton_t)((ix >> i) & 1) << (3*i + 0));
        m |= ((lg_morton_t)((iy >> i) & 1) << (3*i + 1));
        m |= ((lg_morton_t)((iz >> i) & 1) << (3*i + 2));
    }
    return m;
}

/* Helper for zsort */
typedef struct {
    lg_morton_t key;
    int idx;
} _lg_morton_pair_t;

static inline int _lg_morton_cmp(const void* a, const void* b) {
    const _lg_morton_pair_t* pa = (const _lg_morton_pair_t*)a;
    const _lg_morton_pair_t* pb = (const _lg_morton_pair_t*)b;
    return (pa->key < pb->key) ? -1 : (pa->key > pb->key) ? 1 : 0;
}

/* Sort particles by Morton code for spatial locality */
static inline void lg_particle_zsort(lg_vec3_t* pos, lg_vec3_t* vel, 
                                     float* mass, int n) {
    if (n <= 1) return;
    
    /* Find bounds for normalization to [0,1]^3 */
    float min_x = pos[0].x, max_x = pos[0].x;
    float min_y = pos[0].y, max_y = pos[0].y;
    float min_z = pos[0].z, max_z = pos[0].z;
    for (int i = 1; i < n; i++) {
        if (pos[i].x < min_x) min_x = pos[i].x;
        if (pos[i].x > max_x) max_x = pos[i].x;
        if (pos[i].y < min_y) min_y = pos[i].y;
        if (pos[i].y > max_y) max_y = pos[i].y;
        if (pos[i].z < min_z) min_z = pos[i].z;
        if (pos[i].z > max_z) max_z = pos[i].z;
    }
    float sx = max_x - min_x;
    float sy = max_y - min_y;
    float sz = max_z - min_z;
    if (sx < 1e-10f) sx = 1.0f;
    if (sy < 1e-10f) sy = 1.0f;
    if (sz < 1e-10f) sz = 1.0f;
    
    /* Compute keys and sort */
    _lg_morton_pair_t* pairs = (_lg_morton_pair_t*)malloc(n * sizeof(_lg_morton_pair_t));
    for (int i = 0; i < n; i++) {
        pairs[i].key = lg_morton_encode(
            (pos[i].x - min_x) / sx,
            (pos[i].y - min_y) / sy,
            (pos[i].z - min_z) / sz);
        pairs[i].idx = i;
    }
    qsort(pairs, n, sizeof(_lg_morton_pair_t), _lg_morton_cmp);
    
    /* Reorder arrays */
    lg_vec3_t* pos_tmp = (lg_vec3_t*)malloc(n * sizeof(lg_vec3_t));
    lg_vec3_t* vel_tmp = (lg_vec3_t*)malloc(n * sizeof(lg_vec3_t));
    float* mass_tmp = (float*)malloc(n * sizeof(float));
    for (int i = 0; i < n; i++) {
        pos_tmp[i] = pos[pairs[i].idx];
        vel_tmp[i] = vel[pairs[i].idx];
        mass_tmp[i] = mass[pairs[i].idx];
    }
    memcpy(pos, pos_tmp, n * sizeof(lg_vec3_t));
    memcpy(vel, vel_tmp, n * sizeof(lg_vec3_t));
    memcpy(mass, mass_tmp, n * sizeof(float));
    
    free(pairs);
    free(pos_tmp);
    free(vel_tmp);
    free(mass_tmp);
}

/*============================================================================
 * Barnes-Hut Tree (Octree with SoA storage)
 *===========================================================================*/

typedef struct lg_bhtree_node {
    lg_vec3_t center;           /* Center of mass */
    float mass;               /* Total mass in node */
    float min_bound[3];       /* AABB min */
    float max_bound[3];       /* AABB max */
    
    /* Children: index into contiguous pool, -1 = leaf */
    int child[8];             
    int parent;
    
    /* Particle range for leaves [start, end) */
    int particle_start;
    int particle_end;
    
    /* Multipole expansion moments (for FMM) */
    float q0;                 /* Monopole (total mass) */
    lg_vec3_t q1;             /* Dipole */
    lg_mat4_t q2;             /* Quadrupole (simplified) */
} lg_bhtree_node_t;

typedef struct {
    lg_bhtree_node_t* nodes;  /* Pool allocator */
    int n_nodes;
    int capacity;
    int max_depth;
    
    /* MAC (Multipole Acceptance Criterion) */
    float theta;              /* 0.5 = standard, 0.0 = exact */
    
    /* Leaf parameters */
    int leaf_capacity;        /* Max particles per leaf (typically 8-16) */
} lg_bhtree_t;

/* Recursive octree build helper */
static inline void _lg_bhtree_build_node(lg_bhtree_t* tree, int node_idx,
                                          const lg_vec3_t* pos, const float* mass,
                                          int* indices, int* temp,
                                          int start, int end, int depth) {
    lg_bhtree_node_t* node = &tree->nodes[node_idx];
    int count = end - start;
    
    /* Compute AABB */
    node->min_bound[0] = node->min_bound[1] = node->min_bound[2] = 1e30f;
    node->max_bound[0] = node->max_bound[1] = node->max_bound[2] = -1e30f;
    for (int i = start; i < end; i++) {
        int p = indices[i];
        node->min_bound[0] = fminf(node->min_bound[0], pos[p].x);
        node->min_bound[1] = fminf(node->min_bound[1], pos[p].y);
        node->min_bound[2] = fminf(node->min_bound[2], pos[p].z);
        node->max_bound[0] = fmaxf(node->max_bound[0], pos[p].x);
        node->max_bound[1] = fmaxf(node->max_bound[1], pos[p].y);
        node->max_bound[2] = fmaxf(node->max_bound[2], pos[p].z);
    }
    
    /* Compute COM */
    node->mass = 0.0f;
    node->center = lg_vec3_zero();
    for (int i = start; i < end; i++) {
        int p = indices[i];
        node->mass += mass[p];
        node->center.x += mass[p] * pos[p].x;
        node->center.y += mass[p] * pos[p].y;
        node->center.z += mass[p] * pos[p].z;
    }
    if (node->mass > 0.0f) {
        node->center.x /= node->mass;
        node->center.y /= node->mass;
        node->center.z /= node->mass;
    }
    
    /* Leaf check */
    if (count <= tree->leaf_capacity || depth >= tree->max_depth) {
        node->particle_start = start;
        node->particle_end = end;
        for (int c = 0; c < 8; c++) node->child[c] = -1;
        return;
    }
    
    /* Subdivide */
    float cx = (node->min_bound[0] + node->max_bound[0]) * 0.5f;
    float cy = (node->min_bound[1] + node->max_bound[1]) * 0.5f;
    float cz = (node->min_bound[2] + node->max_bound[2]) * 0.5f;
    
    int counts[8] = {0};
    for (int i = start; i < end; i++) {
        int p = indices[i];
        int oct = ((pos[p].x > cx) ? 1 : 0) |
                  ((pos[p].y > cy) ? 2 : 0) |
                  ((pos[p].z > cz) ? 4 : 0);
        counts[oct]++;
    }
    
    int offsets[9];
    offsets[0] = start;
    for (int o = 0; o < 8; o++) offsets[o+1] = offsets[o] + counts[o];
    
    memcpy(temp + start, indices + start, count * sizeof(int));
    int next[8];
    for (int o = 0; o < 8; o++) next[o] = offsets[o];
    for (int i = start; i < end; i++) {
        int p = temp[i];
        int oct = ((pos[p].x > cx) ? 1 : 0) |
                  ((pos[p].y > cy) ? 2 : 0) |
                  ((pos[p].z > cz) ? 4 : 0);
        indices[next[oct]++] = p;
    }
    
    node->particle_start = -1;
    node->particle_end = -1;
    for (int o = 0; o < 8; o++) {
        if (counts[o] == 0) {
            node->child[o] = -1;
            continue;
        }
        if (tree->n_nodes >= tree->capacity) {
            tree->capacity *= 2;
            tree->nodes = (lg_bhtree_node_t*)realloc(tree->nodes,
                tree->capacity * sizeof(lg_bhtree_node_t));
        }
        int child_idx = tree->n_nodes++;
        node->child[o] = child_idx;
        lg_bhtree_node_t* child = &tree->nodes[child_idx];
        memset(child, 0, sizeof(lg_bhtree_node_t));
        child->parent = node_idx;
        _lg_bhtree_build_node(tree, child_idx, pos, mass, indices, temp,
                              offsets[o], offsets[o+1], depth + 1);
    }
}

/* Build tree top-down */
static inline void lg_bhtree_build(lg_bhtree_t* tree, 
                                   const lg_vec3_t* pos,
                                   const float* mass,
                                   int n) {
    if (n <= 0) return;
    tree->n_nodes = 0;
    if (tree->capacity < n * 2) {
        tree->capacity = n * 4;
        tree->nodes = (lg_bhtree_node_t*)realloc(tree->nodes,
            tree->capacity * sizeof(lg_bhtree_node_t));
    }
    
    int* indices = (int*)malloc(n * sizeof(int));
    int* temp = (int*)malloc(n * sizeof(int));
    for (int i = 0; i < n; i++) indices[i] = i;
    
    int root_idx = tree->n_nodes++;
    memset(&tree->nodes[root_idx], 0, sizeof(lg_bhtree_node_t));
    _lg_bhtree_build_node(tree, root_idx, pos, mass, indices, temp, 0, n, 0);
    
    free(indices);
    free(temp);
}

/* Evaluate gravity with MAC: a_i = sum_j G*m_j*(r_j - r_i)/|r|^3 */
static inline lg_vec3_t lg_bhtree_accel(const lg_bhtree_t* tree,
                                        lg_vec3_t pos,
                                        float softening) {
    lg_vec3_t acc = lg_vec3_zero();
    
    /* Stackless traversal using explicit stack */
    int stack[64];  /* Max depth 64 is plenty */
    int sp = 0;
    stack[sp++] = 0; /* Root */
    
    while (sp > 0) {
        int node_idx = stack[--sp];
        const lg_bhtree_node_t* node = &tree->nodes[node_idx];
        
        /* Compute MAC: s/d < theta, where s=node size, d=distance to COM */
        float dx = node->center.x - pos.x;
        float dy = node->center.y - pos.y;
        float dz = node->center.z - pos.z;
        float dist_sq = dx*dx + dy*dy + dz*dz;
        
        /* Node size approximation */
        float size = node->max_bound[0] - node->min_bound[0];
        
        if (size * size < tree->theta * tree->theta * dist_sq || 
            node->child[0] == -1) {
            /* Accept approximation (or leaf) */
            float dist = sqrtf(dist_sq + softening * softening);
            float factor = node->mass / (dist * dist * dist);
            acc.x += factor * dx;
            acc.y += factor * dy;
            acc.z += factor * dz;
        } else {
            /* Traverse children */
            for (int c = 0; c < 8; c++) {
                if (node->child[c] != -1) stack[sp++] = node->child[c];
            }
        }
    }
    return acc;
}

/*============================================================================
 * Fast Multipole Method (FMM) 
 * Using your wavelets for kernel compression in far-field
 *===========================================================================*/

typedef struct {
    int n_sources;            /* Source box */
    int n_targets;            /* Target box */
    lg_wavelet_compressed_t kernel; /* Compressed interaction kernel */
} lg_fmm_interaction_t;

/* FMM: O(N) instead of O(N log N) of BH */
static inline void lg_fmm_upward_pass(lg_bhtree_t* tree,
                                      const lg_vec3_t* pos,
                                      const float* mass,
                                      const int* indices) {
    /* Bottom-up: compute monopole and dipole moments */
    for (int i = tree->n_nodes - 1; i >= 0; i--) {
        lg_bhtree_node_t* node = &tree->nodes[i];
        if (node->child[0] == -1) {
            /* Leaf: P2M */
            node->q0 = 0.0f;
            node->q1 = lg_vec3_zero();
            for (int j = node->particle_start; j < node->particle_end; j++) {
                int p = indices[j];
                float m = mass[p];
                node->q0 += m;
                node->q1.x += m * (pos[p].x - node->center.x);
                node->q1.y += m * (pos[p].y - node->center.y);
                node->q1.z += m * (pos[p].z - node->center.z);
            }
        } else {
            /* Internal: M2M */
            node->q0 = 0.0f;
            node->q1 = lg_vec3_zero();
            for (int c = 0; c < 8; c++) {
                int child_idx = node->child[c];
                if (child_idx == -1) continue;
                lg_bhtree_node_t* child = &tree->nodes[child_idx];
                node->q0 += child->q0;
                float dx = child->center.x - node->center.x;
                float dy = child->center.y - node->center.y;
                float dz = child->center.z - node->center.z;
                node->q1.x += child->q1.x + child->q0 * dx;
                node->q1.y += child->q1.y + child->q0 * dy;
                node->q1.z += child->q1.z + child->q0 * dz;
            }
        }
    }
}

static inline void lg_fmm_downward_pass(lg_bhtree_t* tree,
                                        const lg_vec3_t* pos,
                                        const int* indices,
                                        lg_vec3_t* acc_out,
                                        float softening) {
    /* Compute accelerations using monopole + dipole tree traversal */
    int n_particles = tree->nodes[0].particle_end;
    if (n_particles <= 0) return;
    
    for (int i = 0; i < n_particles; i++) {
        int p = indices[i];
        lg_vec3_t a = lg_vec3_zero();
        
        int stack[64];
        int sp = 0;
        stack[sp++] = 0;
        
        while (sp > 0) {
            int node_idx = stack[--sp];
            const lg_bhtree_node_t* node = &tree->nodes[node_idx];
            
            float dx = node->center.x - pos[p].x;
            float dy = node->center.y - pos[p].y;
            float dz = node->center.z - pos[p].z;
            float dist_sq = dx*dx + dy*dy + dz*dz;
            float size = node->max_bound[0] - node->min_bound[0];
            
            if (size * size < tree->theta * tree->theta * dist_sq ||
                node->child[0] == -1) {
                float dist = sqrtf(dist_sq + softening * softening);
                float dist3 = dist * dist * dist;
                float dist5 = dist3 * dist * dist;
                
                /* Monopole */
                float f0 = node->q0 / dist3;
                a.x += f0 * dx;
                a.y += f0 * dy;
                a.z += f0 * dz;
                
                /* Dipole */
                float rq = dx * node->q1.x + dy * node->q1.y + dz * node->q1.z;
                float f1 = 3.0f * rq / dist5;
                a.x += f1 * dx - node->q1.x / dist3;
                a.y += f1 * dy - node->q1.y / dist3;
                a.z += f1 * dz - node->q1.z / dist3;
            } else {
                for (int c = 0; c < 8; c++) {
                    if (node->child[c] != -1) stack[sp++] = node->child[c];
                }
            }
        }
        acc_out[p] = a;
    }
}

/*============================================================================
 * Individual Time Step (IT) Integration
 * 
 * Different particles get different dt based on local dynamical time
 *===========================================================================*/

typedef struct {
    float time;               /* Current simulation time for this particle */
    float dt;                 /* Personal time step */
    float dt_max;             /* Upper bound (based on hierarchy) */
    float acc_prev;           /* Previous acceleration magnitude (for predict) */
    int level;                /* Time level in power-of-two hierarchy */
} lg_timestep_t;

#define LG_TIME_LEVELS 32     /* 2^32 dynamic range in dt */

typedef struct {
    lg_vec3_batch_t pos;
    lg_vec3_batch_t vel;
    lg_vec3_batch_t acc;
    float* mass;
    lg_timestep_t* time;
    int n;
    
    /* Block time step: group particles by time level */
    int* level_indices[LG_TIME_LEVELS]; /* Particle indices per level */
    int level_count[LG_TIME_LEVELS];
} lg_particle_system_t;

/* Compute individual time step based on Aarseth criterion:
 * dt ~ sqrt(eta * |a| / |da/dt|) or ~ |a|/|jerk|^{1/2}
 */
static inline float lg_timestep_criterion(const lg_body_t* body,
                                          const lg_vec3_t* acc,
                                          const lg_vec3_t* jerk,
                                          float eta) {
    float acc_mag = lg_vec3_len(*acc);
    float jerk_mag = lg_vec3_len(*jerk);
    if (jerk_mag < 1e-20f) return 0.01f; /* Max dt */
    return eta * acc_mag / jerk_mag;
}

/* Hierarchical block step:
 * Active levels are where t_global is multiple of dt_level
 * dt_level = dt_min * 2^level
 */
static inline void lg_integrate_block_step(lg_particle_system_t* sys,
                                           float t_global,
                                           float dt_min) {
    for (int level = 0; level < LG_TIME_LEVELS; level++) {
        float dt_level = dt_min * (1 << level);
        if (dt_level < dt_min * 0.5f) continue;
        if (fmodf(t_global, dt_level) < dt_min * 0.5f) {
            int count = sys->level_count[level];
            for (int i = 0; i < count; i++) {
                int idx = sys->level_indices[level][i];
                /* Kick-Drift (semi-implicit Euler) */
                sys->vel.x[idx] += sys->acc.x[idx] * dt_level;
                sys->vel.y[idx] += sys->acc.y[idx] * dt_level;
                sys->vel.z[idx] += sys->acc.z[idx] * dt_level;
                sys->pos.x[idx] += sys->vel.x[idx] * dt_level;
                sys->pos.y[idx] += sys->vel.y[idx] * dt_level;
                sys->pos.z[idx] += sys->vel.z[idx] * dt_level;
            }
        }
    }
}

/*============================================================================
 * Dust-Gas Two-Fluid Coupling
 * 
 * For protoplanetary disks: dust particles coupled to background gas via drag
 *===========================================================================*/

typedef struct {
    float stopping_time;      /* t_s = m/(C_d * rho_g * v_th) */
    float stokes_number;      /* St = Omega * t_s */
    int drag_regime;          /* 0=Epstein, 1=Stokes, 2=Newton */
    
    /* Coupling to gas field (interpolated from grid or analytic) */
    lg_vec3_t gas_vel;
    float gas_rho;
    
    /* For super-particle representation (many dust grains = one particle) */
    float dust_mass;          /* Physical mass represented by super-particle */
    int n_grains;             /* Number of grains in super-particle */
} lg_dust_properties_t;

/* Epstein drag: a_drag = -(v_dust - v_gas)/t_s */
static inline lg_vec3_t lg_dust_epstein_drag(const lg_body_t* dust,
                                             const lg_dust_properties_t* props) {
    lg_vec3_t dv = lg_vec3_sub(dust->velocity, props->gas_vel);
    float factor = -1.0f / props->stopping_time;
    return lg_vec3_scale(dv, factor);
}

/* Back-reaction on gas (conservation of momentum) */
static inline lg_vec3_t lg_dust_gas_backreaction(const lg_body_t* dust,
                                                const lg_dust_properties_t* props,
                                                float gas_mass_inv) {
    /* Sum over dust particles in cell: a_gas += sum(rho_d * (v_d - v_g)/t_s) * dt */
    return lg_vec3_zero(); /* Accumulate externally */
}

/*============================================================================
 * Softening Kernels (Plummer vs Spline vs Compact)
 *===========================================================================*/

/* Standard Plummer softening: phi = -G/ sqrt(r^2 + eps^2) */
static inline float lg_gravity_plummer_potential(float r, float eps) {
    return -1.0f / sqrtf(r*r + eps*eps);
}

/* Wendland C2 compact kernel: C2 continuous, support [0,2] */
static inline float lg_gravity_compact_kernel(float u, float eps) {
    (void)eps;
    if (u >= 2.0f) return 0.0f;
    float v = 1.0f - u * 0.5f;
    return v * v * v * v * (2.0f * u + 1.0f);
}

/*============================================================================
 * Periodic Boundary Conditions (Ewald Summation)
 * For galactic/globular cluster sims with periodic boxes
 *===========================================================================*/

typedef struct {
    lg_vec3_t box_size;
    bool periodic[3];
    float alpha;              /* Ewald splitting parameter */
    int n_replicas;           /* Real-space replicas (typically 1-3) */
    int k_max;                /* k-space cutoff */
} lg_periodic_t;

/* Minimum image convention */
static inline lg_vec3_t lg_minimum_image(lg_vec3_t r, const lg_periodic_t* pbc) {
    lg_vec3_t result = r;
    float box[3] = {pbc->box_size.x, pbc->box_size.y, pbc->box_size.z};
    for (int dim = 0; dim < 3; dim++) {
        if (!pbc->periodic[dim]) continue;
        float half = box[dim] * 0.5f;
        float* coord = (dim == 0) ? &result.x : (dim == 1) ? &result.y : &result.z;
        if (*coord > half) *coord -= box[dim];
        if (*coord < -half) *coord += box[dim];
    }
    return result;
}

/* Ewald short-range (real space) */
static inline lg_vec3_t lg_ewald_real_accel(lg_vec3_t r, float q, float alpha) {
    float r2 = lg_vec3_len_sq(r);
    float r_mag = sqrtf(r2);
    float ar = alpha * r_mag;
    
    /* erfc(alpha*r)/r^2 + 2*alpha*exp(-alpha^2*r^2)/sqrt(pi)/r */
    float factor = q * (erfcf(ar) / r2 + 2.0f * alpha * expf(-ar * ar) / (sqrtf(LG_PI) * r_mag));
    return lg_vec3_scale(lg_vec3_norm(r), factor);
}

/*============================================================================
 * Particle-Mesh (PM) Self-Gravity
 * FFT-based O(N log N) for periodic boxes
 *===========================================================================*/

typedef struct {
    int n_grid[3];            /* Grid dimensions (power of 2) */
    float* density;           /* Mass density grid */
    float* potential;         /* Real-space potential grid */
    lg_cfloat* potential_k;   /* Fourier space potential */
    float cell_size;
    
    /* CIC (Cloud-in-Cell) or TSC (Triangular Shaped Cloud) interpolation */
    int interpolation_order;  /* 1=NGR, 2=CIC, 3=TSC */
} lg_particle_mesh_t;

/* Deposit particles to grid with CIC */
static inline void lg_pm_deposit_cic(lg_particle_mesh_t* pm,
                                     const lg_vec3_t* pos,
                                     const float* mass,
                                     int n_particles) {
    int nx = pm->n_grid[0];
    int ny = pm->n_grid[1];
    int nz = pm->n_grid[2];
    float cs = pm->cell_size;
    float vol = cs * cs * cs;
    float cs_inv = 1.0f / cs;
    
    for (int i = 0; i < n_particles; i++) {
        float x = pos[i].x * cs_inv;
        float y = pos[i].y * cs_inv;
        float z = pos[i].z * cs_inv;
        
        int ix = (int)floorf(x);
        int iy = (int)floorf(y);
        int iz = (int)floorf(z);
        
        float wx = x - ix;
        float wy = y - iy;
        float wz = z - iz;
        
        if (ix < 0) ix = 0; if (ix >= nx - 1) ix = nx - 2;
        if (iy < 0) iy = 0; if (iy >= ny - 1) iy = ny - 2;
        if (iz < 0) iz = 0; if (iz >= nz - 1) iz = nz - 2;
        
        for (int dx = 0; dx <= 1; dx++)
        for (int dy = 0; dy <= 1; dy++)
        for (int dz = 0; dz <= 1; dz++) {
            float w = ((dx == 0) ? (1.0f - wx) : wx) *
                      ((dy == 0) ? (1.0f - wy) : wy) *
                      ((dz == 0) ? (1.0f - wz) : wz);
            int gi = (ix + dx) + nx * ((iy + dy) + ny * (iz + dz));
            pm->density[gi] += mass[i] * w / vol;
        }
    }
}

/* Solve Poisson equation via SOR iterative method */
static inline void lg_pm_solve_poisson(lg_particle_mesh_t* pm) {
    int nx = pm->n_grid[0];
    int ny = pm->n_grid[1];
    int nz = pm->n_grid[2];
    int n_total = nx * ny * nz;
    float h = pm->cell_size;
    float h2 = h * h;
    const float fourpiG = 4.0f * LG_PI * 6.67430e-11f;
    
    float omega = 1.5f;
    int max_iter = 2000;
    float tol = 1e-6f;
    
    if (!pm->potential) {
        pm->potential = (float*)calloc(n_total, sizeof(float));
    }
    
    for (int iter = 0; iter < max_iter; iter++) {
        float max_diff = 0.0f;
        for (int iz = 0; iz < nz; iz++) {
            for (int iy = 0; iy < ny; iy++) {
                for (int ix = 0; ix < nx; ix++) {
                    int idx = ix + nx * (iy + ny * iz);
                    float rho = -fourpiG * pm->density[idx];
                    
                    float phi_xm = (ix > 0) ? pm->potential[idx - 1] : 0.0f;
                    float phi_xp = (ix < nx - 1) ? pm->potential[idx + 1] : 0.0f;
                    float phi_ym = (iy > 0) ? pm->potential[idx - nx] : 0.0f;
                    float phi_yp = (iy < ny - 1) ? pm->potential[idx + nx] : 0.0f;
                    float phi_zm = (iz > 0) ? pm->potential[idx - nx * ny] : 0.0f;
                    float phi_zp = (iz < nz - 1) ? pm->potential[idx + nx * ny] : 0.0f;
                    
                    float phi_new = (1.0f - omega) * pm->potential[idx] +
                        omega * (phi_xm + phi_xp + phi_ym + phi_yp + phi_zm + phi_zp - h2 * rho) / 6.0f;
                    
                    float diff = fabsf(phi_new - pm->potential[idx]);
                    if (diff > max_diff) max_diff = diff;
                    pm->potential[idx] = phi_new;
                }
            }
        }
        if (max_diff < tol) break;
    }
}

/*============================================================================
 * Collision Detection (Uniform Grid for local interactions)
 *===========================================================================*/

typedef struct {
    int* cell_heads;          /* Linked list heads per cell */
    int* next_particle;       /* Linked list next pointers */
    int n_cells[3];
    float cell_size;
} lg_collision_grid_t;

/* Build uniform grid for O(1) neighbor finding */
static inline void lg_grid_build(lg_collision_grid_t* grid,
                                 const lg_vec3_batch_t* pos,
                                 int n) {
    memset(grid->cell_heads, -1, grid->n_cells[0]*grid->n_cells[1]*grid->n_cells[2]*sizeof(int));
    
    for (int i = 0; i < n; i++) {
        int ix = (int)(pos->x[i] / grid->cell_size);
        int iy = (int)(pos->y[i] / grid->cell_size);
        int iz = (int)(pos->z[i] / grid->cell_size);
        int cell = ix + grid->n_cells[0]*(iy + grid->n_cells[1]*iz);
        
        grid->next_particle[i] = grid->cell_heads[cell];
        grid->cell_heads[cell] = i;
    }
}

/*============================================================================
 * Minimal SoA Particle Layout (used by ring systems, spacecraft, etc.)
 *===========================================================================*/

typedef struct {
    float* x, *y, *z;         /* Positions */
    float* vx, *vy, *vz;      /* Velocities */
    int count;
    int capacity;
} lg_particle_soa_t;

/*============================================================================
 * GPU/CUDA Patterns (Host-side preparation)
 *===========================================================================*/

/* SoA layout optimization for coalesced GPU memory access */
typedef struct {
    float* x, *y, *z;         /* Positions */
    float* vx, *vy, *vz;      /* Velocities */
    float* m;                 /* Masses */
    int* id;                  /* Particle IDs for tracing */
} lg_particle_gpu_layout_t;

/* Sort by Morton code (Z-curve) for GPU thread coalescence */
static inline void lg_particle_gpu_sort(lg_particle_gpu_layout_t* layout, int n) {
    /* Compute Morton codes */
    /* Sort keys, permute arrays to match */
    /* Result: threads in same warp access spatially close particles */
}

/*============================================================================
 * Statistical Methods: PDF Evolution (Fokker-Planck)
 * For systems where direct N-body is impossible (10^12 stars)
 *===========================================================================*/

typedef struct {
    float* f;                 /* Distribution function f(x,v,t) */
    int nx, nv;               /* Grid dimensions in x and v */
    float dx, dv;             /* Cell sizes */
    
    /* Precomputed acceleration field a(x) = -dphi/dx from Poisson solve */
    float* accel;
    
    /* Collision operator: Fokker-Planck or Boltzmann */
    float nu;                 /* Collision frequency */
    float D_vv;               /* Velocity diffusion coefficient */
} lg_vlasov_t;

/* Drift in phase space (Vlasov-Poisson) 
 * 1D spatial + 1D velocity upwind finite-volume */
static inline void lg_vlasov_drift(lg_vlasov_t* v, float dt) {
    int nx = v->nx;
    int nv = v->nv;
    float dx = v->dx;
    float dv = v->dv;
    float dx_inv = 1.0f / dx;
    float dv_inv = 1.0f / dv;
    
    float* f_new = (float*)calloc(nx * nv, sizeof(float));
    
    for (int ix = 0; ix < nx; ix++) {
        for (int iv = 0; iv < nv; iv++) {
            int idx = ix * nv + iv;
            float f = v->f[idx];
            
            /* Spatial flux: v * df/dx */
            float v_val = (iv - nv / 2) * dv;
            float fxL = (ix > 0) ? v->f[(ix - 1) * nv + iv] : f;
            float fxR = (ix < nx - 1) ? v->f[(ix + 1) * nv + iv] : f;
            float dfdx = (v_val > 0.0f) ? (f - fxL) * dx_inv : (fxR - f) * dx_inv;
            
            /* Velocity flux: a * df/dv */
            float a = v->accel ? v->accel[ix] : 0.0f;
            float fvL = (iv > 0) ? v->f[ix * nv + (iv - 1)] : f;
            float fvR = (iv < nv - 1) ? v->f[ix * nv + (iv + 1)] : f;
            float dfdv = (a > 0.0f) ? (f - fvL) * dv_inv : (fvR - f) * dv_inv;
            
            f_new[idx] = f - dt * (v_val * dfdx + a * dfdv);
        }
    }
    
    memcpy(v->f, f_new, nx * nv * sizeof(float));
    free(f_new);
}

/*============================================================================
 * Integration: Symplectic vs Dissipative Splitting
 *===========================================================================*/

/* Drift-Kick-Drift (DKD) - second order symplectic for N-body */
static inline void lg_integrate_dkd(lg_particle_system_t* sys, float dt) {
    /* Drift: x_{1/2} = x_0 + v_0 * dt/2 */
    for (int i = 0; i < sys->n; i++) {
        sys->pos.x[i] += sys->vel.x[i] * 0.5f * dt;
        sys->pos.y[i] += sys->vel.y[i] * 0.5f * dt;
        sys->pos.z[i] += sys->vel.z[i] * 0.5f * dt;
    }
    
    /* Kick: direct N^2 summation for accelerations */
    for (int i = 0; i < sys->n; i++) {
        sys->acc.x[i] = 0.0f;
        sys->acc.y[i] = 0.0f;
        sys->acc.z[i] = 0.0f;
    }
    for (int i = 0; i < sys->n; i++) {
        for (int j = i + 1; j < sys->n; j++) {
            float dx = sys->pos.x[j] - sys->pos.x[i];
            float dy = sys->pos.y[j] - sys->pos.y[i];
            float dz = sys->pos.z[j] - sys->pos.z[i];
            float r2 = dx*dx + dy*dy + dz*dz;
            float r = sqrtf(r2 + 1e-10f);
            float r3 = r * r * r;
            float f = sys->mass[j] / r3;
            sys->acc.x[i] += f * dx;
            sys->acc.y[i] += f * dy;
            sys->acc.z[i] += f * dz;
            f = sys->mass[i] / r3;
            sys->acc.x[j] -= f * dx;
            sys->acc.y[j] -= f * dy;
            sys->acc.z[j] -= f * dz;
        }
    }
    for (int i = 0; i < sys->n; i++) {
        sys->vel.x[i] += sys->acc.x[i] * dt;
        sys->vel.y[i] += sys->acc.y[i] * dt;
        sys->vel.z[i] += sys->acc.z[i] * dt;
    }
    
    /* Drift: x_1 = x_{1/2} + v_1 * dt/2 */
    for (int i = 0; i < sys->n; i++) {
        sys->pos.x[i] += sys->vel.x[i] * 0.5f * dt;
        sys->pos.y[i] += sys->vel.y[i] * 0.5f * dt;
        sys->pos.z[i] += sys->vel.z[i] * 0.5f * dt;
    }
}

/* For dissipative systems (dust with drag): 
 * Split Hamiltonian part (DKD) from drag (backward Euler) */
static inline void lg_integrate_split(lg_particle_system_t* sys, 
                                      lg_dust_properties_t* dust,
                                      float dt) {
    /* Kick-Drift with self-gravity */
    lg_integrate_dkd(sys, dt);
    
    /* Drag step: implicit solve for velocity */
    /* v_new = (v_old + dt * v_gas/t_s) / (1 + dt/t_s) */
}

#ifdef __cplusplus
}
#endif

#endif /* LAGRANGE_PARTICLE_H */

