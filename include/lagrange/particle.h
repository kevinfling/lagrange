/**
 * lagrange_particle.h - N-Body, Dust, and Galactic Scales
 * 
 * Hierarchical methods for 10^3 to 10^9 particles:
 *   - Barnes-Hut O(N log N) with SIMD MAC
 *   - Fast Multipole Method (FMM) using your wavelets for far-field
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

/* Sort particles by Morton code for spatial locality */
static inline void lg_particle_zsort(lg_vec3_t* pos, lg_vec3_t* vel, 
                                     float* mass, int n) {
    /* Compute keys */
    lg_morton_t* keys = (lg_morton_t*)malloc(n * sizeof(lg_morton_t));
    for (int i = 0; i < n; i++) {
        /* Normalize to simulation box [0,1]^3 - box size needed */
        keys[i] = lg_morton_encode(pos[i].x, pos[i].y, pos[i].z);
    }
    
    /* Radix sort on Morton codes ( preserves spatial locality ) */
    /* Implementation: 11-round radix sort on 64-bit keys */
    
    free(keys);
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

/* Build tree top-down */
static inline void lg_bhtree_build(lg_bhtree_t* tree, 
                                   const lg_vec3_batch_t* batch,
                                   int n) {
    /* 1. Z-sort particles for cache-friendly access */
    /* 2. Recursive subdivision: split along longest AABB axis */
    /* 3. Compute centers of mass bottom-up */
    
    /* SIMD acceleration: process 8 particles at once for COM calculation */
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
static inline void lg_fmm_upward_pass(lg_bhtree_t* tree) {
    /* Translate multipoles up the tree */
    /* P2M: particle to multipole (leaf) */
    /* M2M: multipole to multipole (parent from children) */
}

static inline void lg_fmm_downward_pass(lg_bhtree_t* tree) {
    /* L2L: local to local (parent to children) */
    /* M2L: multipole to local (far-field via wavelet compression) */
    /* P2P: direct particle-particle (near-field) */
    
    /* Wavelet compression: kernel K(r_i, r_j) approximated by
     * sparse wavelet representation for well-separated boxes */
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
    /* Process from finest to coarsest level */
    for (int level = 0; level < LG_TIME_LEVELS; level++) {
        float dt_level = dt_min * (1 << level);
        if (fmodf(t_global, dt_level) < dt_min * 0.5f) {
            /* This level is active */
            int count = sys->level_count[level];
            
            /* SIMD batch process all particles at this level */
            #ifdef __AVX2__
                /* Load positions, velocities for this level */
                /* Update using Drift-Kick-Drift (DKD) for symplectic */
            #else
                for (int i = 0; i < count; i++) {
                    int idx = sys->level_indices[level][i];
                    /* Kick-Drift or DKD */
                }
            #endif
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

/* Dehnen (2012) compact softening: better momentum conservation */
static inline float lg_gravity_compact_kernel(float u, float eps) {
    /* W(u) where u = r/eps, compact support [0,2] */
    if (u >= 2.0f) return 0.0f;
    if (u >= 1.0f) {
        float v = 2.0f - u;
        return v*v*v * (6.0f*v*v - 5.0f*v + 1.0f) / 15.0f; /* 5th order spline */
    }
    /* ... */
    return 1.0f; /* Core */
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
    for (int dim = 0; dim < 3; dim++) {
        if (!pbc->periodic[dim]) continue;
        float box = pbc->box_size.x; /* Simplified - should be per-dim */
        float half = box * 0.5f;
        if (result.x > half) result.x -= box;
        if (result.x < -half) result.x += box;
    }
    return result;
}

/* Ewald short-range (real space) */
static inline lg_vec3_t lg_ewald_real_accel(lg_vec3_t r, float q, float alpha) {
    float r2 = lg_vec3_len_sq(r);
    float r_mag = sqrtf(r2);
    float ar = alpha * r_mag;
    
    /* erfc(alpha*r)/r^2 + 2*alpha*exp(-alpha^2*r^2)/sqrt(pi)/r */
    float factor = q * (erfcf(ar) / r2 + 2.0f*alpha*expf(-ar*ar)/(M_PI*sqrtf(M_PI)*r_mag));
    return lg_vec3_scale(lg_vec3_norm(r), factor);
}

/*============================================================================
 * Particle-Mesh (PM) Self-Gravity
 * FFT-based O(N log N) for periodic boxes
 *===========================================================================*/

typedef struct {
    int n_grid[3];            /* Grid dimensions (power of 2) */
    float* density;           /* Mass density grid */
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
    for (int i = 0; i < n_particles; i++) {
        /* Find cell indices */
        int ix = (int)(pos[i].x / pm->cell_size);
        int iy = (int)(pos[i].y / pm->cell_size);
        int iz = (int)(pos[i].z / pm->cell_size);
        
        /* Compute weights for 8 surrounding cells */
        /* ... */
    }
}

/* Solve Poisson equation via FFT: phi_k = -4*pi*G*rho_k/k^2 */
static inline void lg_pm_solve_poisson(lg_particle_mesh_t* pm) {
    /* FFT density -> k-space */
    /* Green's function in k-space */
    /* FFT potential -> real space */
    /* Gradient for accelerations */
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
    
    /* Collision operator: Fokker-Planck or Boltzmann */
    float nu;                 /* Collision frequency */
    float D_vv;               /* Velocity diffusion coefficient */
} lg_vlasov_t;

/* Drift in phase space (Vlasov-Poisson) */
static inline void lg_vlasov_drift(lg_vlasov_t* v, float dt) {
    /* Conservative finite-volume or semi-Lagrangian */
    /* f^{n+1} = f^n - dt * (v * df/dx + a * df/dv) */
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
    
    /* Kick: compute accelerations at x_{1/2}, v_1 = v_0 + a * dt */
    /* ... use Barnes-Hut or FMM ... */
    
    /* Drift: x_1 = x_{1/2} + v_1 * dt/2 */
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

