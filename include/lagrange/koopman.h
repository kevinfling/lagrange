/**
 * lagrange_koopman.h - Koopman Operator Theory for Orbital Mechanics
 * 
 * EDMD (Extended Dynamic Mode Decomposition) with spectral observables.
 * Lifts nonlinear Hamiltonian dynamics to linear operator in DCT basis.
 * 
 * Integrates with: CHEAP (DCT basis), lagrange_math_simd.h (batch modes)
 * Applications:
 *   - Fast trajectory prediction (O(N) vs O(N log N) integration)
 *   - Stability analysis (eigenvalues of Koopman operator)
 *   - Koopman-MPC for station-keeping (convex optimization in lifted space)
 *   - Fractional Koopman (memory effects in spectral domain)
 */

#ifndef LAGRANGE_KOOPMAN_H
#define LAGRANGE_KOOPMAN_H

#include "math.h"
#include "body.h"
#include "integrator.h"
#include "fractional.h"
#include "math_simd.h"
#include <complex.h>
#include <stdlib.h>
#include <string.h>

#if defined(_WIN32)
#include <malloc.h>
#else
#include <alloca.h>
#endif

#ifdef __cplusplus
extern "C" {
#endif

/*============================================================================
 * Koopman Observable Basis
 * 
 * Standard EDMD uses RBFs or polynomials. We use DCT-II basis functions
 * because they diagonalize Toeplitz covariances (CHEAP) and have fast transforms.
 *===========================================================================*/

typedef enum {
    LG_KOOP_DCT,        /* Chebyshev/DCT spectral basis (fast, orbital mechanics) */
    LG_KOOP_RBF,        /* Radial basis functions (generic nonlinear) */
    LG_KOOP_HERMITE,    /* Probabilist's Hermite (energy eigenfunctions) */
    LG_KOOP_FRACTIONAL  /* Fractional power laws (memory/anomalous diffusion) */
} lg_koop_basis_t;

typedef struct {
    lg_koop_basis_t type;
    int n_observables;  /* Dimension of lifted space (D) */
    int n_states;       /* Original state dim (typically 6 for orbital) */
    
    /* DCT-specific: frequency bins for spectral lifting */
    float* frequencies; /* [D] array of frequencies omega_k */
    float* scales;      /* Scaling for each observable */
    
    /* EDMD matrices (allocated on fit) */
    float complex* K;   /* Koopman operator matrix [D x D] */
    float complex* P;   /* Projection matrix [n_states x D] for state recovery */
    float complex* Psi; /* Observable matrix at training points [D x N] */
    
    /* For kernel Koopman (alternative to explicit basis) */
    float* K_kernel;    /* Kernel matrix [N x N] if using kernel EDMD */
    
    /* Fractional extension */
    float hurst;        /* H parameter for fractional Koopman */
    float* memory_kernel; /* Fractional power law weights */
} lg_koopman_t;

/*============================================================================
 * Observable Functions (Lift: R^n -> C^D)
 *===========================================================================*/

/* Hybrid DCT/Linear Observable:
 * First 6 modes are direct state components (linear observables).
 * Higher modes are complex exponentials of state components. */
static inline float complex lg_koop_observable_dct(const lg_vec3_t* pos, 
                                                  const lg_vec3_t* vel,
                                                  float omega, 
                                                  int mode) {
    /* Linear modes: direct state components */
    if (mode < 6) {
        switch (mode) {
            case 0: return pos->x;
            case 1: return pos->y;
            case 2: return pos->z;
            case 3: return vel->x;
            case 4: return vel->y;
            default:
            case 5: return vel->z;
        }
    }
    
    /* Spectral modes: complex exponentials */
    float phase;
    switch (mode % 6) {
        case 0: phase = omega * pos->x; break;
        case 1: phase = omega * pos->y; break;
        case 2: phase = omega * pos->z; break;
        case 3: phase = omega * vel->x; break;
        case 4: phase = omega * vel->y; break;
        default:
        case 5: phase = omega * vel->z; break;
    }
    return cexpf(I * phase);
}

/* Polynomial observable (Hermite-like) for energy moments */
static inline float complex lg_koop_observable_hermite(const lg_vec3_t* pos,
                                                        const lg_vec3_t* vel,
                                                        int order) {
    float r2 = lg_vec3_len_sq(*pos);
    float v2 = lg_vec3_len_sq(*vel);
    
    /* H_0 = 1, H_1 = x, H_2 = x^2 - 1, etc. Simplified here */
    switch(order) {
        case 0: return 1.0f;
        case 1: return r2 + v2;  /* Energy-like */
        case 2: return (r2 + v2) * (r2 + v2) - 2.0f;
        default: return cpowf(r2 + v2, (float)order * 0.5f);
    }
}

/* Fractional/memory observable: incorporates history via power law */
static inline float complex lg_koop_observable_fractional(lg_fractional_t* frac,
                                                          const lg_vec3_t* pos,
                                                          float omega) {
    /* Combine current state with fractional integral of history */
    lg_vec3_t memory = lg_fractional_integral_direct(frac);
    float phase = omega * (pos->x + memory.x);
    return cexpf(I * phase);
}

/*============================================================================
 * EDMD Training (Data-Driven Koopman Approximation)
 *===========================================================================*/

typedef struct {
    float* X;       /* State snapshots [n_states x n_snapshots] */
    float* Y;       /* Evolved states [n_states x n_snapshots] */
    int n_snapshots;
    int n_states;
} lg_koop_data_t;

/* Solve A * X = B for X using Gaussian elimination with partial pivoting.
 * A is n x n (column-major), B is n x nrhs (column-major). Both overwritten. */
static inline bool lg_koopman_solve_complex(int n, float complex* A, int nrhs, float complex* B) {
    int* ipiv = (int*)alloca((size_t)n * sizeof(int));
    
    for (int k = 0; k < n; k++) {
        int max_row = k;
        float max_val = cabsf(A[k + k*n]);
        for (int i = k + 1; i < n; i++) {
            float val = cabsf(A[i + k*n]);
            if (val > max_val) {
                max_val = val;
                max_row = i;
            }
        }
        ipiv[k] = max_row;
        if (max_val < 1e-15f) return false;
        
        if (max_row != k) {
            for (int j = 0; j < n; j++) {
                float complex tmp = A[j + k*n];
                A[j + k*n] = A[j + max_row*n];
                A[j + max_row*n] = tmp;
            }
            for (int j = 0; j < nrhs; j++) {
                float complex tmp = B[k + j*n];
                B[k + j*n] = B[max_row + j*n];
                B[max_row + j*n] = tmp;
            }
        }
        
        for (int i = k + 1; i < n; i++) {
            float complex lik = A[i + k*n] / A[k + k*n];
            A[i + k*n] = lik;
            for (int j = k + 1; j < n; j++) {
                A[i + j*n] -= lik * A[k + j*n];
            }
            for (int j = 0; j < nrhs; j++) {
                B[i + j*n] -= lik * B[k + j*n];
            }
        }
    }
    
    for (int k = n - 1; k >= 0; k--) {
        for (int j = 0; j < nrhs; j++) {
            float complex sum = B[k + j*n];
            for (int i = k + 1; i < n; i++) {
                sum -= A[k + i*n] * B[i + j*n];
            }
            B[k + j*n] = sum / A[k + k*n];
        }
    }
    
    return true;
}

/* Core EDMD matrix algebra: given Psi_X [D x N] and Psi_Y [D x N],
 * compute Koopman operator K [D x D] and projection P [S x D]. */
static inline void _lg_koopman_fit_core(lg_koopman_t* km, const lg_koop_data_t* data,
                                        float complex* Psi_X, float complex* Psi_Y) {
    int D = km->n_observables;
    int N = data->n_snapshots;
    int S = km->n_states;
    
    /* Gram matrix G = Psi_X * Psi_X^H */
    float complex* G = (float complex*)calloc((size_t)D * (size_t)D, sizeof(float complex));
    for (int i = 0; i < D; i++) {
        for (int j = 0; j < D; j++) {
            for (int n = 0; n < N; n++) {
                G[i + j*D] += Psi_X[i + n*D] * conjf(Psi_X[j + n*D]);
            }
        }
    }
    
    /* Ridge regularization */
    for (int i = 0; i < D; i++) G[i + i*D] += 1e-6f;
    
    /* A = Psi_Y * Psi_X^H */
    float complex* A = (float complex*)calloc((size_t)D * (size_t)D, sizeof(float complex));
    for (int i = 0; i < D; i++) {
        for (int j = 0; j < D; j++) {
            for (int n = 0; n < N; n++) {
                A[i + j*D] += Psi_Y[i + n*D] * conjf(Psi_X[j + n*D]);
            }
        }
    }
    
    /* Solve G * K^H = A^H for K */
    float complex* X = (float complex*)calloc((size_t)D * (size_t)D, sizeof(float complex));
    for (int i = 0; i < D; i++) {
        for (int j = 0; j < D; j++) {
            X[j + i*D] = conjf(A[i + j*D]);
        }
    }
    
    float complex* G_copy = (float complex*)calloc((size_t)D * (size_t)D, sizeof(float complex));
    memcpy(G_copy, G, D * D * sizeof(float complex));
    lg_koopman_solve_complex(D, G_copy, D, X);
    
    free(km->K);
    km->K = (float complex*)calloc((size_t)D * (size_t)D, sizeof(float complex));
    for (int i = 0; i < D; i++) {
        for (int j = 0; j < D; j++) {
            km->K[i + j*D] = conjf(X[j + i*D]);
        }
    }
    
    /* Projection matrix P: maps observables back to state space via least squares */
    float complex* Bmat = (float complex*)calloc((size_t)S * (size_t)D, sizeof(float complex));
    for (int i = 0; i < S; i++) {
        for (int j = 0; j < D; j++) {
            for (int n = 0; n < N; n++) {
                Bmat[i + j*S] += data->X[i + n*S] * conjf(Psi_X[j + n*D]);
            }
        }
    }
    
    float complex* X2 = (float complex*)calloc((size_t)D * (size_t)S, sizeof(float complex));
    for (int i = 0; i < S; i++) {
        for (int j = 0; j < D; j++) {
            X2[j + i*D] = conjf(Bmat[i + j*S]);
        }
    }
    
    memcpy(G_copy, G, D * D * sizeof(float complex));
    lg_koopman_solve_complex(D, G_copy, S, X2);
    
    free(km->P);
    km->P = (float complex*)calloc((size_t)S * (size_t)D, sizeof(float complex));
    for (int i = 0; i < S; i++) {
        for (int j = 0; j < D; j++) {
            km->P[i + j*S] = conjf(X2[j + i*D]);
        }
    }
    
    free(G);
    free(A);
    free(X);
    free(G_copy);
    free(Bmat);
    free(X2);
}

/* Fit Koopman operator via EDMD: K = Psi_Y * Psi_X^+ (pseudoinverse) */
static inline void lg_koopman_fit(lg_koopman_t* km, const lg_koop_data_t* data) {
    int D = km->n_observables;
    int N = data->n_snapshots;
    int S = km->n_states;
    
    /* Allocate observable matrices */
    free(km->Psi);
    km->Psi = (float complex*)calloc((size_t)D * (size_t)N, sizeof(float complex));
    float complex* Psi_Y = (float complex*)calloc((size_t)D * (size_t)N, sizeof(float complex));
    
    /* Compute observables at X and Y */
    for (int n = 0; n < N; n++) {
        lg_vec3_t pos = {data->X[0 + n*S], data->X[1 + n*S], data->X[2 + n*S]};
        lg_vec3_t vel = {data->X[3 + n*S], data->X[4 + n*S], data->X[5 + n*S]};
        
        lg_vec3_t pos_y = {data->Y[0 + n*S], data->Y[1 + n*S], data->Y[2 + n*S]};
        lg_vec3_t vel_y = {data->Y[3 + n*S], data->Y[4 + n*S], data->Y[5 + n*S]};
        
        for (int k = 0; k < D; k++) {
            km->Psi[k + n*D] = lg_koop_observable_dct(&pos, &vel, km->frequencies[k], k);
            Psi_Y[k + n*D] = lg_koop_observable_dct(&pos_y, &vel_y, km->frequencies[k], k);
        }
    }
    
    _lg_koopman_fit_core(km, data, km->Psi, Psi_Y);
    free(Psi_Y);
}

/*============================================================================
 * Koopman Prediction (Fast Linear Propagation)
 * 
 * Instead of integrating nonlinear ODEs, apply linear operator K:
 * psi(x_{t+dt}) = K * psi(x_t)
 * x_{t+dt} = argmin ||psi(x) - K*psi(x_t)|| (or linear projection)
 *===========================================================================*/

/* Advance state using Koopman operator (spectral method) */
static inline void lg_koopman_predict(const lg_koopman_t* km,
                                      const lg_vec3_t* pos_in,
                                      const lg_vec3_t* vel_in,
                                      lg_vec3_t* pos_out,
                                      lg_vec3_t* vel_out,
                                      int n_steps) {
    int D = km->n_observables;
    int S = km->n_states;
    
    /* Lift current state to observable space */
    float complex* psi = (float complex*)alloca((size_t)D * sizeof(float complex));
    for (int k = 0; k < D; k++) {
        psi[k] = lg_koop_observable_dct(pos_in, vel_in, km->frequencies[k], k);
    }
    
    /* Apply K^n_steps via repeated matrix multiply */
    float complex* psi_next = (float complex*)alloca((size_t)D * sizeof(float complex));
    
    for (int step = 0; step < n_steps; step++) {
        for (int i = 0; i < D; i++) {
            psi_next[i] = 0;
            for (int j = 0; j < D; j++) {
                psi_next[i] += km->K[i + j*D] * psi[j];
            }
        }
        memcpy(psi, psi_next, D * sizeof(float complex));
    }
    
    /* Project back to state space using learned projection matrix P */
    float complex px = 0, py = 0, pz = 0;
    float complex vx = 0, vy = 0, vz = 0;
    for (int k = 0; k < D; k++) {
        px += km->P[0 + k*S] * psi[k];
        py += km->P[1 + k*S] * psi[k];
        pz += km->P[2 + k*S] * psi[k];
        vx += km->P[3 + k*S] * psi[k];
        vy += km->P[4 + k*S] * psi[k];
        vz += km->P[5 + k*S] * psi[k];
    }
    
    pos_out->x = crealf(px);
    pos_out->y = crealf(py);
    pos_out->z = crealf(pz);
    vel_out->x = crealf(vx);
    vel_out->y = crealf(vy);
    vel_out->z = crealf(vz);
}

/*============================================================================
 * Koopman-MPC (Model Predictive Control in Lifted Space)
 * 
 * Convex optimization in observable space, nonlinear in original.
 * Perfect for station-keeping with collision avoidance.
 *===========================================================================*/

typedef struct {
    float complex* target_psi;  /* Target in lifted space */
    float control_cost;         /* R weight */
    float state_cost;           /* Q weight */
    int horizon;                /* MPC horizon */
} lg_koop_mpc_t;

/* Compute optimal control input via Koopman linearization.
 * In lifted space, the MPC problem reduces to linear feedback:
 *   u = -gain * P_vel * (K^horizon * psi_current - psi_target)
 * where P_vel projects the observable error back to velocity space. */
static inline lg_vec3_t lg_koopman_mpc_control(const lg_koopman_t* km,
                                               const lg_koop_mpc_t* mpc,
                                               const lg_vec3_t* current_pos,
                                               const lg_vec3_t* current_vel) {
    int D = km->n_observables;
    int S = km->n_states;
    
    /* Lift current state */
    float complex* psi = (float complex*)alloca((size_t)D * sizeof(float complex));
    for (int k = 0; k < D; k++) {
        psi[k] = lg_koop_observable_dct(current_pos, current_vel, km->frequencies[k], k);
    }
    
    /* Propagate forward over the MPC horizon: psi_pred = K^horizon * psi */
    float complex* psi_pred = (float complex*)alloca((size_t)D * sizeof(float complex));
    memcpy(psi_pred, psi, D * sizeof(float complex));
    for (int step = 0; step < mpc->horizon; step++) {
        float complex* tmp = (float complex*)alloca((size_t)D * sizeof(float complex));
        for (int i = 0; i < D; i++) {
            tmp[i] = 0;
            for (int j = 0; j < D; j++) {
                tmp[i] += km->K[i + j*D] * psi_pred[j];
            }
        }
        memcpy(psi_pred, tmp, D * sizeof(float complex));
    }
    
    /* Error in lifted space */
    float complex* psi_err = (float complex*)alloca((size_t)D * sizeof(float complex));
    for (int k = 0; k < D; k++) {
        psi_err[k] = mpc->target_psi[k] - psi_pred[k];
    }
    
    /* Project velocity error back to state space */
    float complex vx_err = 0, vy_err = 0, vz_err = 0;
    for (int k = 0; k < D; k++) {
        vx_err += km->P[3 + k*S] * psi_err[k];
        vy_err += km->P[4 + k*S] * psi_err[k];
        vz_err += km->P[5 + k*S] * psi_err[k];
    }
    
    /* Proportional control: balance state tracking vs control effort */
    float gain = mpc->state_cost / (mpc->state_cost + mpc->control_cost + 1e-10f);
    lg_vec3_t control = {
        gain * crealf(vx_err),
        gain * crealf(vy_err),
        gain * crealf(vz_err)
    };
    
    return control;
}

/*============================================================================
 * Fractional Koopman (Memory-Aware Lifting)
 * 
 * Combines fractional calculus with Koopman: the operator now acts on
 * history-dependent observables. Useful for atmospheric drag with memory
 * or viscoelastic tidal dissipation.
 *===========================================================================*/

typedef struct {
    lg_koopman_t base;
    lg_fractional_t* frac;  /* Fractional history for each observable */
} lg_fractional_koopman_t;

/* Fit fractional Koopman (observables include history).
 * Augments standard observables with power-law weighted history features,
 * capturing subdiffusive memory effects in the spectral domain. */
static inline void lg_fractional_koopman_fit(lg_fractional_koopman_t* fkm,
                                              const lg_koop_data_t* data,
                                              float hurst) {
    fkm->base.hurst = hurst;
    
    int D_orig = fkm->base.n_observables;
    int D_frac = 2 * D_orig;
    int N = data->n_snapshots;
    int S = fkm->base.n_states;
    
    float complex* Psi_X = (float complex*)calloc((size_t)D_frac * (size_t)N, sizeof(float complex));
    float complex* Psi_Y = (float complex*)calloc((size_t)D_frac * (size_t)N, sizeof(float complex));
    
    /* Standard observables (first half) */
    for (int n = 0; n < N; n++) {
        lg_vec3_t pos = {data->X[0 + n*S], data->X[1 + n*S], data->X[2 + n*S]};
        lg_vec3_t vel = {data->X[3 + n*S], data->X[4 + n*S], data->X[5 + n*S]};
        lg_vec3_t pos_y = {data->Y[0 + n*S], data->Y[1 + n*S], data->Y[2 + n*S]};
        lg_vec3_t vel_y = {data->Y[3 + n*S], data->Y[4 + n*S], data->Y[5 + n*S]};
        
        for (int k = 0; k < D_orig; k++) {
            float omega = fkm->base.frequencies[k % D_orig];
            Psi_X[k + n*D_frac] = lg_koop_observable_dct(&pos, &vel, omega, k);
            Psi_Y[k + n*D_frac] = lg_koop_observable_dct(&pos_y, &vel_y, omega, k);
        }
    }
    
    /* Fractional observables (second half): power-law weighted history */
    float alpha = 2.0f - 2.0f * hurst;
    if (alpha <= 0.0f) alpha = 0.1f;
    if (alpha > 1.0f) alpha = 1.0f;
    
    for (int n = 0; n < N; n++) {
        for (int k = 0; k < D_orig; k++) {
            float complex frac_sum_x = 0.0f;
            float complex frac_sum_y = 0.0f;
            float w_sum = 0.0f;
            for (int m = 0; m <= n; m++) {
                float w = powf((float)(n - m + 1), alpha - 1.0f);
                frac_sum_x += w * Psi_X[k + m*D_frac];
                frac_sum_y += w * Psi_Y[k + m*D_frac];
                w_sum += w;
            }
            Psi_X[D_orig + k + n*D_frac] = (w_sum > 0.0f) ? frac_sum_x / w_sum : frac_sum_x;
            Psi_Y[D_orig + k + n*D_frac] = (w_sum > 0.0f) ? frac_sum_y / w_sum : frac_sum_y;
        }
    }
    
    /* Update dimension and fit augmented EDMD */
    fkm->base.n_observables = D_frac;
    _lg_koopman_fit_core(&fkm->base, data, Psi_X, Psi_Y);
    
    /* Save augmented observables */
    free(fkm->base.Psi);
    fkm->base.Psi = Psi_X;
    free(Psi_Y);
}

/* Predict with memory effects: automatically captures drag hysteresis.
 * Updates both transform position and body velocity in-place. */
static inline void lg_fractional_koopman_predict(lg_fractional_koopman_t* fkm,
                                                lg_body_t* body,
                                                lg_transform_t* transform,
                                                int n_steps) {
    /* Update fractional history with current velocity */
    lg_fractional_record(fkm->frac, body);
    
    int D = fkm->base.n_observables;
    int S = fkm->base.n_states;
    int D_orig = D / 2;
    
    lg_vec3_t pos = transform->position;
    lg_vec3_t vel = body->velocity;
    
    /* Compute memory term from fractional integral */
    lg_vec3_t memory = lg_fractional_integral_direct(fkm->frac);
    
    /* Lift to augmented observable space */
    float complex* psi = (float complex*)alloca((size_t)D * sizeof(float complex));
    for (int k = 0; k < D_orig; k++) {
        psi[k] = lg_koop_observable_dct(&pos, &vel, fkm->base.frequencies[k], k);
        lg_vec3_t pos_mem = lg_vec3_add(pos, memory);
        psi[D_orig + k] = lg_koop_observable_dct(&pos_mem, &vel, fkm->base.frequencies[k], k);
    }
    
    /* Apply K^n_steps */
    float complex* psi_next = (float complex*)alloca((size_t)D * sizeof(float complex));
    for (int step = 0; step < n_steps; step++) {
        for (int i = 0; i < D; i++) {
            psi_next[i] = 0;
            for (int j = 0; j < D; j++) {
                psi_next[i] += fkm->base.K[i + j*D] * psi[j];
            }
        }
        memcpy(psi, psi_next, D * sizeof(float complex));
    }
    
    /* Project back to state space */
    float complex px = 0, py = 0, pz = 0;
    float complex vx = 0, vy = 0, vz = 0;
    for (int k = 0; k < D; k++) {
        px += fkm->base.P[0 + k*S] * psi[k];
        py += fkm->base.P[1 + k*S] * psi[k];
        pz += fkm->base.P[2 + k*S] * psi[k];
        vx += fkm->base.P[3 + k*S] * psi[k];
        vy += fkm->base.P[4 + k*S] * psi[k];
        vz += fkm->base.P[5 + k*S] * psi[k];
    }
    
    transform->position.x = crealf(px);
    transform->position.y = crealf(py);
    transform->position.z = crealf(pz);
    body->velocity.x = crealf(vx);
    body->velocity.y = crealf(vy);
    body->velocity.z = crealf(vz);
}

/*============================================================================
 * Spectral/CHEAP Integration
 * 
 * Since CHEAP already computes DCT-II transforms, we use the same basis
 * for Koopman observables—zero overhead if already using CHEAP for
 * fractional convolution.
 *===========================================================================*/

/* Initialize Koopman with CHEAP-compatible DCT basis */
static inline void lg_koopman_init_cheap_compatible(lg_koopman_t* km, int n, float hurst) {
    km->type = LG_KOOP_DCT;
    km->n_observables = n;
    km->n_states = 6;
    km->hurst = hurst;
    
    km->frequencies = (float*)malloc(n * sizeof(float));
    km->scales = (float*)malloc(n * sizeof(float));
    
    /* Use same frequency grid as CHEAP's fBm eigenvalues */
    for (int k = 0; k < n; k++) {
        km->frequencies[k] = LG_PI * (k + 0.5f) / n;  /* DCT-II frequencies */
        km->scales[k] = powf(k + 1.0f, -hurst - 0.5f); /* Power law scaling */
    }
}

/* Fast batch prediction using CHEAP's DCT machinery.
 * Processes each state in the batch independently via Koopman operator.
 * O(N * D^2) total — a true O(N log N) DCT-accelerated version would
 * batch the observable evaluation itself. */
static inline void lg_koopman_predict_batch_cheap(const lg_koopman_t* km,
                                                const lg_vec3_batch_t* pos_batch,
                                                const lg_vec3_batch_t* vel_batch,
                                                lg_vec3_batch_t* pos_out,
                                                int n) {
    for (int i = 0; i < n; i++) {
        lg_vec3_t pos = {pos_batch->x[i], pos_batch->y[i], pos_batch->z[i]};
        lg_vec3_t vel = {vel_batch->x[i], vel_batch->y[i], vel_batch->z[i]};
        lg_vec3_t pos_pred, vel_pred;
        
        lg_koopman_predict(km, &pos, &vel, &pos_pred, &vel_pred, 1);
        
        pos_out->x[i] = pos_pred.x;
        pos_out->y[i] = pos_pred.y;
        pos_out->z[i] = pos_pred.z;
    }
}

#ifdef __cplusplus
}
#endif

#endif /* LAGRANGE_KOOPMAN_H */

