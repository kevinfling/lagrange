/*
 * Lagrange Physics Library - Orbit Regularization
 *
 * Coordinate transforms that remove the 1/r² singularity in Newtonian gravity,
 * enabling stable numerical integration through close encounters.
 *
 *   lg_ks_*     — Kustaanheimo-Stiefel (KS) 3D → 4D spinor regularization
 *   lg_lc_*     — Levi-Civita 2D planar regularization
 *   lg_bf_*     — Burdet-Ferrándiz (radial-distance element) regularization
 *
 * All routines use double precision.  Header-only, no external dependencies.
 *
 * References:
 *   Stiefel & Scheifele, "Linear and Regular Celestial Mechanics", Springer 1971.
 *   Burdet, Celest. Mech. 1(1):222-234, 1969.
 *   Ferrándiz, Celest. Mech. 41:343-357, 1988.
 */

#ifndef LAGRANGE_REGULARIZATION_H
#define LAGRANGE_REGULARIZATION_H

#include <math.h>
#include <stdbool.h>
#include <string.h>

#ifdef __cplusplus
extern "C" {
#endif

/*============================================================================
 * Kustaanheimo-Stiefel (KS) Regularization  — 3D version
 *
 * Transforms Cartesian (r, v) ∈ ℝ³×ℝ³ to spinor (u, u') ∈ ℝ⁴×ℝ⁴ with
 * fictitious time s such that  dt/ds = r  (Sundman substitution).
 *
 * KS matrix L(u):
 *   u = [u0, u1, u2, u3]
 *   r = L(u)^T u  →  r[k] = sum_j L_jk u_j
 *
 * Equations of motion in s:
 *   u'' + (E/2) u = (1/2) L^T f_phys   (f_phys = perturbing force / r)
 *   t' = r = ||u||²
 *   E = -mu/2a  (two-body energy; updated each step for perturbed problems)
 *===========================================================================*/

typedef struct {
    double u[4];   /* KS spinor position */
    double up[4];  /* KS spinor velocity  (d/ds) */
    double t;      /* Physical time */
    double s;      /* Fictitious time */
    double E;      /* Two-body orbital energy (negative for bound orbit) */
    double mu;     /* Gravitational parameter */
} lg_ks_state_t;

/* KS L matrix: maps spinor u → physical r.
 * r = L(u)^T * u  where L is 4×4 with the structure below. */
static inline void _lg_ks_L(const double u[4], double L[4][4]) {
    L[0][0] =  u[0]; L[0][1] = -u[1]; L[0][2] = -u[2]; L[0][3] =  u[3];
    L[1][0] =  u[1]; L[1][1] =  u[0]; L[1][2] = -u[3]; L[1][3] = -u[2];
    L[2][0] =  u[2]; L[2][1] =  u[3]; L[2][2] =  u[0]; L[2][3] =  u[1];
    L[3][0] =  u[3]; L[3][1] = -u[2]; L[3][2] =  u[1]; L[3][3] = -u[0];
}

/* Convert physical (r, v) → KS state.
 * A canonical spinor choice: use the convention that u3 = 0 when possible. */
static inline lg_ks_state_t lg_ks_regularize(
    const double r[3], const double v[3],
    double mu, double t0
) {
    lg_ks_state_t st;
    st.mu = mu;
    st.t  = t0;
    st.s  = 0.0;

    double rn = sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);

    /* Choose spinor u such that L(u)^T u = r, ||u||² = r */
    /* Standard initialization: */
    if (r[0] >= 0.0) {
        st.u[0] = sqrt(0.5*(rn + r[0]));
        double denom = (st.u[0] < 1e-14) ? 1.0 : 2.0*st.u[0];
        st.u[1] =  r[1] / denom;
        st.u[2] =  r[2] / denom;
        st.u[3] =  0.0;
    } else {
        st.u[1] = sqrt(0.5*(rn - r[0]));
        double denom = (st.u[1] < 1e-14) ? 1.0 : 2.0*st.u[1];
        st.u[0] =  r[1] / denom;
        st.u[2] =  0.0;
        st.u[3] =  r[2] / denom;
    }

    /* KS velocity: u' = (1/2r) L(u) v_ext, where v_ext = [vx, vy, vz, 0] */
    double L[4][4];
    _lg_ks_L(st.u, L);
    for (int k = 0; k < 4; k++) {
        st.up[k] = 0.5 / rn * (L[k][0]*v[0] + L[k][1]*v[1] + L[k][2]*v[2]);
    }

    /* Energy */
    double v2 = v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
    st.E = 0.5*v2 - mu/rn;

    return st;
}

/* Convert KS state → physical (r, v). */
static inline void lg_ks_deregularize(
    const lg_ks_state_t* st,
    double r[3], double v[3]
) {
    double L[4][4];
    _lg_ks_L(st->u, L);

    /* r = L^T u  (first 3 components; 4th is always 0 by bilinear identity) */
    for (int k = 0; k < 3; k++) {
        r[k] = 0.0;
        for (int j = 0; j < 4; j++)
            r[k] += L[j][k] * st->u[j];
    }

    double rn = sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);
    if (rn < 1e-30) rn = 1e-30;

    /* v = (2/r) L^T u'  (first 3 components) */
    for (int k = 0; k < 3; k++) {
        v[k] = 0.0;
        for (int j = 0; j < 4; j++)
            v[k] += L[j][k] * st->up[j];
        v[k] *= 2.0 / rn;
    }
}

/* KS perturbed equations of motion RHS.
 * Computes  u'' = -(E/2) u + (1/2) L^T f_phys
 * where f_phys = additional perturbing acceleration (not including central body mu/r²).
 * Also integrates t' = ||u||² and E' = 2 <u', L^T f_phys>.
 *
 * State vector layout (for integrator): [u0,u1,u2,u3, up0,up1,up2,up3, t, E] — dim=10.
 * The caller packs/unpacks via lg_ks_pack / lg_ks_unpack. */
static inline void lg_ks_eom(
    double s, const double* y, double* dyds, int dim, void* user_data
) {
    (void)dim;
    /* Unpack: y[0..3]=u, y[4..7]=u', y[8]=t, y[9]=E */
    const double* u  = y;
    const double* up = y + 4;
    double E  = y[9];
    double mu = *(double*)user_data;

    /* r = ||u||^2 */
    double r = u[0]*u[0] + u[1]*u[1] + u[2]*u[2] + u[3]*u[3];

    /* u'' = -(E/2)*u   (unperturbed; add perturbation via f_phys externally) */
    for (int k = 0; k < 4; k++)
        dyds[k + 4] = -0.5 * E * u[k];

    /* u' */
    for (int k = 0; k < 4; k++)
        dyds[k] = up[k];

    /* t' = r */
    dyds[8] = r;

    /* E' = 0 for pure Kepler; for perturbed, the caller can augment */
    dyds[9] = 0.0;

    (void)s; (void)mu;
}

/* Pack KS state into a flat double array for the ODE integrator. */
static inline void lg_ks_pack(const lg_ks_state_t* st, double y[10]) {
    for (int i = 0; i < 4; i++) y[i]   = st->u[i];
    for (int i = 0; i < 4; i++) y[4+i] = st->up[i];
    y[8] = st->t;
    y[9] = st->E;
}

/* Unpack flat array back into KS state. */
static inline void lg_ks_unpack(lg_ks_state_t* st, const double y[10]) {
    for (int i = 0; i < 4; i++) st->u[i]  = y[i];
    for (int i = 0; i < 4; i++) st->up[i] = y[4+i];
    st->t = y[8];
    st->E = y[9];
}

/* Propagate KS orbit from current state by physical time delta_t.
 * Uses a fixed-step symplectic leapfrog in fictitious time s.
 * n_steps: number of fictitious-time steps (more → more accurate).
 * For unperturbed Kepler, any n_steps works; for perturbations use ≥100.
 *
 * Returns the updated state. */
static inline lg_ks_state_t lg_ks_propagate(
    lg_ks_state_t st,
    double delta_t,
    int n_steps
) {
    double mu = st.mu;

    /* Estimate total fictitious time ds from s→s+ds via t' = r = ||u||²  */
    /* Use dt/ds = r ≈ semi-latus rectum / 2  (rough; iterate for accuracy) */
    double r_est  = st.u[0]*st.u[0] + st.u[1]*st.u[1]
                  + st.u[2]*st.u[2] + st.u[3]*st.u[3];
    double ds_est = (r_est > 1e-20) ? delta_t / r_est : 1e-4;

    double ds = ds_est / n_steps;

    /* Leapfrog (Störmer-Verlet) in s:
     *   u_{n+1/2} = u_n + (ds/2) * u'_n
     *   u'_{n+1}  = u'_n - (ds) * (E/2) * u_{n+1/2}
     *   u_{n+1}   = u_{n+1/2} + (ds/2) * u'_{n+1}
     *   t_{n+1}   = t_n + ds * ||u_{n+1/2}||²
     * Stop when accumulated t reaches delta_t. */

    double t_accumulated = 0.0;
    double target_t = delta_t;

    for (int n = 0; n < n_steps && t_accumulated < target_t; n++) {
        double r_half;
        double u_half[4];

        /* Adaptive: shrink last step to hit target_t exactly */
        double r_now = st.u[0]*st.u[0] + st.u[1]*st.u[1]
                     + st.u[2]*st.u[2] + st.u[3]*st.u[3];
        double dt_remaining = target_t - t_accumulated;
        if (ds * r_now > dt_remaining) {
            ds = dt_remaining / r_now;
        }

        /* Half drift */
        for (int k = 0; k < 4; k++)
            u_half[k] = st.u[k] + 0.5*ds * st.up[k];

        r_half = u_half[0]*u_half[0] + u_half[1]*u_half[1]
               + u_half[2]*u_half[2] + u_half[3]*u_half[3];

        /* Full kick: u'' = (E/2)*u, so up += ds*(E/2)*u_half */
        for (int k = 0; k < 4; k++)
            st.up[k] += ds * 0.5 * st.E * u_half[k];

        /* Half drift */
        for (int k = 0; k < 4; k++)
            st.u[k] = u_half[k] + 0.5*ds * st.up[k];

        /* Update fictitious and physical time */
        st.s += ds;
        t_accumulated += ds * r_half;
    }

    st.t += t_accumulated;

    /* Re-derive E from current state to keep it consistent */
    double r[3], v[3];
    lg_ks_deregularize(&st, r, v);
    double rn = sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);
    double v2 = v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
    if (rn > 1e-30) st.E = 0.5*v2 - mu/rn;

    return st;
}

/*============================================================================
 * Levi-Civita 2D Planar Regularization
 *
 * For the planar two-body problem in (x, y):
 *   r = u² (complex squaring: r + 0i = (u1 + i*u2)²)
 *   Maps to harmonic oscillator with fictitious time τ: dt/dτ = r
 *===========================================================================*/

typedef struct {
    double u[2];   /* LC spinor position */
    double up[2];  /* LC spinor velocity (d/dτ) */
    double tau;    /* Fictitious time */
    double t;      /* Physical time */
    double E;      /* Orbital energy */
    double mu;
} lg_lc_state_t;

/* Convert physical planar (x, y, vx, vy) → LC state. */
static inline lg_lc_state_t lg_lc_regularize(
    double x, double y, double vx, double vy,
    double mu, double t0
) {
    lg_lc_state_t st;
    st.mu  = mu;
    st.t   = t0;
    st.tau = 0.0;

    double r = sqrt(x*x + y*y);
    if (r < 1e-30) r = 1e-30;

    /* u such that u1² - u2² = x, 2*u1*u2 = y */
    if (x >= 0.0) {
        st.u[0] = sqrt(0.5*(r + x));
        st.u[1] = (st.u[0] < 1e-14) ? sqrt(r) : 0.5*y / st.u[0];
    } else {
        st.u[1] = sqrt(0.5*(r - x));
        st.u[0] = (st.u[1] < 1e-14) ? sqrt(r) : 0.5*y / st.u[1];
    }

    /* LC matrix A = [[u1, -u2], [u2, u1]], velocity u' = (1/2r) A^T v */
    st.up[0] = 0.5/r * ( st.u[0]*vx + st.u[1]*vy);
    st.up[1] = 0.5/r * (-st.u[1]*vx + st.u[0]*vy);

    double v2 = vx*vx + vy*vy;
    st.E = 0.5*v2 - mu/r;
    return st;
}

/* Convert LC state → physical (x, y, vx, vy). */
static inline void lg_lc_deregularize(
    const lg_lc_state_t* st,
    double* x, double* y, double* vx, double* vy
) {
    *x = st->u[0]*st->u[0] - st->u[1]*st->u[1];
    *y = 2.0 * st->u[0] * st->u[1];
    double r = st->u[0]*st->u[0] + st->u[1]*st->u[1];
    if (r < 1e-30) r = 1e-30;
    *vx = 2.0/r * ( st->u[0]*st->up[0] - st->u[1]*st->up[1]);
    *vy = 2.0/r * ( st->u[1]*st->up[0] + st->u[0]*st->up[1]);
}

/* Leapfrog propagation in LC fictitious time. */
static inline lg_lc_state_t lg_lc_propagate(
    lg_lc_state_t st,
    double delta_t,
    int n_steps
) {
    double r_est = st.u[0]*st.u[0] + st.u[1]*st.u[1];
    double ds    = (r_est > 1e-20) ? (delta_t / r_est) / n_steps : 1e-4;

    double t_acc = 0.0;
    double target_t = delta_t;

    for (int n = 0; n < n_steps && t_acc < target_t; n++) {
        double r_now = st.u[0]*st.u[0] + st.u[1]*st.u[1];
        if (ds * r_now > target_t - t_acc)
            ds = (target_t - t_acc) / r_now;

        double u_half[2];
        u_half[0] = st.u[0] + 0.5*ds * st.up[0];
        u_half[1] = st.u[1] + 0.5*ds * st.up[1];

        double r_h = u_half[0]*u_half[0] + u_half[1]*u_half[1];

        /* u'' = (E/2)*u for LC as well */
        st.up[0] += ds * 0.5 * st.E * u_half[0];
        st.up[1] += ds * 0.5 * st.E * u_half[1];

        st.u[0] = u_half[0] + 0.5*ds * st.up[0];
        st.u[1] = u_half[1] + 0.5*ds * st.up[1];

        st.tau += ds;
        t_acc  += ds * r_h;
    }

    st.t += t_acc;

    /* Update energy */
    double x, y, vx, vy;
    lg_lc_deregularize(&st, &x, &y, &vx, &vy);
    double r = sqrt(x*x + y*y);
    if (r > 1e-30) st.E = 0.5*(vx*vx+vy*vy) - st.mu/r;

    return st;
}

/*============================================================================
 * Burdet-Ferrándiz Regularization
 *
 * Uses radial distance r as an independent element with the time element
 * as a state variable.  State: (r_vec, v_vec, t, r, r') ∈ ℝ^9.
 * EOM expressed in fictitious time s via  dt/ds = r.
 *
 * Equations (unperturbed):
 *   r⃗'' + (1/r³)(mu - r*r'²  + r*h_dot) * r⃗ = 0    ... approximation
 * We use the simpler radial element form:
 *   r'' + ω² r = mu/r²  + ... (angular momentum terms)
 *   where ω² = h²/r⁴  (h = angular momentum magnitude)
 *===========================================================================*/

typedef struct {
    double q[3];   /* position vector */
    double qp[3];  /* velocity vector (d/ds) */
    double s;      /* fictitious time */
    double t;      /* physical time */
    double h_sq;   /* specific angular momentum squared (conserved) */
    double E;      /* energy */
    double mu;
} lg_bf_state_t;

/* Initialize BF state from Cartesian (r, v). */
static inline lg_bf_state_t lg_bf_regularize(
    const double r[3], const double v[3],
    double mu, double t0
) {
    lg_bf_state_t st;
    st.mu = mu;
    st.t  = t0;
    st.s  = 0.0;

    double rn = sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);

    for (int k = 0; k < 3; k++) {
        st.q[k]  = r[k];
        /* BF: q' = r*v (physical velocity scaled by r) */
        st.qp[k] = rn * v[k];
    }

    /* Angular momentum h = r × v */
    double hx = r[1]*v[2] - r[2]*v[1];
    double hy = r[2]*v[0] - r[0]*v[2];
    double hz = r[0]*v[1] - r[1]*v[0];
    st.h_sq = hx*hx + hy*hy + hz*hz;

    double v2 = v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
    st.E = 0.5*v2 - mu/rn;
    return st;
}

/* Convert BF state → physical (r, v). */
static inline void lg_bf_deregularize(
    const lg_bf_state_t* st,
    double r[3], double v[3]
) {
    double rn = sqrt(st->q[0]*st->q[0] + st->q[1]*st->q[1] + st->q[2]*st->q[2]);
    if (rn < 1e-30) rn = 1e-30;
    for (int k = 0; k < 3; k++) {
        r[k] = st->q[k];
        v[k] = st->qp[k] / rn;
    }
}

/* One leapfrog step in BF fictitious time (simplified, unperturbed). */
static inline void _lg_bf_step(lg_bf_state_t* st, double ds) {
    /* Compute current r */
    double rn = sqrt(st->q[0]*st->q[0] + st->q[1]*st->q[1] + st->q[2]*st->q[2]);

    /* BF equations of motion (unperturbed):
     *   q'' = (2E + mu/r) * q   ( = -(omega²)*q with omega² = -2E - mu/r )
     * This is only valid for bound orbits (E < 0). */
    double omega2 = -(2.0*st->E + st->mu/rn);

    /* Half-drift */
    double q_half[3];
    for (int k = 0; k < 3; k++)
        q_half[k] = st->q[k] + 0.5*ds * st->qp[k];

    double rn_half = sqrt(q_half[0]*q_half[0] + q_half[1]*q_half[1] + q_half[2]*q_half[2]);
    double omega2_half = -(2.0*st->E + st->mu/rn_half);
    if (rn_half < 1e-30) omega2_half = omega2;

    /* Full kick */
    for (int k = 0; k < 3; k++)
        st->qp[k] += ds * omega2_half * q_half[k];

    /* Half-drift */
    for (int k = 0; k < 3; k++)
        st->q[k] = q_half[k] + 0.5*ds * st->qp[k];

    /* Physical time increment: dt = r * ds */
    st->t += rn_half * ds;
    st->s += ds;
}

/* Propagate BF orbit by physical time delta_t. */
static inline lg_bf_state_t lg_bf_propagate(
    lg_bf_state_t st,
    double delta_t,
    int n_steps
) {
    double rn_est = sqrt(st.q[0]*st.q[0] + st.q[1]*st.q[1] + st.q[2]*st.q[2]);
    double ds = (rn_est > 1e-20) ? (delta_t / rn_est) / n_steps : 1e-4;

    double t_target = st.t + delta_t;

    for (int n = 0; n < n_steps; n++) {
        double rn = sqrt(st.q[0]*st.q[0] + st.q[1]*st.q[1] + st.q[2]*st.q[2]);
        double dt_left = t_target - st.t;
        if (dt_left <= 0.0) break;
        if (ds * rn > dt_left) ds = dt_left / rn;
        _lg_bf_step(&st, ds);
    }

    /* Update energy from current state */
    double r[3], v[3];
    lg_bf_deregularize(&st, r, v);
    double rn = sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);
    double v2 = v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
    if (rn > 1e-30) st.E = 0.5*v2 - st.mu/rn;

    return st;
}

/*============================================================================
 * Near-singularity detection
 *===========================================================================*/

/* Returns true if Cartesian orbit will approach within r_threshold on the
 * next Kepler half-period (crude periapsis estimate from energy and h). */
static inline bool lg_orbit_near_singular(
    const double r[3], const double v[3],
    double mu, double r_threshold
) {
    double rn = sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);
    double v2 = v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
    double E  = 0.5*v2 - mu/rn;
    if (E >= 0.0) return false; /* Hyperbolic, no bound periapsis */
    double a = -0.5*mu/E;
    /* Angular momentum */
    double hx = r[1]*v[2] - r[2]*v[1];
    double hy = r[2]*v[0] - r[0]*v[2];
    double hz = r[0]*v[1] - r[1]*v[0];
    double h2 = hx*hx + hy*hy + hz*hz;
    double p  = h2 / mu;            /* semi-latus rectum */
    double e  = sqrt(1.0 - p/a);    /* eccentricity */
    if (e > 1.0) return false;
    double r_peri = a * (1.0 - e);  /* periapsis */
    return r_peri < r_threshold;
}

#ifdef __cplusplus
}
#endif

#endif /* LAGRANGE_REGULARIZATION_H */
