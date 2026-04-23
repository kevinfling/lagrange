/*
 * Lagrange Physics Library - Automatic Differentiation & Taylor Propagation
 *
 * Provides:
 *   lg_dual_t          — dual number for forward-mode AD (1st derivative)
 *   lg_dual2_t         — hyper-dual number (2nd derivative)
 *   lg_dual3_t         — 3D dual-vector for state-space Jacobians
 *   lg_stm_dual_*      — State Transition Matrix (STM) via dual lifting
 *   lg_taylor2_cr3bp_step — 2nd-order Taylor propagator for CR3BP
 *   lg_taylor3_cr3bp_step — 3rd-order Taylor propagator for CR3BP
 *
 * All routines are header-only, static inline, C11, no external dependencies.
 * Double precision throughout.
 *
 * References:
 *   Fike & Alonso, "The Development of Hyper-Dual Numbers for Exact
 *     Second-Derivative Calculations", AIAA 2011.
 *   Benettin & Giorgilli, "On the Hamiltonian interpolation...", J. Stat. Phys. 1994.
 */

#ifndef LAGRANGE_AUTODIFF_H
#define LAGRANGE_AUTODIFF_H

#include <math.h>
#include <string.h>
#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif

/*============================================================================
 * Dual Number — Forward-Mode AD, 1st Derivative
 *
 * d = (real, eps) where eps² = 0.
 * f(x + ε) = f(x) + f'(x)ε  →  dual part gives exact 1st derivative.
 *===========================================================================*/

typedef struct {
    double r; /* real part */
    double e; /* dual (epsilon) part */
} lg_dual_t;

static inline lg_dual_t lg_dual(double r, double e) {
    lg_dual_t d = {r, e};
    return d;
}

/* Seed a variable: x with derivative 1 */
static inline lg_dual_t lg_dual_var(double x)      { lg_dual_t d = {x, 1.0}; return d; }
/* Seed a constant: c with derivative 0 */
static inline lg_dual_t lg_dual_const(double c)    { lg_dual_t d = {c, 0.0}; return d; }

static inline lg_dual_t lg_dual_add(lg_dual_t a, lg_dual_t b) {
    lg_dual_t d = {a.r + b.r, a.e + b.e}; return d;
}
static inline lg_dual_t lg_dual_sub(lg_dual_t a, lg_dual_t b) {
    lg_dual_t d = {a.r - b.r, a.e - b.e}; return d;
}
static inline lg_dual_t lg_dual_mul(lg_dual_t a, lg_dual_t b) {
    lg_dual_t d = {a.r*b.r, a.r*b.e + a.e*b.r}; return d;
}
static inline lg_dual_t lg_dual_div(lg_dual_t a, lg_dual_t b) {
    lg_dual_t d = {a.r/b.r, (a.e*b.r - a.r*b.e)/(b.r*b.r)}; return d;
}
static inline lg_dual_t lg_dual_neg(lg_dual_t a) {
    lg_dual_t d = {-a.r, -a.e}; return d;
}
static inline lg_dual_t lg_dual_scale(lg_dual_t a, double s) {
    lg_dual_t d = {a.r*s, a.e*s}; return d;
}
static inline lg_dual_t lg_dual_sqrt(lg_dual_t a) {
    double sr = sqrt(a.r);
    lg_dual_t d = {sr, a.e / (2.0*sr)}; return d;
}
static inline lg_dual_t lg_dual_sq(lg_dual_t a) {
    lg_dual_t d = {a.r*a.r, 2.0*a.r*a.e}; return d;
}
static inline lg_dual_t lg_dual_sin(lg_dual_t a) {
    lg_dual_t d = {sin(a.r), a.e * cos(a.r)}; return d;
}
static inline lg_dual_t lg_dual_cos(lg_dual_t a) {
    lg_dual_t d = {cos(a.r), -a.e * sin(a.r)}; return d;
}
static inline lg_dual_t lg_dual_exp(lg_dual_t a) {
    double ea = exp(a.r);
    lg_dual_t d = {ea, a.e * ea}; return d;
}
static inline lg_dual_t lg_dual_log(lg_dual_t a) {
    lg_dual_t d = {log(a.r), a.e / a.r}; return d;
}
static inline lg_dual_t lg_dual_pow(lg_dual_t a, double p) {
    double rp = pow(a.r, p);
    lg_dual_t d = {rp, a.e * p * pow(a.r, p-1.0)}; return d;
}
static inline lg_dual_t lg_dual_abs(lg_dual_t a) {
    lg_dual_t d = {fabs(a.r), (a.r >= 0.0) ? a.e : -a.e}; return d;
}
static inline lg_dual_t lg_dual_inv(lg_dual_t a) {
    lg_dual_t d = {1.0/a.r, -a.e/(a.r*a.r)}; return d;
}

/*============================================================================
 * Dual 3-Vector — for Jacobian of 3D functions
 *===========================================================================*/

typedef struct {
    lg_dual_t x, y, z;
} lg_dual3_t;

static inline lg_dual3_t lg_dual3(double x, double y, double z,
                                   double dx, double dy, double dz) {
    lg_dual3_t v;
    v.x = lg_dual(x, dx);
    v.y = lg_dual(y, dy);
    v.z = lg_dual(z, dz);
    return v;
}

static inline lg_dual_t lg_dual3_dot(lg_dual3_t a, lg_dual3_t b) {
    return lg_dual_add(lg_dual_add(lg_dual_mul(a.x,b.x), lg_dual_mul(a.y,b.y)), lg_dual_mul(a.z,b.z));
}

static inline lg_dual_t lg_dual3_norm(lg_dual3_t v) {
    return lg_dual_sqrt(lg_dual3_dot(v, v));
}

static inline lg_dual3_t lg_dual3_add(lg_dual3_t a, lg_dual3_t b) {
    lg_dual3_t r;
    r.x = lg_dual_add(a.x, b.x);
    r.y = lg_dual_add(a.y, b.y);
    r.z = lg_dual_add(a.z, b.z);
    return r;
}

static inline lg_dual3_t lg_dual3_sub(lg_dual3_t a, lg_dual3_t b) {
    lg_dual3_t r;
    r.x = lg_dual_sub(a.x, b.x);
    r.y = lg_dual_sub(a.y, b.y);
    r.z = lg_dual_sub(a.z, b.z);
    return r;
}

static inline lg_dual3_t lg_dual3_scale(lg_dual3_t v, lg_dual_t s) {
    lg_dual3_t r;
    r.x = lg_dual_mul(v.x, s);
    r.y = lg_dual_mul(v.y, s);
    r.z = lg_dual_mul(v.z, s);
    return r;
}

/*============================================================================
 * STM (State Transition Matrix) via Dual-Number Lifting
 *
 * For an ODE  dx/dt = f(x),  the STM Phi satisfies  dPhi/dt = Df(x) * Phi.
 * We compute Df(x) * e_k for each unit seed e_k by evaluating f with
 * x perturbed by a dual number in the k-th direction.
 *
 * State dimension n <= LG_STM_MAX_DIM.
 * Phi is stored row-major: Phi[i*n + j] = dxi/dx0_j.
 *===========================================================================*/

#define LG_STM_MAX_DIM 6

/* User-supplied function: dual-lifted version of the RHS.
 * Takes dual state (dim) and writes dual derivative (dim).
 * Called once per column of Phi (n_seeds times per step). */
typedef void (*lg_dual_rhs_t)(double t,
                               const lg_dual_t* x_dual,
                               lg_dual_t* dxdt_dual,
                               int dim,
                               void* user_data);

/* Advance state and STM by dt using RK4 with dual-number Jacobian columns.
 * x[dim]: current state (real)
 * Phi[dim*dim]: current STM (identity at t0)
 * After call, x and Phi are advanced by dt. */
static inline void lg_stm_dual_step_rk4(
    lg_dual_rhs_t f,
    double t, double* x, double* Phi,
    int dim, double dt,
    void* user_data
) {
    /* Advance state with real-valued RK4 first */
    double k1[LG_STM_MAX_DIM], k2[LG_STM_MAX_DIM];
    double k3[LG_STM_MAX_DIM], k4[LG_STM_MAX_DIM];
    double xtmp[LG_STM_MAX_DIM];
    lg_dual_t xd[LG_STM_MAX_DIM], dxdt[LG_STM_MAX_DIM];

    /* k1 */
    for (int i = 0; i < dim; i++) xd[i] = lg_dual_const(x[i]);
    f(t, xd, dxdt, dim, user_data);
    for (int i = 0; i < dim; i++) k1[i] = dxdt[i].r;

    /* k2 */
    for (int i = 0; i < dim; i++) {
        xtmp[i] = x[i] + 0.5*dt*k1[i];
        xd[i]   = lg_dual_const(xtmp[i]);
    }
    f(t + 0.5*dt, xd, dxdt, dim, user_data);
    for (int i = 0; i < dim; i++) k2[i] = dxdt[i].r;

    /* k3 */
    for (int i = 0; i < dim; i++) {
        xtmp[i] = x[i] + 0.5*dt*k2[i];
        xd[i]   = lg_dual_const(xtmp[i]);
    }
    f(t + 0.5*dt, xd, dxdt, dim, user_data);
    for (int i = 0; i < dim; i++) k3[i] = dxdt[i].r;

    /* k4 */
    for (int i = 0; i < dim; i++) {
        xtmp[i] = x[i] + dt*k3[i];
        xd[i]   = lg_dual_const(xtmp[i]);
    }
    f(t + dt, xd, dxdt, dim, user_data);
    for (int i = 0; i < dim; i++) k4[i] = dxdt[i].r;

    /* Update state */
    double x_new[LG_STM_MAX_DIM];
    for (int i = 0; i < dim; i++)
        x_new[i] = x[i] + dt/6.0 * (k1[i] + 2.0*k2[i] + 2.0*k3[i] + k4[i]);

    /* Advance each column of Phi by computing Df(x)*Phi_col via dual numbers.
     * For column j: seed x with dual part = Phi_col_j, propagate one RK4 step. */
    double Phi_new[LG_STM_MAX_DIM * LG_STM_MAX_DIM];

    for (int col = 0; col < dim; col++) {
        /* RK4 for the variational equation Phi_dot = Df(x) * Phi_col */
        /* We compute Df(x)*v for v = Phi_col using one dual evaluation per stage */

        double phi_col[LG_STM_MAX_DIM];
        for (int i = 0; i < dim; i++) phi_col[i] = Phi[i*dim + col];

        /* Stage evaluations at the 4 state midpoints already computed above */
        double stages[4][LG_STM_MAX_DIM]; /* x_stage */
        double x_stages[4][LG_STM_MAX_DIM];
        double t_stages[4];
        {
            t_stages[0] = t;
            t_stages[1] = t + 0.5*dt;
            t_stages[2] = t + 0.5*dt;
            t_stages[3] = t + dt;
            for (int i = 0; i < dim; i++) {
                x_stages[0][i] = x[i];
                x_stages[1][i] = x[i] + 0.5*dt*k1[i];
                x_stages[2][i] = x[i] + 0.5*dt*k2[i];
                x_stages[3][i] = x[i] + dt*k3[i];
            }
        }

        double dphi[4][LG_STM_MAX_DIM];

        /* Stage 0: seed = phi_col */
        for (int i = 0; i < dim; i++)
            xd[i] = lg_dual(x_stages[0][i], phi_col[i]);
        f(t_stages[0], xd, dxdt, dim, user_data);
        for (int i = 0; i < dim; i++) dphi[0][i] = dxdt[i].e;

        /* Stage 1: seed = phi_col + 0.5*dt*dphi[0] */
        for (int i = 0; i < dim; i++)
            xd[i] = lg_dual(x_stages[1][i], phi_col[i] + 0.5*dt*dphi[0][i]);
        f(t_stages[1], xd, dxdt, dim, user_data);
        for (int i = 0; i < dim; i++) dphi[1][i] = dxdt[i].e;

        /* Stage 2: seed = phi_col + 0.5*dt*dphi[1] */
        for (int i = 0; i < dim; i++)
            xd[i] = lg_dual(x_stages[2][i], phi_col[i] + 0.5*dt*dphi[1][i]);
        f(t_stages[2], xd, dxdt, dim, user_data);
        for (int i = 0; i < dim; i++) dphi[2][i] = dxdt[i].e;

        /* Stage 3: seed = phi_col + dt*dphi[2] */
        for (int i = 0; i < dim; i++)
            xd[i] = lg_dual(x_stages[3][i], phi_col[i] + dt*dphi[2][i]);
        f(t_stages[3], xd, dxdt, dim, user_data);
        for (int i = 0; i < dim; i++) dphi[3][i] = dxdt[i].e;

        for (int i = 0; i < dim; i++) {
            stages[0][i] = dphi[0][i];
            stages[1][i] = dphi[1][i];
            stages[2][i] = dphi[2][i];
            stages[3][i] = dphi[3][i];
            Phi_new[i*dim + col] = phi_col[i]
                + dt/6.0*(stages[0][i] + 2.0*stages[1][i] + 2.0*stages[2][i] + stages[3][i]);
        }
    }

    memcpy(x, x_new, sizeof(double)*dim);
    memcpy(Phi, Phi_new, sizeof(double)*dim*dim);
}

/* Initialize the STM to the identity matrix. */
static inline void lg_stm_init(double* Phi, int dim) {
    for (int i = 0; i < dim; i++)
        for (int j = 0; j < dim; j++)
            Phi[i*dim + j] = (i == j) ? 1.0 : 0.0;
}

/* Finite-difference Jacobian for comparison/testing.
 * J[i*dim + j] = d(f_i)/d(x_j) at x. */
typedef void (*lg_real_rhs_t)(double t, const double* x, double* dxdt,
                              int dim, void* user_data);

static inline void lg_fd_jacobian(
    lg_real_rhs_t f,
    double t, const double* x, double* J,
    int dim, double h_fd, void* user_data
) {
    double f0[LG_STM_MAX_DIM], fp[LG_STM_MAX_DIM], fm[LG_STM_MAX_DIM];
    double xp[LG_STM_MAX_DIM], xm[LG_STM_MAX_DIM];
    f(t, x, f0, dim, user_data);
    (void)f0;
    for (int j = 0; j < dim; j++) {
        memcpy(xp, x, sizeof(double)*dim);
        memcpy(xm, x, sizeof(double)*dim);
        xp[j] += h_fd;
        xm[j] -= h_fd;
        f(t, xp, fp, dim, user_data);
        f(t, xm, fm, dim, user_data);
        for (int i = 0; i < dim; i++)
            J[i*dim + j] = (fp[i] - fm[i]) / (2.0*h_fd);
    }
}

/*============================================================================
 * CR3BP — Circular Restricted Three-Body Problem
 *
 * Non-dimensional units: length = L (Earth-Moon distance),
 *   time = 1/n (n = mean motion), mass = M1+M2.
 * mu = M2/(M1+M2), 1-mu = M1/(M1+M2).
 *
 * State: x = [x, y, z, xdot, ydot, zdot]  (6D)
 * EOM:
 *   x'' - 2*ydot = x - (1-mu)*(x+mu)/r1³ - mu*(x-1+mu)/r2³
 *   y'' + 2*xdot = y - (1-mu)*y/r1³ - mu*y/r2³
 *   z''           = -(1-mu)*z/r1³ - mu*z/r2³
 *   r1 = sqrt((x+mu)²+y²+z²)
 *   r2 = sqrt((x-1+mu)²+y²+z²)
 *===========================================================================*/

/* Plain CR3BP RHS — real-valued (for reference / energy checks). */
static inline void lg_cr3bp_rhs(
    double t, const double* y, double* dydt,
    int dim, void* user_data
) {
    (void)t; (void)dim;
    double mu = *(double*)user_data;
    double x  = y[0], yc = y[1], z = y[2];
    double xd = y[3], yd = y[4], zd = y[5];

    double r1 = sqrt((x+mu)*(x+mu) + yc*yc + z*z);
    double r2 = sqrt((x-1.0+mu)*(x-1.0+mu) + yc*yc + z*z);
    double r1_3 = r1*r1*r1;
    double r2_3 = r2*r2*r2;

    dydt[0] = xd;
    dydt[1] = yd;
    dydt[2] = zd;
    dydt[3] =  2.0*yd + x - (1.0-mu)*(x+mu)/r1_3 - mu*(x-1.0+mu)/r2_3;
    dydt[4] = -2.0*xd + yc - (1.0-mu)*yc/r1_3 - mu*yc/r2_3;
    dydt[5] = -(1.0-mu)*z/r1_3 - mu*z/r2_3;
}

/* Dual-number CR3BP RHS for STM computation. */
static inline void lg_cr3bp_rhs_dual(
    double t,
    const lg_dual_t* y, lg_dual_t* dydt,
    int dim, void* user_data
) {
    (void)t; (void)dim;
    double mu_r = *(double*)user_data;
    lg_dual_t mu   = lg_dual_const(mu_r);
    lg_dual_t one  = lg_dual_const(1.0);
    lg_dual_t two  = lg_dual_const(2.0);

    lg_dual_t x  = y[0], yc = y[1], z = y[2];
    lg_dual_t xd = y[3], yd = y[4], zd = y[5];

    /* r1 = sqrt((x+mu)²+y²+z²) */
    lg_dual_t xpmu = lg_dual_add(x, mu);
    lg_dual_t xm1pmu = lg_dual_sub(x, lg_dual_sub(one, mu));

    lg_dual_t r1_2 = lg_dual_add(lg_dual_add(lg_dual_sq(xpmu), lg_dual_sq(yc)), lg_dual_sq(z));
    lg_dual_t r2_2 = lg_dual_add(lg_dual_add(lg_dual_sq(xm1pmu), lg_dual_sq(yc)), lg_dual_sq(z));

    lg_dual_t r1   = lg_dual_sqrt(r1_2);
    lg_dual_t r2   = lg_dual_sqrt(r2_2);
    lg_dual_t r1_3 = lg_dual_mul(r1_2, r1);
    lg_dual_t r2_3 = lg_dual_mul(r2_2, r2);

    lg_dual_t one_m_mu = lg_dual_sub(one, mu);

    dydt[0] = xd;
    dydt[1] = yd;
    dydt[2] = zd;

    /* xdd = 2*yd + x - (1-mu)*(x+mu)/r1³ - mu*(x-1+mu)/r2³ */
    dydt[3] = lg_dual_sub(
                lg_dual_add(lg_dual_mul(two, yd), x),
                lg_dual_add(
                    lg_dual_div(lg_dual_mul(one_m_mu, xpmu), r1_3),
                    lg_dual_div(lg_dual_mul(mu, xm1pmu), r2_3)
                )
              );

    /* ydd = -2*xd + y - (1-mu)*y/r1³ - mu*y/r2³ */
    dydt[4] = lg_dual_sub(
                lg_dual_add(lg_dual_neg(lg_dual_mul(two, xd)), yc),
                lg_dual_add(
                    lg_dual_div(lg_dual_mul(one_m_mu, yc), r1_3),
                    lg_dual_div(lg_dual_mul(mu, yc), r2_3)
                )
              );

    /* zdd = -(1-mu)*z/r1³ - mu*z/r2³ */
    dydt[5] = lg_dual_neg(
                lg_dual_add(
                    lg_dual_div(lg_dual_mul(one_m_mu, z), r1_3),
                    lg_dual_div(lg_dual_mul(mu, z), r2_3)
                )
              );
}

/* Jacobi constant (energy-like conserved quantity in CR3BP). */
static inline double lg_cr3bp_jacobi(
    double mu,
    const double y[6]
) {
    double x = y[0], yc = y[1], z = y[2];
    double xd = y[3], yd = y[4], zd = y[5];
    double r1 = sqrt((x+mu)*(x+mu) + yc*yc + z*z);
    double r2 = sqrt((x-1.0+mu)*(x-1.0+mu) + yc*yc + z*z);
    double v2 = xd*xd + yd*yd + zd*zd;
    double r2_sq = x*x + yc*yc;
    /* C = 2*Omega(x,y) - v²   where Omega = 0.5*(x²+y²) + (1-mu)/r1 + mu/r2 */
    double Omega = 0.5*r2_sq + (1.0-mu)/r1 + mu/r2;
    return 2.0*Omega - v2;
}

/*============================================================================
 * Taylor Series Propagator — Order 2 (CR3BP)
 *
 * Advances state y by dt using 2nd-order Taylor:
 *   y(t+dt) ≈ y(t) + dt*f(t,y) + 0.5*dt²*Df(t,y)*f(t,y)
 *
 * The 2nd derivative term Df*f is computed via dual-number evaluation:
 *   d/dε [f(x + ε*f)] at ε=0 = Df(x)*f(x)
 *===========================================================================*/

static inline void lg_taylor2_cr3bp_step(
    double* y, double mu, double dt
) {
    double f0[6], f1[6];

    /* Compute f(y) */
    lg_cr3bp_rhs(0.0, y, f0, 6, &mu);

    /* Compute Df(y)*f(y) via dual: evaluate f at y + ε*f0 */
    lg_dual_t yd[6], dydt_d[6];
    for (int i = 0; i < 6; i++)
        yd[i] = lg_dual(y[i], f0[i]);
    lg_cr3bp_rhs_dual(0.0, yd, dydt_d, 6, &mu);
    for (int i = 0; i < 6; i++)
        f1[i] = dydt_d[i].e; /* = Df*f */

    /* y_new = y + dt*f0 + 0.5*dt²*f1 */
    for (int i = 0; i < 6; i++)
        y[i] += dt*f0[i] + 0.5*dt*dt*f1[i];
}

/*============================================================================
 * Taylor Series Propagator — Order 3 (CR3BP)
 *
 * y(t+dt) ≈ y + dt*f + (dt²/2)*f' + (dt³/6)*f''
 * where f' = Df*f  and  f'' = D²f*(f,f) + Df*f'
 *
 * We compute f'' via a second dual evaluation: d/dε[Df(x+εf)*f] at ε=0
 * using hyper-dual arithmetic (or nested dual).
 * We use a simpler finite-difference approximation for the 3rd term to
 * avoid the complexity of hyper-dual: f'' ≈ (f'(x+h*f) - f'(x-h*f))/(2h)
 *===========================================================================*/

static inline void lg_taylor3_cr3bp_step(
    double* y, double mu, double dt
) {
    double f0[6], f1[6], f2[6];

    /* f = f(y) */
    lg_cr3bp_rhs(0.0, y, f0, 6, &mu);

    /* f' = Df*f */
    {
        lg_dual_t yd[6], dydt_d[6];
        for (int i = 0; i < 6; i++) yd[i] = lg_dual(y[i], f0[i]);
        lg_cr3bp_rhs_dual(0.0, yd, dydt_d, 6, &mu);
        for (int i = 0; i < 6; i++) f1[i] = dydt_d[i].e;
    }

    /* f'' ≈ D(Df*f)*f  via central difference in the f direction */
    {
        const double h_fd = 1e-7;
        double yp[6], ym[6];
        double f1p[6], f1m[6];
        for (int i = 0; i < 6; i++) {
            yp[i] = y[i] + h_fd * f0[i];
            ym[i] = y[i] - h_fd * f0[i];
        }
        /* Compute Df(yp)*f(yp) */
        {
            double fp[6];
            lg_cr3bp_rhs(0.0, yp, fp, 6, &mu);
            lg_dual_t yd[6], dydt_d[6];
            for (int i = 0; i < 6; i++) yd[i] = lg_dual(yp[i], fp[i]);
            lg_cr3bp_rhs_dual(0.0, yd, dydt_d, 6, &mu);
            for (int i = 0; i < 6; i++) f1p[i] = dydt_d[i].e;
        }
        {
            double fm[6];
            lg_cr3bp_rhs(0.0, ym, fm, 6, &mu);
            lg_dual_t yd[6], dydt_d[6];
            for (int i = 0; i < 6; i++) yd[i] = lg_dual(ym[i], fm[i]);
            lg_cr3bp_rhs_dual(0.0, yd, dydt_d, 6, &mu);
            for (int i = 0; i < 6; i++) f1m[i] = dydt_d[i].e;
        }
        for (int i = 0; i < 6; i++)
            f2[i] = (f1p[i] - f1m[i]) / (2.0 * h_fd);
    }

    /* y_new = y + dt*f + dt²/2 * f' + dt³/6 * f'' */
    for (int i = 0; i < 6; i++)
        y[i] += dt*f0[i] + (dt*dt/2.0)*f1[i] + (dt*dt*dt/6.0)*f2[i];
}

/*============================================================================
 * CR3BP STM — Variational Equations via Dual Lifting
 *
 * Integrates both the state y(t) and the STM Phi(t; t0) simultaneously.
 * Uses the dual-number STM stepper lg_stm_dual_step_rk4.
 *===========================================================================*/

static inline void lg_cr3bp_stm_step(
    double* y, double* Phi,
    double mu, double dt
) {
    lg_stm_dual_step_rk4(lg_cr3bp_rhs_dual, 0.0, y, Phi, 6, dt, &mu);
}

/* Propagate CR3BP state + STM from t0 to t0+T using n_steps fixed steps. */
static inline void lg_cr3bp_stm_propagate(
    double* y, double* Phi,
    double mu, double T, int n_steps
) {
    lg_stm_init(Phi, 6);
    double dt = T / n_steps;
    double t  = 0.0;
    for (int n = 0; n < n_steps; n++) {
        lg_cr3bp_stm_step(y, Phi, mu, dt);
        t += dt;
    }
    (void)t;
}

#ifdef __cplusplus
}
#endif

#endif /* LAGRANGE_AUTODIFF_H */
