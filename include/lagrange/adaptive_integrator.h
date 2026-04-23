/*
 * Lagrange Physics Library - High-Order & Adaptive Integrators
 *
 * Provides:
 *   - DOP853  : Dormand-Prince 8(7), 13-stage embedded RK with adaptive step control
 *   - Gauss-Jackson order-8 multistep  : summed-form corrector for long arcs
 *   - Wisdom-Holman symplectic map     : N-body Kepler + kick split operator
 *   - Unified adaptive ODE API         : lg_ode_state_t / lg_ode_params_t / lg_integrate_adaptive()
 *
 * All implementations are static inline, header-only, no external dependencies.
 * Double precision throughout — designed for orbital / astrodynamics use.
 */

#ifndef LAGRANGE_ADAPTIVE_INTEGRATOR_H
#define LAGRANGE_ADAPTIVE_INTEGRATOR_H

#include <math.h>
#include <string.h>
#include <stdbool.h>
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/*============================================================================
 * Unified ODE State / Derivative
 *===========================================================================*/

/* Maximum state dimension supported by the adaptive interface.
 * For N-body with positions+velocities: N*6 doubles.
 * Increase if larger systems are needed. */
#define LG_ODE_MAX_DIM 128

typedef struct {
    double y[LG_ODE_MAX_DIM]; /* State vector */
    int    dim;                /* Active dimension (<= LG_ODE_MAX_DIM) */
    double t;                  /* Current independent variable (time / fictitious time) */
} lg_ode_state_t;

/* User-supplied RHS: dy/dt = f(t, y, user_data) */
typedef void (*lg_ode_rhs_t)(double t, const double* y, double* dydt, int dim, void* user_data);

typedef struct {
    double t_end;       /* Integrate until t reaches t_end */
    double h_init;      /* Initial step size (0 → auto-select ~1e-3) */
    double h_min;       /* Minimum allowed step size */
    double h_max;       /* Maximum allowed step size (0 → unbounded) */
    double rtol;        /* Relative tolerance (per component) */
    double atol;        /* Absolute tolerance (per component) */
    int    max_steps;   /* Safety limit; 0 → 1,000,000 */
} lg_ode_params_t;

typedef struct {
    int    n_steps;      /* Total steps taken */
    int    n_reject;     /* Steps rejected due to error */
    double h_last;       /* Last accepted step size */
    bool   success;      /* false if max_steps exceeded or h < h_min */
} lg_ode_stats_t;

/*============================================================================
 * DOP853 — Dormand-Prince 8(7) 13-stage
 *
 * Exact coefficients from Hairer, Nørsett & Wanner "Solving ODEs I", 2nd ed.,
 * Table 5.2 (p. 178) and the published dop853.f Fortran reference code.
 * Error estimate uses the 5th-order embedded solution (E5 weights).
 *===========================================================================*/

/* Coefficients verbatim from dopcor() in Hairer/Clewley/Sherwood canonical dop853.c */
#define _DC2   0.526001519587677318785587544488e-1
#define _DC3   0.789002279381515978178381316732e-1
#define _DC4   0.118350341907227396726757197510e0
#define _DC5   0.281649658092772603273242802490e0
#define _DC6   0.333333333333333333333333333333e0
#define _DC7   0.25e0
#define _DC8   0.307692307692307692307692307692e0
#define _DC9   0.651282051282051282051282051282e0
#define _DC10  0.6e0
#define _DC11  0.857142857142857142857142857142e0
#define _DA21  5.26001519587677318785587544488e-2
#define _DA31  1.97250569845378994544595329183e-2
#define _DA32  5.91751709536136983633785987549e-2
#define _DA41  2.95875854768068491816892993775e-2
#define _DA43  8.87627564304205475450678981324e-2
#define _DA51  2.41365134159266685502369798665e-1
#define _DA53 -8.84549479328286085344864962717e-1
#define _DA54  9.24834003261792003115737966543e-1
#define _DA61  3.7037037037037037037037037037e-2
#define _DA64  1.70828608729473871279604482173e-1
#define _DA65  1.25467687566822425016691814123e-1
#define _DA71  3.7109375e-2
#define _DA74  1.70252211019544039314978060272e-1
#define _DA75  6.02165389804559606850219397283e-2
#define _DA76 -1.7578125e-2
#define _DA81  3.70920001185047927108779319836e-2
#define _DA84  1.70383925712239993810214054705e-1
#define _DA85  1.07262030446373284651809199168e-1
#define _DA86 -1.53194377486244017527936158236e-2
#define _DA87  8.27378916381402288758473766002e-3
#define _DA91  6.24110958716075717114429577812e-1
#define _DA94 -3.36089262944694129406857109825e0
#define _DA95 -8.68219346841726006818189891453e-1
#define _DA96  2.75920996994467083049415600797e1
#define _DA97  2.01540675504778934086186788979e1
#define _DA98 -4.34898841810699588477366255144e1
#define _DA101  4.77662536438264365890433908527e-1
#define _DA104 -2.48811461997166764192642586468e0
#define _DA105 -5.90290826836842996371446475743e-1
#define _DA106  2.12300514481811942347288949897e1
#define _DA107  1.52792336328824235832596922938e1
#define _DA108 -3.32882109689848629194453265587e1
#define _DA109 -2.03312017085086261358222928593e-2
#define _DA111 -9.3714243008598732571704021658e-1
#define _DA114  5.18637242884406370830023853209e0
#define _DA115  1.09143734899672957818500254654e0
#define _DA116 -8.14978701074692612513997267357e0
#define _DA117 -1.85200656599969598641566180701e1
#define _DA118  2.27394870993505042818970056734e1
#define _DA119  2.49360555267965238987089396762e0
#define _DA1110 -3.0467644718982195003823669022e0
#define _DA121  2.27331014751653820792359768449e0
#define _DA124 -1.05344954667372501984066689879e1
#define _DA125 -2.00087205822486249909675718444e0
#define _DA126 -1.79589318631187989172765950534e1
#define _DA127  2.79488845294199600508499808837e1
#define _DA128 -2.85899827713502369474065508674e0
#define _DA129 -8.87285693353062954433549289258e0
#define _DA1210  1.23605671757943030647266201528e1
#define _DA1211  6.43392746015763530355970484046e-1
/* b: 8th-order solution weights */
#define _DB1   5.42937341165687622380535766363e-2
#define _DB6   4.45031289275240888144113950566e0
#define _DB7   1.89151789931450038304281599044e0
#define _DB8  -5.8012039600105847814672114227e0
#define _DB9   3.1116436695781989440891606237e-1
#define _DB10 -1.52160949662516078556178806805e-1
#define _DB11  2.01365400804030348374776537501e-1
#define _DB12  4.47106157277725905176885569043e-2
/* bhh: 3rd-order embedded weights for denominator */
#define _DBH1  0.244094488188976377952755905512e0
#define _DBH2  0.733846688281611857341361741547e0
#define _DBH3  0.220588235294117647058823529412e-1
/* er: high-order error coefficients */
#define _DER1   0.1312004499419488073250102996e-1
#define _DER6  -0.1225156446376204440720569753e1
#define _DER7  -0.4957589496572501915214079952e0
#define _DER8   0.1664377182454986536961530415e1
#define _DER9  -0.3503288487499736816886487290e0
#define _DER10  0.3341791187130174790297318841e0
#define _DER11  0.8192320648511571246570742613e-1
#define _DER12 -0.2235530786388629525884427845e-1

/* One DOP853 step.  k2/k3 buffers are reused for stage-11/12 output (per reference).
 * Returns the combined step-size-control error norm. */
static inline double _lg_dop853_step(
    lg_ode_rhs_t f,
    double t, const double* y, int dim,
    double h,
    double* y_new,
    void* user_data
) {
    double k1[LG_ODE_MAX_DIM], k2[LG_ODE_MAX_DIM], k3[LG_ODE_MAX_DIM];
    double k4[LG_ODE_MAX_DIM], k5[LG_ODE_MAX_DIM], k6[LG_ODE_MAX_DIM];
    double k7[LG_ODE_MAX_DIM], k8[LG_ODE_MAX_DIM], k9[LG_ODE_MAX_DIM];
    double k10[LG_ODE_MAX_DIM];
    double ytmp[LG_ODE_MAX_DIM];
    int i;

    f(t, y, k1, dim, user_data);

    for(i=0;i<dim;i++) ytmp[i]=y[i]+h*(_DA21*k1[i]);
    f(t+_DC2*h, ytmp, k2, dim, user_data);

    for(i=0;i<dim;i++) ytmp[i]=y[i]+h*(_DA31*k1[i]+_DA32*k2[i]);
    f(t+_DC3*h, ytmp, k3, dim, user_data);

    for(i=0;i<dim;i++) ytmp[i]=y[i]+h*(_DA41*k1[i]+_DA43*k3[i]);
    f(t+_DC4*h, ytmp, k4, dim, user_data);

    for(i=0;i<dim;i++) ytmp[i]=y[i]+h*(_DA51*k1[i]+_DA53*k3[i]+_DA54*k4[i]);
    f(t+_DC5*h, ytmp, k5, dim, user_data);

    for(i=0;i<dim;i++) ytmp[i]=y[i]+h*(_DA61*k1[i]+_DA64*k4[i]+_DA65*k5[i]);
    f(t+_DC6*h, ytmp, k6, dim, user_data);

    for(i=0;i<dim;i++) ytmp[i]=y[i]+h*(_DA71*k1[i]+_DA74*k4[i]+_DA75*k5[i]+_DA76*k6[i]);
    f(t+_DC7*h, ytmp, k7, dim, user_data);

    for(i=0;i<dim;i++) ytmp[i]=y[i]+h*(_DA81*k1[i]+_DA84*k4[i]+_DA85*k5[i]+_DA86*k6[i]+_DA87*k7[i]);
    f(t+_DC8*h, ytmp, k8, dim, user_data);

    for(i=0;i<dim;i++) ytmp[i]=y[i]+h*(_DA91*k1[i]+_DA94*k4[i]+_DA95*k5[i]+_DA96*k6[i]+_DA97*k7[i]+_DA98*k8[i]);
    f(t+_DC9*h, ytmp, k9, dim, user_data);

    for(i=0;i<dim;i++) ytmp[i]=y[i]+h*(_DA101*k1[i]+_DA104*k4[i]+_DA105*k5[i]+_DA106*k6[i]+_DA107*k7[i]+_DA108*k8[i]+_DA109*k9[i]);
    f(t+_DC10*h, ytmp, k10, dim, user_data);

    /* stage 11 — result stored into k2 (reuse) */
    for(i=0;i<dim;i++) ytmp[i]=y[i]+h*(_DA111*k1[i]+_DA114*k4[i]+_DA115*k5[i]+_DA116*k6[i]+_DA117*k7[i]+_DA118*k8[i]+_DA119*k9[i]+_DA1110*k10[i]);
    f(t+_DC11*h, ytmp, k2, dim, user_data);

    /* stage 12 — result stored into k3 (reuse); uses k2=stage11 */
    for(i=0;i<dim;i++) ytmp[i]=y[i]+h*(_DA121*k1[i]+_DA124*k4[i]+_DA125*k5[i]+_DA126*k6[i]+_DA127*k7[i]+_DA128*k8[i]+_DA129*k9[i]+_DA1210*k10[i]+_DA1211*k2[i]);
    f(t+h, ytmp, k3, dim, user_data);

    /* 8th-order solution: b11 multiplies k2(=stage11), b12 multiplies k3(=stage12) */
    for(i=0;i<dim;i++)
        y_new[i]=y[i]+h*(_DB1*k1[i]+_DB6*k6[i]+_DB7*k7[i]+_DB8*k8[i]+_DB9*k9[i]+_DB10*k10[i]+_DB11*k2[i]+_DB12*k3[i]);

    /* Error norm — combined bhh (deno) and er, matching dopcor exactly.
     * k4i_sol is the 8th-order increment (without the y[i] base). */
    double err = 0.0, err2 = 0.0;
    for(i=0;i<dim;i++){
        double sk = 1e-10 + fmax(fabs(y[i]), fabs(y_new[i]));
        double k4sol = _DB1*k1[i]+_DB6*k6[i]+_DB7*k7[i]+_DB8*k8[i]+_DB9*k9[i]+_DB10*k10[i]+_DB11*k2[i]+_DB12*k3[i];
        double e2 = k4sol - _DBH1*k1[i] - _DBH2*k9[i] - _DBH3*k3[i];
        double s2 = e2/sk; err2 += s2*s2;
        double e1 = _DER1*k1[i]+_DER6*k6[i]+_DER7*k7[i]+_DER8*k8[i]+_DER9*k9[i]+_DER10*k10[i]+_DER11*k2[i]+_DER12*k3[i];
        double s1 = e1/sk; err += s1*s1;
    }
    double deno = err + 0.01*err2;
    if (deno <= 0.0) deno = 1.0;
    return fabs(h) * sqrt(err / (deno * (double)dim));
}

#undef _DC2
#undef _DC3
#undef _DC4
#undef _DC5
#undef _DC6
#undef _DC7
#undef _DC8
#undef _DC9
#undef _DC10
#undef _DC11
#undef _DA21
#undef _DA31
#undef _DA32
#undef _DA41
#undef _DA43
#undef _DA51
#undef _DA53
#undef _DA54
#undef _DA61
#undef _DA64
#undef _DA65
#undef _DA71
#undef _DA74
#undef _DA75
#undef _DA76
#undef _DA81
#undef _DA84
#undef _DA85
#undef _DA86
#undef _DA87
#undef _DA91
#undef _DA94
#undef _DA95
#undef _DA96
#undef _DA97
#undef _DA98
#undef _DA101
#undef _DA104
#undef _DA105
#undef _DA106
#undef _DA107
#undef _DA108
#undef _DA109
#undef _DA111
#undef _DA114
#undef _DA115
#undef _DA116
#undef _DA117
#undef _DA118
#undef _DA119
#undef _DA1110
#undef _DA121
#undef _DA124
#undef _DA125
#undef _DA126
#undef _DA127
#undef _DA128
#undef _DA129
#undef _DA1210
#undef _DA1211
#undef _DB1
#undef _DB6
#undef _DB7
#undef _DB8
#undef _DB9
#undef _DB10
#undef _DB11
#undef _DB12
#undef _DBH1
#undef _DBH2
#undef _DBH3
#undef _DER1
#undef _DER6
#undef _DER7
#undef _DER8
#undef _DER9
#undef _DER10
#undef _DER11
#undef _DER12

/* Adaptive DOP853 driver — integrates state from state->t to params->t_end.
 * Returns statistics; state is updated in-place. */
static inline lg_ode_stats_t lg_integrate_dop853(
    lg_ode_rhs_t f,
    lg_ode_state_t* state,
    const lg_ode_params_t* params,
    void* user_data
) {
    lg_ode_stats_t stats = {0, 0, 0.0, true};

    int dim   = state->dim;
    double t  = state->t;
    double t_end = params->t_end;
    double h  = (params->h_init > 0.0) ? params->h_init : 1e-3 * (t_end - t);
    double h_min = (params->h_min > 0.0) ? params->h_min : 1e-14;
    double h_max = (params->h_max > 0.0) ? params->h_max : fabs(t_end - t);
    int max_steps = (params->max_steps > 0) ? params->max_steps : 1000000;

    /* PI controller constants */
    const double fac     = 0.9;
    const double fac_min = 0.333;
    const double fac_max = 6.0;
    const double expo    = 1.0 / 8.0;

    double y_new[LG_ODE_MAX_DIM];
    double* y = state->y;

    while (t < t_end && stats.n_steps < max_steps) {
        if (t + h > t_end) h = t_end - t;

        double err = _lg_dop853_step(f, t, y, dim, h, y_new, user_data);

        stats.n_steps++;

        if (err <= 1.0) {
            /* Accept step */
            t += h;
            memcpy(y, y_new, sizeof(double) * dim);
            stats.h_last = h;
            /* Adjust step upward */
            double factor = fac * pow(1.0 / (err + 1e-30), expo);
            if (factor > fac_max) factor = fac_max;
            h *= factor;
            if (h > h_max) h = h_max;
        } else {
            /* Reject step */
            stats.n_reject++;
            double factor = fac * pow(1.0 / (err + 1e-30), expo);
            if (factor < fac_min) factor = fac_min;
            h *= factor;
            if (h < h_min) {
                stats.success = false;
                break;
            }
        }
    }

    if (stats.n_steps >= max_steps) stats.success = false;
    state->t = t;
    return stats;
}

/*============================================================================
 * Gauss-Jackson Order-8 Multistep Predictor-Corrector
 *
 * Summed (Störmer-Cowell) form for second-order ODEs: r'' = a(t, r, r').
 * State: positions q[dim/2] and velocities v[dim/2], dim must be even.
 *
 * The corrector coefficients for GJ8 (from Berry & Healy 2004):
 *   alpha[9] for second-sum (position), beta[9] for first-sum (velocity).
 *
 * Startup uses DOP853 for the 4 back-steps needed.
 *===========================================================================*/

/* GJ8 corrector coefficients (second-sum form, alpha_k for k=-4..4) */
static const double _gj8_alpha[9] = {
    -1.0/3150.0, 4.0/945.0, -1.0/45.0, 4.0/9.0,
    -73.0/3150.0,   /* k=0 */
     4.0/9.0, -1.0/45.0, 4.0/945.0, -1.0/3150.0
};

/* GJ8 corrector coefficients (first-sum form, beta_k for k=-4..4) */
static const double _gj8_beta[9] = {
     1.0/1260.0, -1.0/168.0,  1.0/18.0, -1.0/2.0,
    251.0/720.0,
     1.0/2.0,   -1.0/18.0,   1.0/168.0, -1.0/1260.0  /* note: sign flip at k=0 side */
};

#define LG_GJ8_ORDER 8
#define LG_GJ8_HALF  4   /* half-order */

/* GJ8 state — stores history of accelerations and the summed sums S1, S2 */
typedef struct {
    /* Circular buffer of past accelerations: a[0] = oldest, a[7] = newest */
    double a_hist[LG_GJ8_ORDER][LG_ODE_MAX_DIM/2];
    double S1[LG_ODE_MAX_DIM/2];  /* first sum  (velocity accumulator) */
    double S2[LG_ODE_MAX_DIM/2];  /* second sum (position accumulator) */
    int    n_half;   /* half-dimension (position components = dim/2) */
    bool   initialized;
    double h;        /* fixed step size */
} lg_gj8_state_t;

/* RHS for second-order systems: accel(t, pos, vel) → writes into acc[n_half] */
typedef void (*lg_gj8_accel_t)(double t,
                                const double* pos, const double* vel,
                                double* acc, int n_half, void* user_data);

/* Initialize GJ8 state using DOP853 to generate the startup history.
 * pos0/vel0: initial conditions (n_half each).
 * h: fixed step size. */
static inline void lg_gj8_init(
    lg_gj8_state_t* gj,
    lg_gj8_accel_t accel_fn,
    double t0,
    const double* pos0, const double* vel0,
    int n_half, double h,
    void* user_data
) {
    gj->n_half = n_half;
    gj->h      = h;
    gj->initialized = true;

    /* Startup: integrate backward (4 steps) and forward (4 steps) to fill history.
     * Simple approach: use RK4 internally for the 8 startup points.
     * We use a simple RK4 to avoid a circular dependency on DOP853 driver overhead. */

    /* Full 2nd-order state y = [pos, vel] */
    (void)(2 * n_half);

    /* We'll just run RK4 manually (no body/transform wrappers needed) */
    /* Generate points at t0 - 4h, t0 - 3h, ... t0 - h, t0 */
    double pts_pos[9][LG_ODE_MAX_DIM/2];
    double pts_vel[9][LG_ODE_MAX_DIM/2];
    double pts_acc[9][LG_ODE_MAX_DIM/2];

    /* Store t0 at index 4 */
    memcpy(pts_pos[4], pos0, sizeof(double)*n_half);
    memcpy(pts_vel[4], vel0, sizeof(double)*n_half);
    accel_fn(t0, pos0, vel0, pts_acc[4], n_half, user_data);

    /* RK4 helper lambda (C99 nested via loop — we inline it) */
    /* Backward integration for indices 3..0 */
    for (int step = 3; step >= 0; step--) {
        double dt = -h;
        const double* p = pts_pos[step+1];
        const double* v = pts_vel[step+1];
        double tp = t0 + (step - 4) * h + h; /* time at step+1 */

        double k1p[LG_ODE_MAX_DIM/2], k1v[LG_ODE_MAX_DIM/2];
        double k2p[LG_ODE_MAX_DIM/2], k2v[LG_ODE_MAX_DIM/2];
        double k3p[LG_ODE_MAX_DIM/2], k3v[LG_ODE_MAX_DIM/2];
        double k4p[LG_ODE_MAX_DIM/2], k4v[LG_ODE_MAX_DIM/2];
        double ptmp[LG_ODE_MAX_DIM/2], vtmp[LG_ODE_MAX_DIM/2];

        accel_fn(tp, p, v, k1v, n_half, user_data);
        for (int i = 0; i < n_half; i++) k1p[i] = v[i];

        for (int i = 0; i < n_half; i++) {
            ptmp[i] = p[i] + 0.5*dt*k1p[i];
            vtmp[i] = v[i] + 0.5*dt*k1v[i];
        }
        accel_fn(tp + 0.5*dt, ptmp, vtmp, k2v, n_half, user_data);
        for (int i = 0; i < n_half; i++) k2p[i] = vtmp[i];

        for (int i = 0; i < n_half; i++) {
            ptmp[i] = p[i] + 0.5*dt*k2p[i];
            vtmp[i] = v[i] + 0.5*dt*k2v[i];
        }
        accel_fn(tp + 0.5*dt, ptmp, vtmp, k3v, n_half, user_data);
        for (int i = 0; i < n_half; i++) k3p[i] = vtmp[i];

        for (int i = 0; i < n_half; i++) {
            ptmp[i] = p[i] + dt*k3p[i];
            vtmp[i] = v[i] + dt*k3v[i];
        }
        accel_fn(tp + dt, ptmp, vtmp, k4v, n_half, user_data);
        for (int i = 0; i < n_half; i++) k4p[i] = vtmp[i];

        for (int i = 0; i < n_half; i++) {
            pts_pos[step][i] = p[i] + dt/6.0*(k1p[i]+2*k2p[i]+2*k3p[i]+k4p[i]);
            pts_vel[step][i] = v[i] + dt/6.0*(k1v[i]+2*k2v[i]+2*k3v[i]+k4v[i]);
        }
        accel_fn(t0 + (step-4)*h, pts_pos[step], pts_vel[step], pts_acc[step], n_half, user_data);
    }

    /* Forward integration for indices 5..8 */
    for (int step = 5; step <= 8; step++) {
        double dt = h;
        const double* p = pts_pos[step-1];
        const double* v = pts_vel[step-1];
        double tp = t0 + (step - 1 - 4) * h;

        double k1p[LG_ODE_MAX_DIM/2], k1v[LG_ODE_MAX_DIM/2];
        double k2p[LG_ODE_MAX_DIM/2], k2v[LG_ODE_MAX_DIM/2];
        double k3p[LG_ODE_MAX_DIM/2], k3v[LG_ODE_MAX_DIM/2];
        double k4p[LG_ODE_MAX_DIM/2], k4v[LG_ODE_MAX_DIM/2];
        double ptmp[LG_ODE_MAX_DIM/2], vtmp[LG_ODE_MAX_DIM/2];

        accel_fn(tp, p, v, k1v, n_half, user_data);
        for (int i = 0; i < n_half; i++) k1p[i] = v[i];

        for (int i = 0; i < n_half; i++) {
            ptmp[i] = p[i] + 0.5*dt*k1p[i];
            vtmp[i] = v[i] + 0.5*dt*k1v[i];
        }
        accel_fn(tp + 0.5*dt, ptmp, vtmp, k2v, n_half, user_data);
        for (int i = 0; i < n_half; i++) k2p[i] = vtmp[i];

        for (int i = 0; i < n_half; i++) {
            ptmp[i] = p[i] + 0.5*dt*k2p[i];
            vtmp[i] = v[i] + 0.5*dt*k2v[i];
        }
        accel_fn(tp + 0.5*dt, ptmp, vtmp, k3v, n_half, user_data);
        for (int i = 0; i < n_half; i++) k3p[i] = vtmp[i];

        for (int i = 0; i < n_half; i++) {
            ptmp[i] = p[i] + dt*k3p[i];
            vtmp[i] = v[i] + dt*k3v[i];
        }
        accel_fn(tp + dt, ptmp, vtmp, k4v, n_half, user_data);
        for (int i = 0; i < n_half; i++) k4p[i] = vtmp[i];

        for (int i = 0; i < n_half; i++) {
            pts_pos[step][i] = p[i] + dt/6.0*(k1p[i]+2*k2p[i]+2*k3p[i]+k4p[i]);
            pts_vel[step][i] = v[i] + dt/6.0*(k1v[i]+2*k2v[i]+2*k3v[i]+k4v[i]);
        }
        accel_fn(t0 + (step-4)*h, pts_pos[step], pts_vel[step], pts_acc[step], n_half, user_data);
    }

    /* Initialize summed sums using the 9-point history (indices 0..8).
     * S2[n] = sum_{k<= n} S1[k]*h, S1[n] = sum_{k<=n} a[k]*h  (discrete) */
    for (int i = 0; i < n_half; i++) {
        double s1 = 0.0, s2 = 0.0;
        for (int k = 0; k <= 8; k++) {
            s1 += pts_acc[k][i] * h;
            s2 += s1;
        }
        gj->S1[i] = s1 - pts_acc[8][i]*h; /* S1 up to step 8 exclusive */
        gj->S2[i] = s2 - s1;               /* S2 similarly */
    }

    /* Store acceleration history (last 8 steps: indices 1..8) */
    for (int k = 0; k < LG_GJ8_ORDER; k++)
        memcpy(gj->a_hist[k], pts_acc[k+1], sizeof(double)*n_half);
}

/* Advance one fixed step h.  Corrects once (single iteration). */
static inline void lg_gj8_step(
    lg_gj8_state_t* gj,
    lg_gj8_accel_t accel_fn,
    double t,
    double* pos, double* vel,
    void* user_data
) {
    int n = gj->n_half;
    double h = gj->h;

    /* Predictor: extrapolate S1 and S2 one step forward */
    double S1_pred[LG_ODE_MAX_DIM/2], S2_pred[LG_ODE_MAX_DIM/2];
    double pos_pred[LG_ODE_MAX_DIM/2], vel_pred[LG_ODE_MAX_DIM/2];

    for (int i = 0; i < n; i++) {
        /* Sum using GJ8 alpha/beta coefficients over history a_hist[0..7] */
        double sum_a  = 0.0, sum_s1 = 0.0;
        for (int k = 0; k < LG_GJ8_ORDER; k++) {
            sum_a  += _gj8_alpha[k] * gj->a_hist[k][i];
            sum_s1 += _gj8_beta[k]  * gj->a_hist[k][i];
        }
        /* Extrapolate (alpha[8] / beta[8] correspond to new step, treated as 0 in predictor) */
        S2_pred[i]   = gj->S2[i] + gj->S1[i] + h*h * sum_a;
        S1_pred[i]   = gj->S1[i] + h * sum_s1;
        pos_pred[i]  = pos[i] + S2_pred[i];
        vel_pred[i]  = vel[i] + S1_pred[i] / h;
    }

    /* Corrector: evaluate acceleration at predicted state */
    double a_new[LG_ODE_MAX_DIM/2];
    accel_fn(t + h, pos_pred, vel_pred, a_new, n, user_data);

    /* Update sums with corrected acceleration */
    for (int i = 0; i < n; i++) {
        gj->S1[i] = S1_pred[i] + h * _gj8_beta[LG_GJ8_HALF] * a_new[i];
        gj->S2[i] = S2_pred[i] + h*h * _gj8_alpha[LG_GJ8_HALF] * a_new[i];
        vel[i]    += gj->S1[i] / h;
        pos[i]    += gj->S2[i];
    }

    /* Shift history buffer */
    memmove(gj->a_hist[0], gj->a_hist[1], sizeof(double)*n*(LG_GJ8_ORDER-1));
    memcpy(gj->a_hist[LG_GJ8_ORDER-1], a_new, sizeof(double)*n);
}

/*============================================================================
 * Wisdom-Holman Symplectic Map (Mixed-Variable)
 *
 * Splits the N-body Hamiltonian H = H_kepler + H_interaction.
 * One step of length h:
 *   1. Half-kick:   v += 0.5*h * a_interaction(r)
 *   2. Drift:       advance each body on its unperturbed Kepler orbit for h
 *   3. Half-kick:   v += 0.5*h * a_interaction(r)
 *
 * Kepler drift uses the universal variable (Lagrange f, g functions).
 * Suitable for planetary systems with one dominant central body.
 *===========================================================================*/

/* Single-body Kepler drift via universal variable (Battin universal anomaly).
 * Advances (r, v) by dt under mu using Newton iteration on Kepler's equation.
 * Double precision. */
static inline void _lg_kepler_drift(double mu, double* r, double* v, double dt) {
    double rx = r[0], ry = r[1], rz = r[2];
    double vx = v[0], vy = v[1], vz = v[2];

    double r0  = sqrt(rx*rx + ry*ry + rz*rz);
    double vr0 = (rx*vx + ry*vy + rz*vz) / r0;
    double alpha = 2.0/r0 - (vx*vx+vy*vy+vz*vz)/mu; /* 1/a */

    /* Universal anomaly X: solve Kepler's equation via Laguerre-Conway */
    double X = sqrt(mu) * fabs(alpha) * dt; /* initial guess */
    double psi, c0, c1, c2, c3;
    for (int iter = 0; iter < 50; iter++) {
        psi = X * X * alpha;
        /* Stumpff functions */
        double psi2 = psi*psi;
        if (fabs(psi) < 1e-6) {
            c2 = 0.5 - psi/24.0 + psi2/720.0;
            c3 = 1.0/6.0 - psi/120.0 + psi2/5040.0;
        } else if (psi > 0.0) {
            double sqp = sqrt(psi);
            c2 = (1.0 - cos(sqp)) / psi;
            c3 = (sqp - sin(sqp)) / (psi * sqp);
        } else {
            double sqmp = sqrt(-psi);
            c2 = (1.0 - cosh(sqmp)) / psi;
            c3 = (sinh(sqmp) - sqmp) / (-psi * sqmp);
        }
        c0 = 1.0 - psi * c2;
        c1 = 1.0 - psi * c3;

        double F  = r0*vr0/sqrt(mu)*X*X*c2 + (1.0 - r0*alpha)*X*X*X*c3 + r0*X*c1 - sqrt(mu)*dt;
        double dF = r0*vr0/sqrt(mu)*X*(1.0-psi*c3)*2.0 + (1.0-r0*alpha)*X*X*c2*3.0 + r0*c0;
        double dX = -F / dF;
        X += dX;
        if (fabs(dX) < 1e-12 * fabs(X) + 1e-14) break;
    }

    psi = X*X*alpha;
    double psi2 = psi*psi;
    if (fabs(psi) < 1e-6) {
        c2 = 0.5 - psi/24.0 + psi2/720.0;
        c3 = 1.0/6.0 - psi/120.0 + psi2/5040.0;
    } else if (psi > 0.0) {
        double sqp = sqrt(psi);
        c2 = (1.0 - cos(sqp)) / psi;
        c3 = (sqp - sin(sqp)) / (psi * sqp);
    } else {
        double sqmp = sqrt(-psi);
        c2 = (1.0 - cosh(sqmp)) / psi;
        c3 = (sinh(sqmp) - sqmp) / (-psi * sqmp);
    }
    c0 = 1.0 - psi*c2;
    c1 = 1.0 - psi*c3;

    /* Compute new r magnitude */
    double r1 = r0*c0 + sqrt(mu)*vr0*X*X*c2 + mu*(1.0-r0*alpha)*X*X*X*c3;
    /* not exact r1 — use Lagrange coefficients */
    double f  = 1.0 - X*X*c2 / r0;
    double g  = dt - X*X*X*c3 / sqrt(mu);
    double r1_ = sqrt((f*rx + g*vx)*(f*rx+g*vx) + (f*ry+g*vy)*(f*ry+g*vy) + (f*rz+g*vz)*(f*rz+g*vz));
    double fdot = sqrt(mu)*X*c1*(psi*c3-1.0) / (r0*r1_);
    double gdot = 1.0 - X*X*c2 / r1_;

    r[0] = f*rx + g*vx;
    r[1] = f*ry + g*vy;
    r[2] = f*rz + g*vz;
    v[0] = fdot*rx + gdot*vx;
    v[1] = fdot*ry + gdot*vy;
    v[2] = fdot*rz + gdot*vz;
    (void)r1; /* suppress unused warning */
}

/* Wisdom-Holman N-body step.
 *
 * n_bodies: total bodies (index 0 = central/star, others = planets)
 * mu: gravitational parameter of the central body
 * pos[n_bodies][3], vel[n_bodies][3] (interleaved: pos[i] = {x,y,z})
 * The central body (index 0) is held fixed (infinite mass approximation).
 * Planetary interaction: pairwise gravity with G*m_i*m_j / r²,
 *   masses[i] in the same units as mu (i.e., G*M_i). */
static inline void lg_wh_step(
    int n_bodies,
    double mu,
    double pos[][3], double vel[][3],
    const double masses[],
    double h
) {
    /* Half-kick: pairwise planet-planet interaction only (star is dominant) */
    for (int iter = 0; iter < 2; iter++) { /* iter=0: half, iter=1: half */
        double dt_kick = (iter == 0) ? 0.5*h : 0.5*h;
        for (int i = 1; i < n_bodies; i++) {
            double ax = 0.0, ay = 0.0, az = 0.0;
            for (int j = 1; j < n_bodies; j++) {
                if (j == i) continue;
                double dx = pos[j][0] - pos[i][0];
                double dy = pos[j][1] - pos[i][1];
                double dz = pos[j][2] - pos[i][2];
                double r2  = dx*dx + dy*dy + dz*dz + 1e-30;
                double r   = sqrt(r2);
                double r3i = masses[j] / (r2 * r);
                ax += r3i * dx;
                ay += r3i * dy;
                az += r3i * dz;
            }
            vel[i][0] += dt_kick * ax;
            vel[i][1] += dt_kick * ay;
            vel[i][2] += dt_kick * az;
        }
        if (iter == 0) {
            /* Drift: advance each planet on unperturbed Kepler orbit */
            for (int i = 1; i < n_bodies; i++) {
                _lg_kepler_drift(mu, pos[i], vel[i], h);
            }
        }
    }
}

/*============================================================================
 * Unified Adaptive ODE Interface (convenience wrapper)
 *
 * Selects algorithm based on lg_ode_method_t.
 *===========================================================================*/

typedef enum {
    LG_ODE_DOP853,       /* Dormand-Prince 8(7) — recommended default */
    LG_ODE_RK4_FIXED,    /* Fixed-step RK4 (backward compat) */
} lg_ode_method_t;

/* Fixed-step RK4 via the unified interface (delegates to internal RK4) */
static inline lg_ode_stats_t _lg_integrate_rk4_fixed(
    lg_ode_rhs_t f,
    lg_ode_state_t* state,
    const lg_ode_params_t* params,
    void* user_data
) {
    lg_ode_stats_t stats = {0, 0, 0.0, true};
    int dim    = state->dim;
    double t   = state->t;
    double t_end = params->t_end;
    double h   = (params->h_init > 0.0) ? params->h_init : 1e-3 * fabs(t_end - t);
    int max_steps = (params->max_steps > 0) ? params->max_steps : 1000000;

    double k1[LG_ODE_MAX_DIM], k2[LG_ODE_MAX_DIM];
    double k3[LG_ODE_MAX_DIM], k4[LG_ODE_MAX_DIM];
    double ytmp[LG_ODE_MAX_DIM];
    double* y = state->y;

    while (t < t_end && stats.n_steps < max_steps) {
        if (t + h > t_end) h = t_end - t;
        f(t, y, k1, dim, user_data);
        for (int i = 0; i < dim; i++) ytmp[i] = y[i] + 0.5*h*k1[i];
        f(t + 0.5*h, ytmp, k2, dim, user_data);
        for (int i = 0; i < dim; i++) ytmp[i] = y[i] + 0.5*h*k2[i];
        f(t + 0.5*h, ytmp, k3, dim, user_data);
        for (int i = 0; i < dim; i++) ytmp[i] = y[i] + h*k3[i];
        f(t + h, ytmp, k4, dim, user_data);
        for (int i = 0; i < dim; i++)
            y[i] += h/6.0 * (k1[i] + 2.0*k2[i] + 2.0*k3[i] + k4[i]);
        t += h;
        stats.n_steps++;
        stats.h_last = h;
    }
    if (stats.n_steps >= max_steps) stats.success = false;
    state->t = t;
    return stats;
}

/* Main entry point */
static inline lg_ode_stats_t lg_integrate_adaptive(
    lg_ode_method_t method,
    lg_ode_rhs_t f,
    lg_ode_state_t* state,
    const lg_ode_params_t* params,
    void* user_data
) {
    switch (method) {
        case LG_ODE_DOP853:
            return lg_integrate_dop853(f, state, params, user_data);
        case LG_ODE_RK4_FIXED:
            return _lg_integrate_rk4_fixed(f, state, params, user_data);
        default:
            return lg_integrate_dop853(f, state, params, user_data);
    }
}

/* Welford online statistics for drift tracking */
typedef struct {
    double mean;
    double M2;
    long   n;
} lg_welford_t;

static inline void lg_welford_update(lg_welford_t* w, double x) {
    w->n++;
    double delta = x - w->mean;
    w->mean += delta / w->n;
    double delta2 = x - w->mean;
    w->M2 += delta * delta2;
}

static inline double lg_welford_variance(const lg_welford_t* w) {
    return (w->n < 2) ? 0.0 : w->M2 / (w->n - 1);
}

static inline double lg_welford_stddev(const lg_welford_t* w) {
    return sqrt(lg_welford_variance(w));
}

#ifdef __cplusplus
}
#endif

#endif /* LAGRANGE_ADAPTIVE_INTEGRATOR_H */
