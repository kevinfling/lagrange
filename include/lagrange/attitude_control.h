/**
 * lagrange_attitude_control.h - Reaction Wheels, CMGs, and Gyroscopic Dynamics
 * 
 * Momentum exchange devices for spacecraft attitude control:
 *   - Reaction Wheels (RW): variable speed, fixed axis, 1 DOF control per wheel
 *   - Control Moment Gyros (CMG): constant speed, gimballed axis, 2 DOF torque amplification
 *   - Variable Speed CMG (VSCMG): combined (not yet common in flight)
 *   - Gyroscopic precession coupling in rigid body dynamics
 *   - Momentum management (saturation avoidance, external torque desaturation)
 * 
 * Integrates with: lagrange_body.h (angular momentum, torque application),
 *                  lagrange_integrator.h (gyroscopic terms in equations of motion),
 *                  lagrange_math.h (quaternion kinematics, rotation matrices)
 */

#ifndef LAGRANGE_ATTITUDE_CONTROL_H
#define LAGRANGE_ATTITUDE_CONTROL_H

#include "math.h"
#include "body.h"
#include "integrator.h"
#include <string.h>

#ifdef __cplusplus
extern "C" {
#endif

/*============================================================================
 * Gyroscopic Fundamentals (Euler's Equations with Internal Angular Momentum)
 * 
 * Standard rigid body: I*omega_dot + omega x (I*omega) = tau_ext
 * With reaction wheels: h_total = I*omega + A*h_wheels
 * Euler eq: I*omega_dot + omega x (I*omega + A*h_w) = tau_ext - A*tau_w
 * 
 * Where A is Jacobian (axis alignment), tau_w is wheel motor torque
 *===========================================================================*/

typedef struct {
    lg_vec3_t h_total;        /* Total system angular momentum (body + wheels) */
    lg_vec3_t h_body;         /* Body contribution I*omega */
    lg_vec3_t h_wheels;       /* Wheels contribution (sum of all devices) */
    float h_mag;              /* |h_total| for saturation monitoring */
} lg_gyro_state_t;

/* Compute gyroscopic torque: tau_gyro = -omega x h_wheels (cross-coupling) */
static inline lg_vec3_t lg_gyroscopic_coupling(const lg_vec3_t* omega,
                                              const lg_vec3_t* h_wheels) {
    return lg_vec3_cross(*omega, *h_wheels);
}

/* Modified Euler equation with internal angular momentum */
static inline lg_vec3_t lg_euler_gyroscopic(const lg_body_t* body,
                                           const lg_vec3_t* h_wheels,
                                           const lg_vec3_t* tau_ext) {
    /* tau = I*omega_dot + omega x (I*omega + h_w) */
    /* Solve for omega_dot: */
    
    lg_vec3_t I_omega = lg_vec3_mul(body->angular_velocity, body->inertia);
    lg_vec3_t h_total = lg_vec3_add(I_omega, *h_wheels);
    
    /* Gyroscopic term */
    lg_vec3_t gyro = lg_vec3_cross(body->angular_velocity, h_total);
    
    /* Net torque on body */
    lg_vec3_t net_tau = lg_vec3_sub(*tau_ext, gyro);
    
    /* Motor torques (from wheels) would be added here as -A*tau_motor */
    
    /* Return angular acceleration: inv(I) * net_tau */
    return lg_vec3_mul(net_tau, body->inv_inertia);
}

/*============================================================================
 * Reaction Wheel Model
 * 
 * Axis fixed in body frame, angular momentum h = J*omega_wheel
 * Torque on body = -axis * tau_motor (Newton's 3rd law)
 *===========================================================================*/

typedef struct {
    lg_vec3_t axis;           /* Spin axis in body frame (normalized) */
    float J;                  /* Moment of inertia about spin axis (kg*m^2) */
    float omega;            /* Current spin rate (rad/s) */
    float h;                /* Angular momentum J*omega */
    float max_omega;        /* Saturation speed (rad/s) */
    float max_torque;       /* Motor torque limit (Nm) */
    
    /* Control */
    float h_command;        /* Target momentum */
    float tau_command;      /* Current motor torque command */
    
    /* State */
    float momentum_saturation; /* 0..1, 1 = fully saturated */
} lg_reaction_wheel_t;

/* Initialize RW with common specs (e.g., 0.5 kg*m^2, 6000 RPM max) */
static inline lg_reaction_wheel_t lg_rw_init(lg_vec3_t axis_body, float J, float max_rpm) {
    lg_reaction_wheel_t rw;
    rw.axis = lg_vec3_norm(axis_body);
    rw.J = J;
    rw.omega = 0.0f;
    rw.h = 0.0f;
    rw.max_omega = max_rpm * 2.0f * LG_PI / 60.0f;
    rw.max_torque = 0.1f; /* Default 0.1 Nm */
    rw.h_command = 0.0f;
    rw.tau_command = 0.0f;
    rw.momentum_saturation = 0.0f;
    return rw;
}

/* Standard NASA specs: 34.3 kg wheel, 0.65 kg*m^2, 6000 RPM */
static inline lg_reaction_wheel_t lg_rw_nasa_std(lg_vec3_t axis) {
    return lg_rw_init(axis, 0.65f, 6000.0f);
}

/* Update wheel dynamics (spin acceleration) */
static inline void lg_rw_update(lg_reaction_wheel_t* rw, float dt) {
    /* Clamp command torque */
    float tau = fmaxf(-rw->max_torque, fminf(rw->tau_command, rw->max_torque));
    
    /* Integrate spin rate */
    float omega_dot = tau / rw->J;
    rw->omega += omega_dot * dt;
    
    /* Saturation limits */
    if (rw->omega > rw->max_omega) {
        rw->omega = rw->max_omega;
        omega_dot = 0.0f; /* Can't accelerate further */
    } else if (rw->omega < -rw->max_omega) {
        rw->omega = -rw->max_omega;
        omega_dot = 0.0f;
    }
    
    /* Update momentum */
    rw->h = rw->J * rw->omega;
    rw->momentum_saturation = fabsf(rw->h) / (rw->J * rw->max_omega);
}

/* Get angular momentum vector in body frame */
static inline lg_vec3_t lg_rw_angular_momentum(const lg_reaction_wheel_t* rw) {
    return lg_vec3_scale(rw->axis, rw->h);
}

/* Get torque on spacecraft body (negative of wheel acceleration torque) */
static inline lg_vec3_t lg_rw_body_torque(const lg_reaction_wheel_t* rw) {
    /* tau_body = -axis * tau_motor (reaction torque) */
    return lg_vec3_scale(rw->axis, -rw->tau_command);
}

/*============================================================================
 * Reaction Wheel Array (3+ wheels for 3-axis control)
 *===========================================================================*/

#define LG_MAX_WHEELS 6

typedef struct {
    lg_reaction_wheel_t wheels[LG_MAX_WHEELS];
    int n_wheels;
    float total_momentum_capacity;
    
    /* Control distribution matrix (pseudo-inverse of axis matrix) */
    float W[3][LG_MAX_WHEELS]; /* Weighting matrix for distribution */
    bool singularity_warn;       /* Near degenerate configuration */
} lg_rw_array_t;

/* Common configurations */
static inline lg_rw_array_t lg_rw_array_pyramid(void) {
    /* 4-wheel pyramid (NASA standard): each wheel at 45deg to base, 
     * base wheels at 90deg to each other */
    lg_rw_array_t arr;
    arr.n_wheels = 4;
    
    float angle = 45.0f * LG_PI / 180.0f;
    float base_angles[4] = {45, 135, 225, 315}; /* deg */
    
    for (int i = 0; i < 4; i++) {
        float phi = base_angles[i] * LG_PI / 180.0f;
        lg_vec3_t axis = {
            sinf(angle) * cosf(phi),
            sinf(angle) * sinf(phi),
            cosf(angle)
        };
        arr.wheels[i] = lg_rw_nasa_std(axis);
    }
    
    arr.total_momentum_capacity = 4.0f * arr.wheels[0].J * arr.wheels[0].max_omega;
    arr.singularity_warn = false;
    return arr;
}

/* Skewed 3-wheel configuration (minimum for full 3D) */
static inline lg_rw_array_t lg_rw_array_triad(void) {
    lg_rw_array_t arr;
    arr.n_wheels = 3;
    
    /* Mutually orthogonal but rotated to avoid singularity */
    arr.wheels[0] = lg_rw_nasa_std(lg_vec3(1, 0, 0));
    arr.wheels[1] = lg_rw_nasa_std(lg_vec3(0, 1, 0));
    arr.wheels[2] = lg_rw_nasa_std(lg_vec3(0, 0, 1));
    
    arr.total_momentum_capacity = 3.0f * arr.wheels[0].J * arr.wheels[0].max_omega;
    return arr;
}

/* Compute total stored momentum */
static inline lg_vec3_t lg_rw_array_momentum(const lg_rw_array_t* arr) {
    lg_vec3_t h = lg_vec3_zero();
    for (int i = 0; i < arr->n_wheels; i++) {
        h = lg_vec3_add(h, lg_rw_angular_momentum(&arr->wheels[i]));
    }
    return h;
}

/* Distribute torque command to wheels (simple pseudo-inverse) */
static inline void lg_rw_array_command_torque(lg_rw_array_t* arr,
                                             const lg_vec3_t* tau_cmd) {
    /* Build Jacobian A (3 x n) where column i is wheel i axis */
    /* tau_body = A * tau_motors */
    /* Minimize energy: tau_motors = A^T * (A*A^T)^{-1} * tau_cmd */

    float A[3][LG_MAX_WHEELS];
    for (int i = 0; i < arr->n_wheels; i++) {
        A[0][i] = arr->wheels[i].axis.x;
        A[1][i] = arr->wheels[i].axis.y;
        A[2][i] = arr->wheels[i].axis.z;
    }

    /* Compute A*A^T (3x3) */
    float AAT[3][3] = {0};
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            for (int k = 0; k < arr->n_wheels; k++) {
                AAT[i][j] += A[i][k] * A[j][k];
            }
        }
    }

    /* Convert to lg_mat3_t (column-major) and invert */
    lg_mat3_t m;
    m.m[0] = AAT[0][0]; m.m[1] = AAT[1][0]; m.m[2] = AAT[2][0];
    m.m[3] = AAT[0][1]; m.m[4] = AAT[1][1]; m.m[5] = AAT[2][1];
    m.m[6] = AAT[0][2]; m.m[7] = AAT[1][2]; m.m[8] = AAT[2][2];

    float det = lg_mat3_determinant(m);
    arr->singularity_warn = fabsf(det) < 1e-6f;

    lg_mat3_t inv = lg_mat3_inverse(m);
    lg_vec3_t temp = lg_mat3_mul_vec3(inv, *tau_cmd);

    /* Distribute: tau_motor = A^T * temp */
    for (int i = 0; i < arr->n_wheels; i++) {
        arr->wheels[i].tau_command = A[0][i]*temp.x + A[1][i]*temp.y + A[2][i]*temp.z;
    }
}

/* Update entire array */
static inline void lg_rw_array_update(lg_rw_array_t* arr, float dt) {
    for (int i = 0; i < arr->n_wheels; i++) {
        lg_rw_update(&arr->wheels[i], dt);
    }
}

/* Check for momentum saturation (all wheels near max speed) */
static inline bool lg_rw_array_saturated(const lg_rw_array_t* arr, float threshold) {
    for (int i = 0; i < arr->n_wheels; i++) {
        if (arr->wheels[i].momentum_saturation < threshold) return false;
    }
    return true;
}

/* Desaturation using external torque (magnetic torque coils, thrusters, gravity gradient) */
static inline void lg_rw_array_desaturate(lg_rw_array_t* arr,
                                         const lg_vec3_t* external_torque_available,
                                         float strength) {
    /* Apply external torque opposite to net wheel momentum */
    lg_vec3_t h_net = lg_rw_array_momentum(arr);
    (void)h_net;
    
    /* Command would mix with attitude control... simplified here */
}

/*============================================================================
 * Control Moment Gyro (CMG)
 * 
 * Constant high-speed rotor, gimballed axis. Torque = h x gimbal_rate
 * "Torque amplification": small gimbal torque -> large output torque
 *===========================================================================*/

typedef struct {
    lg_vec3_t spin_axis;      /* Current orientation of spin vector (body frame) */
    lg_vec3_t gimbal_axis;    /* Gimbal rotation axis (perpendicular to spin) */
    float h;                  /* Constant angular momentum magnitude (J*omega) */
    float gimbal_angle;       /* Current gimbal position (rad) */
    float gimbal_rate;        /* Current gimbal rate (rad/s) */
    float max_gimbal_rate;    /* Mechanical limit */
    float max_gimbal_angle;   /* Typically +/- 90 deg or less (hardware stop) */
    
    /* Control */
    float gimbal_rate_cmd;    /* Commanded gimbal rate */
    
    /* State */
    bool singular;            /* Near gimbal lock (rank deficiency) */
} lg_cmg_t;

/* Initialize CMG with spin axis along x, gimbal about z (can vary) */
static inline lg_cmg_t lg_cmg_init(lg_vec3_t initial_spin_axis, 
                                    lg_vec3_t gimbal_axis,
                                    float momentum,
                                    float max_rate_deg_s) {
    lg_cmg_t cmg;
    cmg.spin_axis = lg_vec3_norm(initial_spin_axis);
    cmg.gimbal_axis = lg_vec3_norm(gimbal_axis);
    cmg.h = momentum;
    cmg.gimbal_angle = 0.0f;
    cmg.gimbal_rate = 0.0f;
    cmg.max_gimbal_rate = max_rate_deg_s * LG_PI / 180.0f;
    cmg.max_gimbal_angle = LG_PI / 2.0f; /* +/- 90 deg */
    cmg.gimbal_rate_cmd = 0.0f;
    cmg.singular = false;
    return cmg;
}

/* High-power CMG (e.g., 5000 RPM, 100 kg rotor, 1m radius) */
static inline lg_cmg_t lg_cmg_high_power(lg_vec3_t spin, lg_vec3_t gimbal) {
    /* J = 0.5 * m * r^2 for disk, omega = 5000*2pi/60 */
    float J = 0.5f * 100.0f * 1.0f; /* 50 kg*m^2 */
    float omega = 5000.0f * 2.0f * LG_PI / 60.0f;
    return lg_cmg_init(spin, gimbal, J * omega, 10.0f); /* 10 deg/s gimbal */
}

/* Update gimbal kinematics */
static inline void lg_cmg_update(lg_cmg_t* cmg, float dt) {
    /* Clamp command */
    float rate = fmaxf(-cmg->max_gimbal_rate, 
                      fminf(cmg->gimbal_rate_cmd, cmg->max_gimbal_rate));
    
    /* Integrate angle */
    cmg->gimbal_angle += rate * dt;
    
    /* Hardware stops */
    if (cmg->gimbal_angle > cmg->max_gimbal_angle) {
        cmg->gimbal_angle = cmg->max_gimbal_angle;
        rate = 0.0f;
    } else if (cmg->gimbal_angle < -cmg->max_gimbal_angle) {
        cmg->gimbal_angle = -cmg->max_gimbal_angle;
        rate = 0.0f;
    }
    
    cmg->gimbal_rate = rate;
    
    /* Update spin axis orientation via rotation about gimbal axis */
    lg_quat_t delta_q = lg_quat_from_axis_angle(cmg->gimbal_axis, rate * dt);
    lg_quat_t spin_q = lg_quat(cmg->spin_axis.x, cmg->spin_axis.y, cmg->spin_axis.z, 0.0f);
    lg_quat_t new_spin_q = lg_quat_mul(lg_quat_mul(delta_q, spin_q), lg_quat_conj(delta_q));
    cmg->spin_axis = lg_vec3_norm(lg_vec3(new_spin_q.x, new_spin_q.y, new_spin_q.z));
}

/* Output torque on spacecraft: tau = h x gimbal_rate (cross product) */
/* Actually: tau = gimbal_rate * (gimbal_axis x h_vector) */
static inline lg_vec3_t lg_cmg_output_torque(const lg_cmg_t* cmg) {
    /* h_vector = spin_axis * h */
    lg_vec3_t h_vec = lg_vec3_scale(cmg->spin_axis, cmg->h);
    
    /* tau = gimbal_rate * (gimbal_axis x h_vec) */
    lg_vec3_t torque_dir = lg_vec3_cross(cmg->gimbal_axis, h_vec);
    return lg_vec3_scale(torque_dir, cmg->gimbal_rate);
}

/* Angular momentum stored in CMG (constant magnitude, varying direction) */
static inline lg_vec3_t lg_cmg_angular_momentum(const lg_cmg_t* cmg) {
    return lg_vec3_scale(cmg->spin_axis, cmg->h);
}

/*============================================================================
 * CMG Array and Steering Laws (Singularity Avoidance)
 *===========================================================================*/

#define LG_MAX_CMGS 6

typedef struct {
    lg_cmg_t cmgs[LG_MAX_CMGS];
    int n_cmgs;
    
    /* Jacobian: tau = C * delta (where delta is gimbal rate vector) */
    /* C[i][j] = d(tau_i)/d(delta_j) = (gimbal_axis_j x h_j)_i */
    float Jacobian[3][LG_MAX_CMGS];
    
    /* Singularity measure */
    float singularity_index; /* det(C*C^T), small = near singular */
} lg_cmg_array_t;

/* Standard 4-CMG pyramid (similar to ISS) */
static inline lg_cmg_array_t lg_cmg_array_iss_like(void) {
    lg_cmg_array_t arr;
    arr.n_cmgs = 4;
    
    float beta = 54.75f * LG_PI / 180.0f; /* Skew angle */
    
    for (int i = 0; i < 4; i++) {
        float angle = i * LG_PI / 2.0f;
        lg_vec3_t spin = {
            cosf(beta) * sinf(angle),
            cosf(beta) * cosf(angle),
            sinf(beta) * (i % 2 == 0 ? 1.0f : -1.0f) /* Alternate up/down */
        };
        lg_vec3_t gimbal = {
            -sinf(angle),
            cosf(angle),
            0.0f
        };
        
        arr.cmgs[i] = lg_cmg_high_power(spin, gimbal);
    }
    
    return arr;
}

/* Compute Jacobian matrix C (3 x n) */
static inline void lg_cmg_compute_jacobian(lg_cmg_array_t* arr) {
    for (int j = 0; j < arr->n_cmgs; j++) {
        lg_vec3_t h_vec = lg_vec3_scale(arr->cmgs[j].spin_axis, arr->cmgs[j].h);
        lg_vec3_t col = lg_vec3_cross(arr->cmgs[j].gimbal_axis, h_vec);
        
        arr->Jacobian[0][j] = col.x;
        arr->Jacobian[1][j] = col.y;
        arr->Jacobian[2][j] = col.z;
    }
    
    /* Compute singularity index: det(C*C^T)
     * Zero determinant → CMG array is in a singularity (cannot produce
     * torque in at least one direction).
     * C*C^T is symmetric; A is stored column-major: m[col*3 + row]. */
    lg_mat3_t A = {{0}};
    for (int j = 0; j < arr->n_cmgs; j++) {
        A.m[0] += arr->Jacobian[0][j] * arr->Jacobian[0][j];
        A.m[1] += arr->Jacobian[1][j] * arr->Jacobian[0][j];
        A.m[2] += arr->Jacobian[2][j] * arr->Jacobian[0][j];
        A.m[4] += arr->Jacobian[1][j] * arr->Jacobian[1][j];
        A.m[5] += arr->Jacobian[2][j] * arr->Jacobian[1][j];
        A.m[8] += arr->Jacobian[2][j] * arr->Jacobian[2][j];
    }
    A.m[3] = A.m[1];
    A.m[6] = A.m[2];
    A.m[7] = A.m[5];
    arr->singularity_index = lg_mat3_determinant(A);
}

/* Pseudo-inverse steering (Moore-Penrose) with singularity robustness */
static inline void lg_cmg_steering_pinv(lg_cmg_array_t* arr,
                                       const lg_vec3_t* tau_cmd,
                                       float dt) {
    (void)dt;
    lg_cmg_compute_jacobian(arr);

    /* delta = C^T * (C*C^T)^{-1} * tau_cmd */
    /* Add damping near singularity: (C*C^T + lambda*I)^{-1} */

    float CCT[3][3] = {0};
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            for (int k = 0; k < arr->n_cmgs; k++) {
                CCT[i][j] += arr->Jacobian[i][k] * arr->Jacobian[j][k];
            }
        }
    }

    /* Convert to lg_mat3_t (column-major) */
    lg_mat3_t m;
    m.m[0] = CCT[0][0]; m.m[1] = CCT[1][0]; m.m[2] = CCT[2][0];
    m.m[3] = CCT[0][1]; m.m[4] = CCT[1][1]; m.m[5] = CCT[2][1];
    m.m[6] = CCT[0][2]; m.m[7] = CCT[1][2]; m.m[8] = CCT[2][2];

    arr->singularity_index = lg_mat3_determinant(m);

    /* Add damping if singular */
    float lambda = 0.0f;
    if (arr->singularity_index < 0.01f) {
        lambda = 0.1f; /* Damping factor */
    }
    for (int i = 0; i < 3; i++) {
        m.m[i*3 + i] += lambda;
    }

    lg_mat3_t inv = lg_mat3_inverse(m);
    lg_vec3_t temp = lg_mat3_mul_vec3(inv, *tau_cmd);

    /* Compute gimbal rates: delta = C^T * temp */
    for (int j = 0; j < arr->n_cmgs; j++) {
        float rate = arr->Jacobian[0][j] * temp.x
                   + arr->Jacobian[1][j] * temp.y
                   + arr->Jacobian[2][j] * temp.z;
        arr->cmgs[j].gimbal_rate_cmd = rate;
    }
}

/* Null motion steering (singularity escape without producing torque) */
static inline void lg_cmg_null_motion(lg_cmg_array_t* arr, float strength) {
    /* Move along null space of Jacobian toward preferred gimbal angles */
    /* Only works for n > 3 (redundant CMGs) */
    if (arr->n_cmgs <= 3) return;

    lg_cmg_compute_jacobian(arr);

    /* Compute C*C^T and invert */
    float CCT[3][3] = {0};
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            for (int k = 0; k < arr->n_cmgs; k++) {
                CCT[i][j] += arr->Jacobian[i][k] * arr->Jacobian[j][k];
            }
        }
    }

    lg_mat3_t m;
    m.m[0] = CCT[0][0]; m.m[1] = CCT[1][0]; m.m[2] = CCT[2][0];
    m.m[3] = CCT[0][1]; m.m[4] = CCT[1][1]; m.m[5] = CCT[2][1];
    m.m[6] = CCT[0][2]; m.m[7] = CCT[1][2]; m.m[8] = CCT[2][2];

    lg_mat3_t inv = lg_mat3_inverse(m);

    /* Preferred gimbal angle is zero (home position) */
    float delta_bias[LG_MAX_CMGS];
    for (int j = 0; j < arr->n_cmgs; j++) {
        delta_bias[j] = -strength * arr->cmgs[j].gimbal_angle;
    }

    /* y = C * delta_bias */
    lg_vec3_t y = lg_vec3_zero();
    for (int j = 0; j < arr->n_cmgs; j++) {
        y.x += arr->Jacobian[0][j] * delta_bias[j];
        y.y += arr->Jacobian[1][j] * delta_bias[j];
        y.z += arr->Jacobian[2][j] * delta_bias[j];
    }

    /* z = inv(CCT) * y */
    lg_vec3_t z = lg_mat3_mul_vec3(inv, y);

    /* w = C^T * z */
    float w[LG_MAX_CMGS] = {0};
    for (int j = 0; j < arr->n_cmgs; j++) {
        w[j] = arr->Jacobian[0][j] * z.x
             + arr->Jacobian[1][j] * z.y
             + arr->Jacobian[2][j] * z.z;
    }

    /* delta_null = delta_bias - w */
    for (int j = 0; j < arr->n_cmgs; j++) {
        arr->cmgs[j].gimbal_rate_cmd += delta_bias[j] - w[j];
    }
}

/* Update all CMGs */
static inline void lg_cmg_array_update(lg_cmg_array_t* arr, float dt) {
    for (int i = 0; i < arr->n_cmgs; i++) {
        lg_cmg_update(&arr->cmgs[i], dt);
    }
}

/* Total stored momentum */
static inline lg_vec3_t lg_cmg_array_momentum(const lg_cmg_array_t* arr) {
    lg_vec3_t h = lg_vec3_zero();
    for (int i = 0; i < arr->n_cmgs; i++) {
        h = lg_vec3_add(h, lg_cmg_angular_momentum(&arr->cmgs[i]));
    }
    return h;
}

/* Total output torque */
static inline lg_vec3_t lg_cmg_array_torque(const lg_cmg_array_t* arr) {
    lg_vec3_t tau = lg_vec3_zero();
    for (int i = 0; i < arr->n_cmgs; i++) {
        tau = lg_vec3_add(tau, lg_cmg_output_torque(&arr->cmgs[i]));
    }
    return tau;
}

/*============================================================================
 * Gyroscopic Precession (Free Rotation of Rigid Body with Angular Momentum)
 * Demonstration: torque-free motion of axisymmetric body
 *===========================================================================*/

typedef struct {
    float I_axial;            /* Moment of inertia about symmetry axis */
    float I_transverse;       /* Moment about perpendicular axes */
    lg_vec3_t omega_body;   /* Angular velocity in body frame */
    lg_vec3_t H_body;       /* Constant angular momentum in inertial frame (mapped to body) */
    lg_quat_t attitude;       /* Current orientation */
} lg_precession_state_t;

/* Initialize with spin about axis, some nutation */
static inline lg_precession_state_t lg_precession_init(float I_axial, float I_trans, 
                                                      float spin_rate, float nutation_angle) {
    lg_precession_state_t p;
    p.I_axial = I_axial;
    p.I_transverse = I_trans;
    p.attitude = lg_quat_identity();
    
    /* Body omega: mostly spin_z plus small transverse component causing wobble */
    p.omega_body = lg_vec3(
        spin_rate * sinf(nutation_angle), /* Small x component */
        0.0f,
        spin_rate * cosf(nutation_angle)  /* Main spin */
    );
    
    /* H is constant in inertial frame, compute initial */
    lg_vec3_t H = lg_vec3_scale(p.omega_body, I_trans); /* Rough approx */
    H.z = p.omega_body.z * I_axial;
    p.H_body = H; /* Actually this rotates with body... store inertial separately */
    
    return p;
}

/* Torque-free Euler equations for axisymmetric body */
/* omega_dot_x = (I_trans - I_axial)/I_trans * omega_y * omega_z */
/* omega_dot_y = (I_axial - I_trans)/I_trans * omega_x * omega_z */
/* omega_dot_z = 0 */
static inline void lg_precession_step(lg_precession_state_t* p, float dt) {
    (void)p; /* float k = (p->I_transverse - p->I_axial) / p->I_transverse; */
    
    /* Exact solution for constant precession rate */
    float omega_p = p->omega_body.z * (p->I_axial - p->I_transverse) / p->I_transverse;
    
    /* Rotate omega_xy in body frame at rate omega_p */
    float wx = p->omega_body.x;
    float wy = p->omega_body.y;
    p->omega_body.x = wx * cosf(omega_p * dt) - wy * sinf(omega_p * dt);
    p->omega_body.y = wx * sinf(omega_p * dt) + wy * cosf(omega_p * dt);
    
    /* Update attitude quaternion */
    lg_vec3_t omega_world = lg_quat_rotate(p->attitude, p->omega_body);
    lg_quat_t delta = lg_quat_from_axis_angle(omega_world, dt);
    p->attitude = lg_quat_norm(lg_quat_mul(delta, p->attitude));
}

/*============================================================================
 * Integration with Rigid Body Dynamics
 * Apply gyroscopic coupling in integrator
 *===========================================================================*/

/* Modified semi-implicit Euler with RW/CMG angular momentum */
static inline void lg_integrate_gyroscopic(lg_body_t* body,
                                          lg_rw_array_t* rws,
                                          lg_cmg_array_t* cmgs,
                                          lg_transform_t* xform,
                                          float dt) {
    /* Total internal angular momentum */
    lg_vec3_t h_rw = rws ? lg_rw_array_momentum(rws) : lg_vec3_zero();
    lg_vec3_t h_cmg = cmgs ? lg_cmg_array_momentum(cmgs) : lg_vec3_zero();
    lg_vec3_t h_int = lg_vec3_add(h_rw, h_cmg);
    
    /* Euler equation: I*alpha + omega x (I*omega + h_int) = tau_ext - tau_reaction */
    lg_vec3_t I_omega = lg_vec3_mul(body->angular_velocity, body->inertia);
    lg_vec3_t h_total = lg_vec3_add(I_omega, h_int);
    lg_vec3_t gyro = lg_vec3_cross(body->angular_velocity, h_total);
    
    /* External torques plus reaction motor torques */
    lg_vec3_t tau_motors = rws ? lg_vec3_neg(lg_rw_body_torque(&rws->wheels[0])) : lg_vec3_zero(); /* Simplified */
    lg_vec3_t tau_cmg = cmgs ? lg_cmg_array_torque(cmgs) : lg_vec3_zero();
    
    lg_vec3_t net_tau = lg_vec3_sub(lg_vec3_add(body->torque, tau_cmg), gyro);
    net_tau = lg_vec3_add(net_tau, tau_motors);
    
    /* Angular acceleration */
    lg_vec3_t alpha = lg_vec3_mul(net_tau, body->inv_inertia);
    
    /* Semi-implicit update */
    body->angular_velocity = lg_vec3_add(body->angular_velocity, lg_vec3_scale(alpha, dt));
    
    /* Quaternion update */
    float angle = lg_vec3_len(body->angular_velocity) * dt;
    if (angle > 1e-6f) {
        lg_vec3_t axis = lg_vec3_norm(body->angular_velocity);
        lg_quat_t delta = lg_quat_from_axis_angle(axis, angle);
        xform->rotation = lg_quat_norm(lg_quat_mul(delta, xform->rotation));
    }
    
    /* Update wheel speeds */
    if (rws) lg_rw_array_update(rws, dt);
    if (cmgs) lg_cmg_array_update(cmgs, dt);
    
    /* Clear external torques */
    lg_body_clear_forces(body);
}

/*============================================================================
 * Gravity Gradient Torque
 * Torque on a spacecraft with principal inertia I due to gravity gradient.
 * tau = (3*mu / |r|^5) * r x (I * r)
 *===========================================================================*/
static inline lg_vec3_t lg_gravity_gradient_torque(lg_vec3_t r,
                                                    lg_vec3_t inertia,
                                                    float mu) {
    float r2 = lg_vec3_len_sq(r);
    float r5 = r2 * r2 * sqrtf(r2);
    if (r5 < 1e-12f) return lg_vec3_zero();
    
    /* I * r for diagonal inertia tensor */
    lg_vec3_t Ir = {
        inertia.x * r.x,
        inertia.y * r.y,
        inertia.z * r.z
    };
    
    float factor = 3.0f * mu / r5;
    return lg_vec3_scale(lg_vec3_cross(r, Ir), factor);
}

/*============================================================================
 * Momentum Management (keeping total H within envelope)
 *===========================================================================*/
static inline void lg_momentum_manage(lg_rw_array_t* rws,
                                      lg_cmg_array_t* cmgs,
                                      const lg_vec3_t* gravity_gradient_torque,
                                      float dt) {
    (void)dt;
    /* Dump momentum using environmental torques when approaching saturation */
    lg_vec3_t h_total = rws ? lg_rw_array_momentum(rws) : lg_vec3_zero();
    if (cmgs) h_total = lg_vec3_add(h_total, lg_cmg_array_momentum(cmgs));

    float h_mag = lg_vec3_len(h_total);
    if (h_mag < 1e-6f) return;

    /* Compute total momentum capacity */
    float capacity = 0.0f;
    if (rws) capacity += rws->total_momentum_capacity;
    if (cmgs) {
        for (int i = 0; i < cmgs->n_cmgs; i++) {
            capacity += cmgs->cmgs[i].h;
        }
    }

    float threshold = 0.5f * capacity;
    if (h_mag <= threshold) return;

    /* Desaturation direction: opposite to total momentum */
    lg_vec3_t desat_dir = lg_vec3_scale(h_total, -1.0f / h_mag);

    /* Scale by available external torque magnitude (up to 10% of it) */
    float ext_mag = gravity_gradient_torque ? lg_vec3_len(*gravity_gradient_torque) : 0.0f;
    float dump_strength = fminf(ext_mag * 0.1f, (h_mag - threshold) * 0.1f);
    if (dump_strength < 1e-6f) return;

    lg_vec3_t tau_dump = lg_vec3_scale(desat_dir, dump_strength);

    /* Command reaction wheels to dump momentum */
    if (rws) {
        lg_rw_array_command_torque(rws, &tau_dump);
    }
}

#ifdef __cplusplus
}
#endif

#endif /* LAGRANGE_ATTITUDE_CONTROL_H */

