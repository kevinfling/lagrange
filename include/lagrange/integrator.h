/*
 * Lagrange Physics Library - Integration Methods
 * Numerical integration for physics simulation
 */

#ifndef LAGRANGE_INTEGRATOR_H
#define LAGRANGE_INTEGRATOR_H

#include "types.h"
#include "math.h"
#include "body.h"
#include "transform.h"

#ifdef __cplusplus
extern "C" {
#endif

/*============================================================================
 * Integration Types
 *===========================================================================*/

typedef enum {
    LG_INTEGRATOR_EXPLICIT_EULER,
    LG_INTEGRATOR_SEMI_IMPLICIT_EULER,
    LG_INTEGRATOR_VERLET,
    LG_INTEGRATOR_VELOCITY_VERLET,
    LG_INTEGRATOR_RK4
} lg_integrator_type_t;

/*============================================================================
 * Explicit Euler (First Order)
 * v' = v + a * dt
 * x' = x + v * dt
 *===========================================================================*/
static inline void lg_integrate_explicit_euler(
    lg_body_t* body,
    lg_transform_t* transform,
    float dt
) {
    if (body->type != LG_BODY_DYNAMIC || body->is_sleeping) return;
    
    /* Update velocity */
    lg_vec3_t accel = lg_vec3_scale(body->force, body->inv_mass);
    body->velocity = lg_vec3_add(body->velocity, lg_vec3_scale(accel, dt));
    
    /* Apply damping */
    body->velocity = lg_vec3_scale(body->velocity, 1.0f - body->linear_damping * dt);
    
    /* Update position */
    transform->position = lg_vec3_add(transform->position, lg_vec3_scale(body->velocity, dt));
    
    /* Angular velocity */
    lg_vec3_t ang_accel = lg_vec3_mul(body->torque, body->inv_inertia);
    body->angular_velocity = lg_vec3_add(body->angular_velocity, lg_vec3_scale(ang_accel, dt));
    body->angular_velocity = lg_vec3_scale(body->angular_velocity, 1.0f - body->angular_damping * dt);
    
    /* Update rotation (simplified) */
    float angle = lg_vec3_len(body->angular_velocity) * dt;
    if (angle > 1e-6f) {
        lg_vec3_t axis = lg_vec3_norm(body->angular_velocity);
        lg_quat_t delta = lg_quat_from_axis_angle(axis, angle);
        transform->rotation = lg_quat_mul(delta, transform->rotation);
    }
}

/*============================================================================
 * Semi-Implicit Euler (Symplectic Euler)
 * v' = v + a * dt
 * x' = x + v' * dt
 * More stable than explicit Euler, preserves energy better
 *===========================================================================*/
static inline void lg_integrate_semi_implicit_euler(
    lg_body_t* body,
    lg_transform_t* transform,
    float dt
) {
    if (body->type != LG_BODY_DYNAMIC || body->is_sleeping) return;
    
    /* Store previous position for velocity verlet if needed */
    body->prev_position = transform->position;
    
    /* Update velocity first */
    lg_vec3_t accel = lg_vec3_scale(body->force, body->inv_mass);
    body->velocity = lg_vec3_add(body->velocity, lg_vec3_scale(accel, dt));
    body->velocity = lg_vec3_scale(body->velocity, 1.0f - body->linear_damping * dt);
    
    /* Then update position with new velocity */
    transform->position = lg_vec3_add(transform->position, lg_vec3_scale(body->velocity, dt));
    
    /* Angular part */
    lg_vec3_t ang_accel = lg_vec3_mul(body->torque, body->inv_inertia);
    body->angular_velocity = lg_vec3_add(body->angular_velocity, lg_vec3_scale(ang_accel, dt));
    body->angular_velocity = lg_vec3_scale(body->angular_velocity, 1.0f - body->angular_damping * dt);
    
    /* Update rotation */
    float angle = lg_vec3_len(body->angular_velocity) * dt;
    if (angle > 1e-6f) {
        lg_vec3_t axis = lg_vec3_norm(body->angular_velocity);
        lg_quat_t delta = lg_quat_from_axis_angle(axis, angle);
        transform->rotation = lg_quat_norm(lg_quat_mul(delta, transform->rotation));
    }
}

/*============================================================================
 * Velocity Verlet
 * x' = x + v * dt + 0.5 * a * dt²
 * v' = v + 0.5 * (a + a') * dt
 * Good energy conservation, time-reversible
 *===========================================================================*/
typedef struct {
    lg_vec3_t prev_accel;
    bool first_step;
} lg_verlet_state_t;

static inline void lg_integrate_velocity_verlet(
    lg_body_t* body,
    lg_transform_t* transform,
    lg_verlet_state_t* state,
    float dt
) {
    if (body->type != LG_BODY_DYNAMIC || body->is_sleeping) return;
    
    float dt_sq = dt * dt;
    
    /* Current acceleration */
    lg_vec3_t accel = lg_vec3_scale(body->force, body->inv_mass);
    
    if (state->first_step) {
        /* First step: use semi-implicit Euler */
        state->prev_accel = accel;
        state->first_step = false;
        lg_integrate_semi_implicit_euler(body, transform, dt);
        return;
    }
    
    /* Position update */
    lg_vec3_t pos_change = lg_vec3_add(
        lg_vec3_scale(body->velocity, dt),
        lg_vec3_scale(state->prev_accel, 0.5f * dt_sq)
    );
    transform->position = lg_vec3_add(transform->position, pos_change);
    
    /* Velocity update (half step with old accel) */
    body->velocity = lg_vec3_add(body->velocity, lg_vec3_scale(state->prev_accel, 0.5f * dt));
    
    /* Apply damping */
    body->velocity = lg_vec3_scale(body->velocity, 1.0f - body->linear_damping * dt);
    
    /* Velocity update (half step with new accel) */
    body->velocity = lg_vec3_add(body->velocity, lg_vec3_scale(accel, 0.5f * dt));
    
    /* Store acceleration for next step */
    state->prev_accel = accel;
    
    /* Angular part (simplified, using semi-implicit) */
    lg_vec3_t ang_accel = lg_vec3_mul(body->torque, body->inv_inertia);
    body->angular_velocity = lg_vec3_add(body->angular_velocity, lg_vec3_scale(ang_accel, dt));
    body->angular_velocity = lg_vec3_scale(body->angular_velocity, 1.0f - body->angular_damping * dt);
    
    float angle = lg_vec3_len(body->angular_velocity) * dt;
    if (angle > 1e-6f) {
        lg_vec3_t axis = lg_vec3_norm(body->angular_velocity);
        lg_quat_t delta = lg_quat_from_axis_angle(axis, angle);
        transform->rotation = lg_quat_norm(lg_quat_mul(delta, transform->rotation));
    }
}

/*============================================================================
 * RK4 (Runge-Kutta 4th Order)
 * Higher accuracy, but more expensive
 *===========================================================================*/

typedef struct {
    lg_vec3_t position;
    lg_vec3_t velocity;
} lg_rk4_state_t;

typedef struct {
    lg_vec3_t accel;
} lg_rk4_derivative_t;

static inline lg_rk4_derivative_t lg_rk4_eval(
    const lg_rk4_state_t* state,
    const lg_body_t* body,
    float dt,
    const lg_rk4_derivative_t* derivative
) {
    lg_rk4_state_t s = *state;
    
    if (derivative) {
        s.velocity = lg_vec3_add(s.velocity, lg_vec3_scale(derivative->accel, dt));
    }
    
    lg_rk4_derivative_t result;
    result.accel = lg_vec3_scale(body->force, body->inv_mass);
    
    return result;
}

static inline void lg_integrate_rk4(
    lg_body_t* body,
    lg_transform_t* transform,
    float dt
) {
    if (body->type != LG_BODY_DYNAMIC || body->is_sleeping) return;
    
    lg_rk4_state_t state = {
        .position = transform->position,
        .velocity = body->velocity
    };
    
    /* Four evaluations */
    lg_rk4_derivative_t a = lg_rk4_eval(&state, body, 0.0f, NULL);
    lg_rk4_derivative_t b = lg_rk4_eval(&state, body, dt * 0.5f, &a);
    lg_rk4_derivative_t c = lg_rk4_eval(&state, body, dt * 0.5f, &b);
    lg_rk4_derivative_t d = lg_rk4_eval(&state, body, dt, &c);
    
    /* Combine */
    lg_vec3_t dv = lg_vec3_scale(
        lg_vec3_add(
            lg_vec3_add(a.accel, d.accel),
            lg_vec3_scale(lg_vec3_add(b.accel, c.accel), 2.0f)
        ),
        1.0f / 6.0f * dt
    );
    
    /* Apply velocity change */
    body->velocity = lg_vec3_add(body->velocity, dv);
    body->velocity = lg_vec3_scale(body->velocity, 1.0f - body->linear_damping * dt);
    
    /* Position update using new velocity */
    transform->position = lg_vec3_add(transform->position, lg_vec3_scale(body->velocity, dt));
    
    /* Angular part (simplified, using semi-implicit) */
    lg_vec3_t ang_accel = lg_vec3_mul(body->torque, body->inv_inertia);
    body->angular_velocity = lg_vec3_add(body->angular_velocity, lg_vec3_scale(ang_accel, dt));
    body->angular_velocity = lg_vec3_scale(body->angular_velocity, 1.0f - body->angular_damping * dt);
    
    float angle = lg_vec3_len(body->angular_velocity) * dt;
    if (angle > 1e-6f) {
        lg_vec3_t axis = lg_vec3_norm(body->angular_velocity);
        lg_quat_t delta = lg_quat_from_axis_angle(axis, angle);
        transform->rotation = lg_quat_norm(lg_quat_mul(delta, transform->rotation));
    }
}

/*============================================================================
 * Unified Integration Interface
 *===========================================================================*/

static inline void lg_integrate(
    lg_integrator_type_t type,
    lg_body_t* body,
    lg_transform_t* transform,
    void* state,  /* Type-specific state (e.g., lg_verlet_state_t) */
    float dt
) {
    switch (type) {
        case LG_INTEGRATOR_EXPLICIT_EULER:
            lg_integrate_explicit_euler(body, transform, dt);
            break;
        case LG_INTEGRATOR_SEMI_IMPLICIT_EULER:
            lg_integrate_semi_implicit_euler(body, transform, dt);
            break;
        case LG_INTEGRATOR_VELOCITY_VERLET:
            lg_integrate_velocity_verlet(body, transform, (lg_verlet_state_t*)state, dt);
            break;
        case LG_INTEGRATOR_RK4:
            lg_integrate_rk4(body, transform, dt);
            break;
        default:
            lg_integrate_semi_implicit_euler(body, transform, dt);
            break;
    }
}

/*============================================================================
 * Utility Functions
 *===========================================================================*/

/* Check if body should sleep based on velocity magnitude */
static inline bool lg_integrator_check_sleep(
    lg_body_t* body,
    float threshold_sq,  /* Squared velocity threshold */
    float dt
) {
    if (!body->allow_sleep || body->type != LG_BODY_DYNAMIC) return false;
    
    float linear_sq = lg_vec3_len_sq(body->velocity);
    float angular_sq = lg_vec3_len_sq(body->angular_velocity);
    
    if (linear_sq < threshold_sq && angular_sq < threshold_sq) {
        body->sleep_timer += dt;
        return body->sleep_timer > 0.5f; /* Sleep after 0.5s of low activity */
    } else {
        body->sleep_timer = 0.0f;
        return false;
    }
}

#ifdef __cplusplus
}
#endif

#endif /* LAGRANGE_INTEGRATOR_H */
