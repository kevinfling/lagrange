/*
 * Lagrange Physics Library - Simulation
 * Physics step, collision detection, and constraint solving
 */

#ifndef LAGRANGE_SIM_H
#define LAGRANGE_SIM_H

#include "world.h"
#include "gravity.h"
#include <stdio.h>  /* DEBUG */

#ifdef __cplusplus
extern "C" {
#endif

/*============================================================================
 * Collision Detection
 *===========================================================================*/

typedef struct {
    lg_entity_t entity_a;
    lg_entity_t entity_b;
    lg_vec3_t point;        /* Contact point in world space */
    lg_vec3_t normal;       /* Normal pointing from A to B */
    float penetration;      /* Penetration depth (positive = overlap) */
    float restitution;      /* Combined restitution */
    float friction;         /* Combined friction */
} lg_contact_t;

/* Sphere-sphere collision */
static inline bool lg_collide_spheres(
    lg_vec3_t pos_a, float radius_a,
    lg_vec3_t pos_b, float radius_b,
    lg_contact_t* out_contact
) {
    lg_vec3_t delta = lg_vec3_sub(pos_b, pos_a);
    float dist_sq = lg_vec3_len_sq(delta);
    float radius_sum = radius_a + radius_b;
    
    if (dist_sq > radius_sum * radius_sum) {
        return false; /* No collision */
    }
    
    float dist = sqrtf(dist_sq);
    if (dist < 1e-6f) {
        /* Centers coincide, arbitrary normal */
        out_contact->normal = lg_vec3(0.0f, 1.0f, 0.0f);
        out_contact->penetration = radius_sum;
    } else {
        out_contact->normal = lg_vec3_scale(delta, 1.0f / dist);
        out_contact->penetration = radius_sum - dist;
    }
    
    /* Contact point is on the surface of A */
    out_contact->point = lg_vec3_add(pos_a, lg_vec3_scale(out_contact->normal, radius_a));
    
    return true;
}

/* Sphere-plane collision */
static inline bool lg_collide_sphere_plane(
    lg_vec3_t sphere_pos, float sphere_radius,
    lg_vec3_t plane_normal, float plane_distance,
    lg_contact_t* out_contact
) {
    plane_normal = lg_vec3_norm(plane_normal);
    
    /* Distance from sphere center to plane */
    float dist = lg_vec3_dot(sphere_pos, plane_normal) - plane_distance;
    
    if (dist > sphere_radius) {
        return false; /* Too far from plane */
    }
    
    out_contact->normal = plane_normal;
    out_contact->penetration = sphere_radius - dist;
    out_contact->point = lg_vec3_sub(sphere_pos, lg_vec3_scale(plane_normal, dist));
    
    return true;
}

/* Box-plane collision (simplified - just check vertices) */
static inline bool lg_collide_box_plane(
    lg_vec3_t box_pos, lg_vec3_t box_half_extents,
    lg_vec3_t plane_normal, float plane_distance,
    lg_contact_t* out_contact
) {
    plane_normal = lg_vec3_norm(plane_normal);
    
    /* Find deepest point on box */
    float min_dist = 1e30f;
    lg_vec3_t deepest_point;
    
    /* Check all 8 corners */
    for (int i = 0; i < 8; i++) {
        lg_vec3_t corner = box_pos;
        corner.x += ((i & 1) ? 1.0f : -1.0f) * box_half_extents.x;
        corner.y += ((i & 2) ? 1.0f : -1.0f) * box_half_extents.y;
        corner.z += ((i & 4) ? 1.0f : -1.0f) * box_half_extents.z;
        
        float dist = lg_vec3_dot(corner, plane_normal) - plane_distance;
        if (dist < min_dist) {
            min_dist = dist;
            deepest_point = corner;
        }
    }
    
    if (min_dist > 0.0f) {
        return false; /* No collision */
    }
    
    out_contact->normal = plane_normal;
    out_contact->penetration = -min_dist;
    out_contact->point = deepest_point;
    
    return true;
}

/* Sphere-box collision (AABB-based) */
static inline bool lg_collide_sphere_box(
    lg_vec3_t sphere_pos, float sphere_radius,
    lg_vec3_t box_pos, lg_vec3_t box_half_extents,
    lg_contact_t* out_contact
) {
    lg_vec3_t closest;
    closest.x = fmaxf(box_pos.x - box_half_extents.x, fminf(sphere_pos.x, box_pos.x + box_half_extents.x));
    closest.y = fmaxf(box_pos.y - box_half_extents.y, fminf(sphere_pos.y, box_pos.y + box_half_extents.y));
    closest.z = fmaxf(box_pos.z - box_half_extents.z, fminf(sphere_pos.z, box_pos.z + box_half_extents.z));

    lg_vec3_t delta = lg_vec3_sub(sphere_pos, closest);
    float dist_sq = lg_vec3_len_sq(delta);
    float radius_sq = sphere_radius * sphere_radius;

    if (dist_sq > radius_sq) {
        return false;
    }

    float dist = sqrtf(dist_sq);
    if (dist < 1e-6f) {
        out_contact->normal = lg_vec3(1.0f, 0.0f, 0.0f);
        out_contact->penetration = sphere_radius;
    } else {
        out_contact->normal = lg_vec3_scale(delta, 1.0f / dist);
        out_contact->penetration = sphere_radius - dist;
    }
    out_contact->point = closest;
    return true;
}

/* Box-box collision (AABB overlap) */
static inline bool lg_collide_box_box(
    lg_vec3_t pos_a, lg_vec3_t half_a,
    lg_vec3_t pos_b, lg_vec3_t half_b,
    lg_contact_t* out_contact
) {
    float dx = half_a.x + half_b.x - fabsf(pos_a.x - pos_b.x);
    float dy = half_a.y + half_b.y - fabsf(pos_a.y - pos_b.y);
    float dz = half_a.z + half_b.z - fabsf(pos_a.z - pos_b.z);

    if (dx <= 0.0f || dy <= 0.0f || dz <= 0.0f) {
        return false;
    }

    if (dx <= dy && dx <= dz) {
        out_contact->penetration = dx;
        float sign = (pos_a.x > pos_b.x) ? 1.0f : -1.0f;
        out_contact->normal = lg_vec3(sign, 0.0f, 0.0f);
        float x = (pos_a.x - sign * half_a.x + pos_b.x + sign * half_b.x) * 0.5f;
        float min_y = fmaxf(pos_a.y - half_a.y, pos_b.y - half_b.y);
        float max_y = fminf(pos_a.y + half_a.y, pos_b.y + half_b.y);
        float min_z = fmaxf(pos_a.z - half_a.z, pos_b.z - half_b.z);
        float max_z = fminf(pos_a.z + half_a.z, pos_b.z + half_b.z);
        out_contact->point = lg_vec3(x, (min_y + max_y) * 0.5f, (min_z + max_z) * 0.5f);
    } else if (dy <= dx && dy <= dz) {
        out_contact->penetration = dy;
        float sign = (pos_a.y > pos_b.y) ? 1.0f : -1.0f;
        out_contact->normal = lg_vec3(0.0f, sign, 0.0f);
        float y = (pos_a.y - sign * half_a.y + pos_b.y + sign * half_b.y) * 0.5f;
        float min_x = fmaxf(pos_a.x - half_a.x, pos_b.x - half_b.x);
        float max_x = fminf(pos_a.x + half_a.x, pos_b.x + half_b.x);
        float min_z = fmaxf(pos_a.z - half_a.z, pos_b.z - half_b.z);
        float max_z = fminf(pos_a.z + half_a.z, pos_b.z + half_b.z);
        out_contact->point = lg_vec3((min_x + max_x) * 0.5f, y, (min_z + max_z) * 0.5f);
    } else {
        out_contact->penetration = dz;
        float sign = (pos_a.z > pos_b.z) ? 1.0f : -1.0f;
        out_contact->normal = lg_vec3(0.0f, 0.0f, sign);
        float z = (pos_a.z - sign * half_a.z + pos_b.z + sign * half_b.z) * 0.5f;
        float min_x = fmaxf(pos_a.x - half_a.x, pos_b.x - half_b.x);
        float max_x = fminf(pos_a.x + half_a.x, pos_b.x + half_b.x);
        float min_y = fmaxf(pos_a.y - half_a.y, pos_b.y - half_b.y);
        float max_y = fminf(pos_a.y + half_a.y, pos_b.y + half_b.y);
        out_contact->point = lg_vec3((min_x + max_x) * 0.5f, (min_y + max_y) * 0.5f, z);
    }
    return true;
}

/* Sphere-capsule collision */
static inline bool lg_collide_sphere_capsule(
    lg_vec3_t sphere_pos, float sphere_radius,
    lg_vec3_t cap_pos, lg_vec3_t cap_axis, float cap_half_height, float cap_radius,
    lg_contact_t* out_contact
) {
    lg_vec3_t a = lg_vec3_sub(cap_pos, lg_vec3_scale(cap_axis, cap_half_height));
    lg_vec3_t b = lg_vec3_add(cap_pos, lg_vec3_scale(cap_axis, cap_half_height));
    lg_vec3_t ab = lg_vec3_sub(b, a);
    lg_vec3_t ap = lg_vec3_sub(sphere_pos, a);
    float ab_len_sq = lg_vec3_len_sq(ab);
    float t = 0.0f;
    if (ab_len_sq > 1e-6f) {
        t = fmaxf(0.0f, fminf(1.0f, lg_vec3_dot(ap, ab) / ab_len_sq));
    }
    lg_vec3_t closest = lg_vec3_add(a, lg_vec3_scale(ab, t));
    lg_vec3_t delta = lg_vec3_sub(sphere_pos, closest);
    float dist_sq = lg_vec3_len_sq(delta);
    float radius_sum = sphere_radius + cap_radius;

    if (dist_sq > radius_sum * radius_sum) {
        return false;
    }

    float dist = sqrtf(dist_sq);
    if (dist < 1e-6f) {
        out_contact->normal = cap_axis;
        out_contact->penetration = radius_sum;
    } else {
        out_contact->normal = lg_vec3_scale(delta, 1.0f / dist);
        out_contact->penetration = radius_sum - dist;
    }
    out_contact->point = closest;
    return true;
}

/* Capsule-capsule collision */
static inline bool lg_collide_capsule_capsule(
    lg_vec3_t pos_a, lg_vec3_t axis_a, float half_h_a, float radius_a,
    lg_vec3_t pos_b, lg_vec3_t axis_b, float half_h_b, float radius_b,
    lg_contact_t* out_contact
) {
    lg_vec3_t p1 = lg_vec3_sub(pos_a, lg_vec3_scale(axis_a, half_h_a));
    lg_vec3_t p2 = lg_vec3_add(pos_a, lg_vec3_scale(axis_a, half_h_a));
    lg_vec3_t q1 = lg_vec3_sub(pos_b, lg_vec3_scale(axis_b, half_h_b));
    lg_vec3_t q2 = lg_vec3_add(pos_b, lg_vec3_scale(axis_b, half_h_b));

    lg_vec3_t u = lg_vec3_sub(p2, p1);
    lg_vec3_t v = lg_vec3_sub(q2, q1);
    lg_vec3_t w = lg_vec3_sub(p1, q1);

    float a = lg_vec3_dot(u, u);
    float b = lg_vec3_dot(u, v);
    float c = lg_vec3_dot(v, v);
    float d = lg_vec3_dot(u, w);
    float e = lg_vec3_dot(v, w);
    float D = a * c - b * b;
    float sD = D, tD = D;
    float sN = 0.0f, tN = 0.0f;
    const float SMALL_NUM = 1e-6f;

    if (D < SMALL_NUM) {
        sN = 0.0f;
        sD = 1.0f;
        tN = e;
        tD = c;
    } else {
        sN = b * e - c * d;
        tN = a * e - b * d;
        if (sN < 0.0f) {
            sN = 0.0f;
            tN = e;
            tD = c;
        } else if (sN > sD) {
            sN = sD;
            tN = e + b;
            tD = c;
        }
    }

    if (tN < 0.0f) {
        tN = 0.0f;
        if (-d < 0.0f) {
            sN = 0.0f;
        } else if (-d > a) {
            sN = sD;
        } else {
            sN = -d;
            sD = a;
        }
    } else if (tN > tD) {
        tN = tD;
        if ((-d + b) < 0.0f) {
            sN = 0.0f;
        } else if ((-d + b) > a) {
            sN = sD;
        } else {
            sN = (-d + b);
            sD = a;
        }
    }

    float sc = (fabsf(sN) < SMALL_NUM) ? 0.0f : sN / sD;
    float tc = (fabsf(tN) < SMALL_NUM) ? 0.0f : tN / tD;

    lg_vec3_t c1 = lg_vec3_add(p1, lg_vec3_scale(u, sc));
    lg_vec3_t c2 = lg_vec3_add(q1, lg_vec3_scale(v, tc));
    lg_vec3_t delta = lg_vec3_sub(c1, c2);
    float dist_sq = lg_vec3_len_sq(delta);
    float radius_sum = radius_a + radius_b;

    if (dist_sq > radius_sum * radius_sum) {
        return false;
    }

    float dist = sqrtf(dist_sq);
    if (dist < 1e-6f) {
        out_contact->normal = lg_vec3_cross(axis_a, axis_b);
        if (lg_vec3_len_sq(out_contact->normal) < 1e-6f) {
            out_contact->normal = lg_vec3(1.0f, 0.0f, 0.0f);
        } else {
            out_contact->normal = lg_vec3_norm(out_contact->normal);
        }
        out_contact->penetration = radius_sum;
    } else {
        out_contact->normal = lg_vec3_scale(delta, 1.0f / dist);
        out_contact->penetration = radius_sum - dist;
    }
    out_contact->point = lg_vec3_add(c2, lg_vec3_scale(out_contact->normal, radius_b));
    return true;
}

/* Sphere-cylinder collision */
static inline bool lg_collide_sphere_cylinder(
    lg_vec3_t sphere_pos, float sphere_radius,
    lg_vec3_t cyl_pos, lg_vec3_t cyl_axis, float cyl_half_height, float cyl_radius,
    lg_contact_t* out_contact
) {
    cyl_axis = lg_vec3_norm(cyl_axis);
    
    lg_vec3_t to_sphere = lg_vec3_sub(sphere_pos, cyl_pos);
    float axial = lg_vec3_dot(to_sphere, cyl_axis);
    lg_vec3_t axial_vec = lg_vec3_scale(cyl_axis, axial);
    lg_vec3_t radial_vec = lg_vec3_sub(to_sphere, axial_vec);
    float radial_dist = lg_vec3_len(radial_vec);
    
    lg_vec3_t normal;
    float penetration;
    lg_vec3_t point;
    
    if (axial >= -cyl_half_height && axial <= cyl_half_height) {
        /* Side collision */
        if (radial_dist < 1e-6f) {
            normal = lg_vec3(1.0f, 0.0f, 0.0f);
        } else {
            normal = lg_vec3_scale(radial_vec, 1.0f / radial_dist);
        }
        penetration = (sphere_radius + cyl_radius) - radial_dist;
        if (penetration < 0.0f) return false;
        point = lg_vec3_add(cyl_pos, lg_vec3_add(axial_vec, lg_vec3_scale(normal, cyl_radius)));
    } else {
        float cap_sign = (axial > 0.0f) ? 1.0f : -1.0f;
        float axial_from_cap = fabsf(axial) - cyl_half_height;
        lg_vec3_t cap_center = lg_vec3_add(cyl_pos, lg_vec3_scale(cyl_axis, cap_sign * cyl_half_height));
        
        if (radial_dist <= cyl_radius) {
            /* Cap face collision */
            normal = lg_vec3_scale(cyl_axis, cap_sign);
            penetration = sphere_radius - axial_from_cap;
            if (penetration < 0.0f) return false;
            point = lg_vec3_add(cap_center, radial_vec);
        } else {
            /* Edge collision */
            float radial_from_edge = radial_dist - cyl_radius;
            float dist = sqrtf(radial_from_edge * radial_from_edge + axial_from_cap * axial_from_cap);
            if (dist > sphere_radius) return false;
            if (dist < 1e-6f) {
                normal = lg_vec3_scale(cyl_axis, cap_sign);
            } else {
                lg_vec3_t edge_point = lg_vec3_add(cap_center, lg_vec3_scale(radial_vec, cyl_radius / radial_dist));
                normal = lg_vec3_norm(lg_vec3_sub(sphere_pos, edge_point));
            }
            penetration = sphere_radius - dist;
            point = lg_vec3_add(sphere_pos, lg_vec3_scale(normal, -sphere_radius));
        }
    }
    
    out_contact->normal = normal;
    out_contact->penetration = penetration;
    out_contact->point = point;
    return true;
}

/*============================================================================
 * Collision Response
 *===========================================================================*/

static inline void lg_resolve_contact(lg_world_t* world, const lg_contact_t* contact) {
    lg_body_t* body_a = lg_get_body(world, contact->entity_a);
    lg_body_t* body_b = lg_get_body(world, contact->entity_b);
    lg_transform_t* trans_a = lg_get_transform(world, contact->entity_a);
    lg_transform_t* trans_b = lg_get_transform(world, contact->entity_b);
    

    
    if (!body_a && !body_b) return; /* Both static */
    
    /* Contact vectors from center of mass to contact point */
    lg_vec3_t r_a = lg_vec3_zero(), r_b = lg_vec3_zero();
    if (body_a && trans_a) {
        r_a = lg_vec3_sub(contact->point, trans_a->position);
    }
    if (body_b && trans_b) {
        r_b = lg_vec3_sub(contact->point, trans_b->position);
    }
    
    /* Relative velocity at contact point including angular motion */
    lg_vec3_t vel_a = body_a ? body_a->velocity : lg_vec3_zero();
    lg_vec3_t vel_b = body_b ? body_b->velocity : lg_vec3_zero();
    
    if (body_a) {
        vel_a = lg_vec3_add(vel_a, lg_vec3_cross(body_a->angular_velocity, r_a));
    }
    if (body_b) {
        vel_b = lg_vec3_add(vel_b, lg_vec3_cross(body_b->angular_velocity, r_b));
    }
    
    /* Relative velocity of A with respect to B */
    lg_vec3_t rel_vel = lg_vec3_sub(vel_a, vel_b);
    float vel_along_normal = lg_vec3_dot(rel_vel, contact->normal);
    
    /* Don't resolve if velocities are separating (A moving away from B along normal) */
    if (vel_along_normal > 0.0f) {
        return;
    }
    
    /* Mass and inertia properties */
    float inv_mass_a = body_a ? body_a->inv_mass : 0.0f;
    float inv_mass_b = body_b ? body_b->inv_mass : 0.0f;
    lg_vec3_t inv_inertia_a = body_a ? body_a->inv_inertia : lg_vec3_zero();
    lg_vec3_t inv_inertia_b = body_b ? body_b->inv_inertia : lg_vec3_zero();
    
    /* Compute effective mass along normal (including angular effects) */
    float eff_mass_normal = inv_mass_a + inv_mass_b;
    
    /* Add angular contribution: (r x n)^2 / I for each body */
    if (body_a) {
        lg_vec3_t r_cross_n = lg_vec3_cross(r_a, contact->normal);
        eff_mass_normal += r_cross_n.x * r_cross_n.x * inv_inertia_a.x
                        + r_cross_n.y * r_cross_n.y * inv_inertia_a.y
                        + r_cross_n.z * r_cross_n.z * inv_inertia_a.z;
    }
    if (body_b) {
        lg_vec3_t r_cross_n = lg_vec3_cross(r_b, contact->normal);
        eff_mass_normal += r_cross_n.x * r_cross_n.x * inv_inertia_b.x
                        + r_cross_n.y * r_cross_n.y * inv_inertia_b.y
                        + r_cross_n.z * r_cross_n.z * inv_inertia_b.z;
    }
    
    if (eff_mass_normal < 1e-10f) return; /* Avoid division by zero */
    
    /*=== NORMAL IMPULSE (bounce) ===*/
    float e = contact->restitution;
    /* 
     * Standard impulse formula: j = -(1+e) * v_rel / (1/m1 + 1/m2)
     * Since body A receives velocity -= impulse/mass, we need impulse 
     * in direction of normal to push A away from B when approaching.
     * vel_along_normal < 0 means approaching, so j should be positive
     * to push A along normal direction (away from B).
     */
    float j_normal = -(1.0f + e) * vel_along_normal / eff_mass_normal;
    /* Impulse should repel A from B along normal direction */
    if (j_normal < 0.0f) return; /* Don't attract bodies together */
    
    lg_vec3_t impulse_normal = lg_vec3_scale(contact->normal, j_normal);
    
    /*=== FRICTION IMPULSE (tangential) ===*/
    /* Compute tangent direction */
    lg_vec3_t vel_tangent = lg_vec3_sub(rel_vel, lg_vec3_scale(contact->normal, vel_along_normal));
    float vel_tangent_len = lg_vec3_len(vel_tangent);
    
    lg_vec3_t impulse_friction = lg_vec3_zero();
    if (vel_tangent_len > 1e-6f) {
        lg_vec3_t tangent = lg_vec3_scale(vel_tangent, -1.0f / vel_tangent_len);
        
        /* Compute effective mass along tangent */
        float eff_mass_tangent = inv_mass_a + inv_mass_b;
        
        if (body_a) {
            lg_vec3_t r_cross_t = lg_vec3_cross(r_a, tangent);
            eff_mass_tangent += r_cross_t.x * r_cross_t.x * inv_inertia_a.x
                             + r_cross_t.y * r_cross_t.y * inv_inertia_a.y
                             + r_cross_t.z * r_cross_t.z * inv_inertia_a.z;
        }
        if (body_b) {
            lg_vec3_t r_cross_t = lg_vec3_cross(r_b, tangent);
            eff_mass_tangent += r_cross_t.x * r_cross_t.x * inv_inertia_b.x
                             + r_cross_t.y * r_cross_t.y * inv_inertia_b.y
                             + r_cross_t.z * r_cross_t.z * inv_inertia_b.z;
        }
        
        if (eff_mass_tangent > 1e-10f) {
            /* Friction magnitude: clamped by Coulomb model */
            float j_friction = vel_tangent_len / eff_mass_tangent;
            float max_friction = contact->friction * j_normal;
            j_friction = fminf(j_friction, max_friction);
            
            impulse_friction = lg_vec3_scale(tangent, j_friction);
        }
    }
    
    /*=== APPLY IMPULSES ===*/
    lg_vec3_t total_impulse = lg_vec3_add(impulse_normal, impulse_friction);
    
    /* Linear velocity changes 
     * Standard convention: impulse points from B to A (along normal)
     * Body A gets pushed along normal (+= impulse)
     * Body B gets pushed opposite to normal (-= impulse)
     */
    if (body_a) {
        body_a->velocity = lg_vec3_add(body_a->velocity, 
            lg_vec3_scale(total_impulse, inv_mass_a));
    }
    if (body_b) {
        body_b->velocity = lg_vec3_sub(body_b->velocity, 
            lg_vec3_scale(total_impulse, inv_mass_b));
    }
    
    /* Angular velocity changes: tau = r x impulse, delta_omega = I^-1 * tau */
    if (body_a) {
        lg_vec3_t torque_a = lg_vec3_cross(r_a, total_impulse);
        lg_vec3_t delta_omega_a = lg_vec3_mul(torque_a, inv_inertia_a);
        body_a->angular_velocity = lg_vec3_add(body_a->angular_velocity, delta_omega_a);
    }
    if (body_b) {
        lg_vec3_t torque_b = lg_vec3_cross(r_b, total_impulse);
        lg_vec3_t delta_omega_b = lg_vec3_mul(torque_b, inv_inertia_b);
        body_b->angular_velocity = lg_vec3_sub(body_b->angular_velocity, delta_omega_b);
    }
    
    /*=== POSITIONAL CORRECTION ===*/
    float percent = 0.4f; /* Penetration percentage to correct */
    float slop = 0.01f;   /* Penetration allowance */
    float correction_mag = fmaxf(contact->penetration - slop, 0.0f) / 
                          (inv_mass_a + inv_mass_b) * percent;
    /* Correction is along normal (from B toward A)
     * Body A should move along normal (+= correction)
     * Body B should move opposite to normal (-= correction) */
    lg_vec3_t correction = lg_vec3_scale(contact->normal, correction_mag);
    
    if (body_a) trans_a->position = lg_vec3_add(trans_a->position, 
        lg_vec3_scale(correction, inv_mass_a));
    if (body_b) trans_b->position = lg_vec3_sub(trans_b->position, 
        lg_vec3_scale(correction, inv_mass_b));
}

/*============================================================================
 * Broad Phase (Simple Brute Force)
 *===========================================================================*/

static inline void lg_broad_phase(lg_world_t* world, lg_contact_t* contacts, int max_contacts, int* out_count) {
    *out_count = 0;
    
    lg_storage_t* s = &world->storage;
    

    
    for (size_t i = 0; i < s->count && *out_count < max_contacts; i++) {
        for (size_t j = i + 1; j < s->count && *out_count < max_contacts; j++) {
            lg_collider_t* ca = &s->colliders[i];
            lg_collider_t* cb = &s->colliders[j];
            lg_transform_t* ta = &s->transforms[i];
            lg_transform_t* tb = &s->transforms[j];
            
            /* Skip trigger-only pairs for now */
            if (ca->is_trigger || cb->is_trigger) continue;
            
            lg_contact_t contact = {0};
            contact.entity_a = s->entities[i];
            contact.entity_b = s->entities[j];
            
            /* Material combination */
            lg_material_t* ma = &s->materials[i];
            lg_material_t* mb = &s->materials[j];
            contact.restitution = fminf(ma->restitution, mb->restitution);
            contact.friction = sqrtf(ma->friction * mb->friction);
            
            bool collided = false;
            
            /* Sphere-sphere */
            if (ca->type == LG_SHAPE_SPHERE && cb->type == LG_SHAPE_SPHERE) {
                collided = lg_collide_spheres(
                    ta->position, ca->sphere.radius + ca->margin,
                    tb->position, cb->sphere.radius + cb->margin,
                    &contact
                );
            }
            /* Sphere-plane (plane-sphere handled by reversing) */
            else if (ca->type == LG_SHAPE_SPHERE && cb->type == LG_SHAPE_PLANE) {
                collided = lg_collide_sphere_plane(
                    ta->position, ca->sphere.radius + ca->margin,
                    fabsf(tb->rotation.x) < 1e-6f && fabsf(tb->rotation.z) < 1e-6f ? 
                        lg_vec3(0, 1, 0) : lg_quat_rotate(tb->rotation, lg_vec3_up()),
                    cb->plane.distance,
                    &contact
                );
            }
            /* Plane-sphere (swap entities for consistent normal direction) */
            else if (ca->type == LG_SHAPE_PLANE && cb->type == LG_SHAPE_SPHERE) {
                collided = lg_collide_sphere_plane(
                    tb->position, cb->sphere.radius + cb->margin,
                    fabsf(ta->rotation.x) < 1e-6f && fabsf(ta->rotation.z) < 1e-6f ? 
                        lg_vec3(0, 1, 0) : lg_quat_rotate(ta->rotation, lg_vec3_up()),
                    ca->plane.distance,
                    &contact
                );
                /* The contact now has sphere as A, plane as B, normal points from sphere to plane (up)
                 * We want: A = sphere (dynamic), B = plane (static), normal from B to A (up)
                 * The collision function sets normal from sphere to plane (up), which is what we want!
                 * Just need to swap entity IDs.
                 */
                contact.entity_a = s->entities[j];  /* Sphere becomes A */
                contact.entity_b = s->entities[i];  /* Plane becomes B */
            }
            /* Box-plane */
            else if (ca->type == LG_SHAPE_BOX && cb->type == LG_SHAPE_PLANE) {
                collided = lg_collide_box_plane(
                    ta->position, ca->box.half_extents,
                    cb->plane.normal, cb->plane.distance,
                    &contact
                );

            }
            /* Plane-box (swap entities) */
            else if (ca->type == LG_SHAPE_PLANE && cb->type == LG_SHAPE_BOX) {
                contact.entity_a = s->entities[j];  /* Swap: box becomes A */
                contact.entity_b = s->entities[i];  /* Plane becomes B */
                collided = lg_collide_box_plane(
                    tb->position, cb->box.half_extents,
                    ca->plane.normal, ca->plane.distance,
                    &contact
                );
                /* Collision function returns normal = plane_normal (points from plane toward box)
                 * With A=box, B=plane, we want normal from B to A, which is what we have.
                 * No need to flip.
                 */

            }
            /* Sphere-box */
            else if (ca->type == LG_SHAPE_SPHERE && cb->type == LG_SHAPE_BOX) {
                collided = lg_collide_sphere_box(
                    ta->position, ca->sphere.radius + ca->margin,
                    tb->position, cb->box.half_extents,
                    &contact
                );
            }
            /* Box-sphere (swap entities) */
            else if (ca->type == LG_SHAPE_BOX && cb->type == LG_SHAPE_SPHERE) {
                contact.entity_a = s->entities[j];  /* Swap: sphere becomes A */
                contact.entity_b = s->entities[i];  /* Box becomes B */
                collided = lg_collide_sphere_box(
                    tb->position, cb->sphere.radius + cb->margin,
                    ta->position, ca->box.half_extents,
                    &contact
                );
            }
            /* Box-box */
            else if (ca->type == LG_SHAPE_BOX && cb->type == LG_SHAPE_BOX) {
                collided = lg_collide_box_box(
                    ta->position, ca->box.half_extents,
                    tb->position, cb->box.half_extents,
                    &contact
                );
            }
            /* Sphere-capsule */
            else if (ca->type == LG_SHAPE_SPHERE && cb->type == LG_SHAPE_CAPSULE) {
                collided = lg_collide_sphere_capsule(
                    ta->position, ca->sphere.radius + ca->margin,
                    tb->position, lg_quat_rotate(tb->rotation, lg_vec3_up()),
                    cb->capsule.half_height, cb->capsule.radius + cb->margin,
                    &contact
                );
            }
            /* Capsule-sphere (swap entities) */
            else if (ca->type == LG_SHAPE_CAPSULE && cb->type == LG_SHAPE_SPHERE) {
                contact.entity_a = s->entities[j];  /* Swap: sphere becomes A */
                contact.entity_b = s->entities[i];  /* Capsule becomes B */
                collided = lg_collide_sphere_capsule(
                    tb->position, cb->sphere.radius + cb->margin,
                    ta->position, lg_quat_rotate(ta->rotation, lg_vec3_up()),
                    ca->capsule.half_height, ca->capsule.radius + ca->margin,
                    &contact
                );
            }
            /* Capsule-capsule */
            else if (ca->type == LG_SHAPE_CAPSULE && cb->type == LG_SHAPE_CAPSULE) {
                collided = lg_collide_capsule_capsule(
                    ta->position, lg_quat_rotate(ta->rotation, lg_vec3_up()),
                    ca->capsule.half_height, ca->capsule.radius + ca->margin,
                    tb->position, lg_quat_rotate(tb->rotation, lg_vec3_up()),
                    cb->capsule.half_height, cb->capsule.radius + cb->margin,
                    &contact
                );
            }
            /* Sphere-cylinder */
            else if (ca->type == LG_SHAPE_SPHERE && cb->type == LG_SHAPE_CYLINDER) {
                collided = lg_collide_sphere_cylinder(
                    ta->position, ca->sphere.radius + ca->margin,
                    tb->position, lg_quat_rotate(tb->rotation, lg_vec3_up()),
                    cb->cylinder.half_height, cb->cylinder.radius + cb->margin,
                    &contact
                );
            }
            /* Cylinder-sphere (swap) */
            else if (ca->type == LG_SHAPE_CYLINDER && cb->type == LG_SHAPE_SPHERE) {
                contact.entity_a = s->entities[j];
                contact.entity_b = s->entities[i];
                collided = lg_collide_sphere_cylinder(
                    tb->position, cb->sphere.radius + cb->margin,
                    ta->position, lg_quat_rotate(ta->rotation, lg_vec3_up()),
                    ca->cylinder.half_height, ca->cylinder.radius + ca->margin,
                    &contact
                );
            }

            if (collided) {
                contacts[(*out_count)++] = contact;
            }
        }
    }
}

/*============================================================================
 * Simulation Step
 *===========================================================================*/

static inline void lg_world_step(lg_world_t* world) {
    if (!world) return;
    
    float dt = world->config.time_step;
    lg_storage_t* s = &world->storage;
    
    /* Store previous transforms for interpolation */
    if (s->count > 0) {
        memcpy(s->prev_transforms, s->transforms, s->count * sizeof(lg_transform_t));
    }
    
    /* Apply gravity */
    if (world->use_gravity) {
        for (size_t i = 0; i < s->count; i++) {
            lg_body_t* b = &s->bodies[i];
            if (b->type == LG_BODY_DYNAMIC) {
                lg_vec3_t gravity_force = lg_vec3_scale(world->gravity, b->mass);
                b->force = lg_vec3_add(b->force, gravity_force);
            }
        }
    }
    
    /* Integrate motion */
    for (size_t i = 0; i < s->count; i++) {
        lg_integrate(world->integrator, &s->bodies[i], &s->transforms[i], 
                     &s->verlet_states[i], dt);
    }
    
    /* Collision detection and response */
    #define MAX_CONTACTS 1024
    lg_contact_t contacts[MAX_CONTACTS];
    int contact_count;
    
    lg_broad_phase(world, contacts, MAX_CONTACTS, &contact_count);
    
    /* Invoke collision callbacks once per step */
    if (world->collision_callback) {
        for (int i = 0; i < contact_count; i++) {
            world->collision_callback(world, contacts[i].entity_a, contacts[i].entity_b,
                                      world->collision_callback_user_data);
        }
    }
    
    /* Iterative contact resolution */
    for (int iter = 0; iter < world->config.collision_iterations; iter++) {
        for (int i = 0; i < contact_count; i++) {
            lg_resolve_contact(world, &contacts[i]);
        }
    }
    
    /* Clear forces */
    for (size_t i = 0; i < s->count; i++) {
        lg_body_clear_forces(&s->bodies[i]);
    }
    
    /* Update time */
    world->time += dt;
    world->step_count++;
}

/* Variable timestep - steps multiple fixed steps */
static inline void lg_world_update(lg_world_t* world, float dt) {
    if (!world) return;
    
    world->accumulator += dt;
    float step = world->config.time_step;
    
    /* Clamp to prevent spiral of death */
    if (world->accumulator > 0.25f) {
        world->accumulator = 0.25f;
    }
    
    while (world->accumulator >= step) {
        lg_world_step(world);
        world->accumulator -= step;
    }
}

/* Get interpolation factor for visual smoothing */
static inline float lg_world_get_alpha(lg_world_t* world) {
    if (!world) return 0.0f;
    return world->accumulator / world->config.time_step;
}

/* Get interpolated transform for rendering */
static inline lg_transform_t lg_get_interpolated_transform(
    lg_world_t* world, 
    lg_entity_t entity
) {
    int idx = lg_storage_find(&world->storage, entity);
    if (idx < 0) return lg_transform();
    
    float alpha = lg_world_get_alpha(world);
    lg_transform_t* curr = &world->storage.transforms[idx];
    lg_transform_t* prev = &world->storage.prev_transforms[idx];
    
    lg_transform_t result;
    result.position = lg_vec3_lerp(prev->position, curr->position, alpha);
    result.rotation = lg_quat_slerp(prev->rotation, curr->rotation, alpha);
    result.scale = lg_vec3_lerp(prev->scale, curr->scale, alpha);
    
    return result;
}

#ifdef __cplusplus
}
#endif

#endif /* LAGRANGE_SIM_H */
