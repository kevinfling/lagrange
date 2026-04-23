/*
 * Lagrange Example: Bouncing Ball
 * Simple rigid body physics demonstration
 * 
 * Build: gcc -o bouncing_ball bouncing_ball.c -I../include -lm
 */

#define LAGRANGE_IMPLEMENTATION
#include "lagrange.h"

#include <stdio.h>
#include <math.h>

int main(void) {
    printf("=== Lagrange: Bouncing Ball ===\n\n");
    
    /* Create physics world with standard Earth gravity */
    lg_world_config_t config = LG_WORLD_CONFIG_DEFAULT;
    config.time_step = 1.0f / 60.0f;
    config.gravity[0] = 0.0f;
    config.gravity[1] = -9.81f;
    config.gravity[2] = 0.0f;
    
    lg_world_t* world = lg_world_create(&config);
    if (!world) {
        printf("Failed to create world!\n");
        return 1;
    }
    
    /* Create a ball */
    lg_entity_t ball = lg_entity_create(world);
    
    /* Set up transform (start at y=10) */
    lg_set_position(world, ball, lg_vec3(0.0f, 10.0f, 0.0f));
    lg_set_rotation(world, ball, lg_quat_identity());
    
    /* Set up rigid body */
    lg_body_t body = lg_body(1.0f);  /* 1 kg mass */
    body.linear_damping = 0.01f;
    lg_set_body(world, ball, &body);
    
    /* Set up collider (sphere with radius 0.5m) */
    lg_collider_t collider = lg_collider_sphere(0.5f);
    lg_set_collider(world, ball, &collider);
    
    /* Set up material */
    lg_material_t material = LG_MATERIAL_DEFAULT;
    material.restitution = 0.8f;
    lg_set_material(world, ball, &material);
    
    /* Create ground plane */
    lg_entity_t ground = lg_entity_create(world);
    lg_set_position(world, ground, lg_vec3(0.0f, 0.0f, 0.0f));
    
    lg_body_t ground_body = lg_body_static();
    lg_set_body(world, ground, &ground_body);
    
    lg_collider_t ground_collider = lg_collider_plane(lg_vec3(0.0f, 1.0f, 0.0f), 0.0f);
    lg_set_collider(world, ground, &ground_collider);
    
    /* Simulate */
    printf("Simulating bouncing ball...\n");
    printf("Initial: y=10m, v=0, mass=1kg, radius=0.5m, restitution=0.8\n\n");
    
    printf("Time (s) | Height (m) | Velocity (m/s)\n");
    printf("---------|------------|----------------\n");
    
    float max_height = 10.0f;
    int bounces = 0;
    
    for (int frame = 0; frame < 600; frame++) {  /* 10 seconds at 60 FPS */
        float time = frame * config.time_step;
        
        /* Step physics */
        lg_world_step(world);
        
        /* Get state */
        lg_vec3_t pos = lg_get_position(world, ball);
        lg_vec3_t vel = lg_get_velocity(world, ball);
        
        /* Track max height after bounces */
        if (vel.y > 0.0f && pos.y > max_height) {
            max_height = pos.y;
        }
        
        /* Detect bounce (velocity changes from negative to positive) */
        static float prev_vel_y = 0.0f;
        if (prev_vel_y < -0.1f && vel.y > 0.1f) {
            bounces++;
        }
        prev_vel_y = vel.y;
        
        /* Print every 60 frames (1 second) */
        if (frame % 60 == 0) {
            printf("%8.2f | %10.3f | (%6.2f, %6.2f, %6.2f)\n",
                   time, pos.y, vel.x, vel.y, vel.z);
        }
        
        /* Stop if ball comes to rest */
        if (frame > 300 && fabsf(vel.y) < 0.01f && pos.y < 0.51f) {
            printf("\nBall came to rest at t=%.2fs\n", time);
            break;
        }
    }
    
    printf("\n=== Results ===\n");
    printf("Total bounces: %d\n", bounces);
    printf("Max height: %.2fm\n", max_height);
    printf("Final height: %.3fm\n", lg_get_position(world, ball).y);
    
    /* Cleanup */
    lg_world_destroy(world);
    
    return 0;
}
