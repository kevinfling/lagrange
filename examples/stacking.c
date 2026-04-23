/*
 * Lagrange Example: Stacking
 * Demonstrates collision and stacking stability
 * 
 * Build: gcc -o stacking stacking.c -I../include -lm
 */

#define LAGRANGE_IMPLEMENTATION
#include "lagrange.h"

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

/* Random float between min and max */
float randf(float min, float max) {
    return min + (float)rand() / (float)RAND_MAX * (max - min);
}

int main(void) {
    srand((unsigned)time(NULL));
    
    printf("=== Lagrange: Stacking ===\n\n");
    
    /* Create physics world */
    lg_world_config_t config = LG_WORLD_CONFIG_DEFAULT;
    config.time_step = 1.0f / 60.0f;
    config.velocity_iterations = 8;
    config.position_iterations = 3;
    config.collision_iterations = 4;
    
    lg_world_t* world = lg_world_create(&config);
    if (!world) {
        printf("Failed to create world!\n");
        return 1;
    }
    
    /* Create ground */
    lg_entity_t ground = lg_entity_create(world);
    lg_body_t ground_body = lg_body_static();
    lg_set_body(world, ground, &ground_body);
    lg_collider_t ground_collider = lg_collider_plane(lg_vec3(0.0f, 1.0f, 0.0f), 0.0f);
    lg_set_collider(world, ground, &ground_collider);
    
    /* Create a stack of boxes */
    printf("Creating stack of boxes...\n\n");
    
    int num_boxes = 10;
    float box_size = 1.0f;
    
    lg_entity_t boxes[20];
    
    for (int i = 0; i < num_boxes && i < 20; i++) {
        boxes[i] = lg_entity_create(world);
        
        /* Stack vertically with slight offset for realism */
        float x_offset = randf(-0.02f, 0.02f);
        float z_offset = randf(-0.02f, 0.02f);
        float y_pos = 0.5f * box_size + i * box_size;  /* Half height + stack */
        
        lg_set_position(world, boxes[i], lg_vec3(x_offset, y_pos, z_offset));
        
        /* Create box body */
        lg_body_t body = lg_body(1.0f);  /* 1 kg each */
        body.linear_damping = 0.1f;
        body.angular_damping = 0.1f;
        lg_set_body(world, boxes[i], &body);
        
        /* Create box collider */
        lg_collider_t collider = lg_collider_box(box_size * 0.5f, box_size * 0.5f, box_size * 0.5f);
        lg_set_collider(world, boxes[i], &collider);
        
        /* Set material */
        lg_material_t material = LG_MATERIAL_DEFAULT;
        material.restitution = 0.1f;
        material.friction = 0.7f;
        lg_set_material(world, boxes[i], &material);
    }
    
    printf("Stack of %d boxes created.\n", num_boxes);
    printf("Box size: %.2fm x %.2fm x %.2fm\n", box_size, box_size, box_size);
    printf("Simulating %d frames (%.1f seconds)...\n\n", 600, 600.0f / 60.0f);
    
    printf("Frame | Top Box Y | Stack Height | Status\n");
    printf("------|-----------|--------------|--------\n");
    
    /* Simulate */
    int stable_frames = 0;
    float prev_height = 0.0f;
    
    for (int frame = 0; frame < 600; frame++) {
        /* Step physics */
        lg_world_step(world);
        
        /* Get top box position */
        lg_vec3_t top_pos = lg_get_position(world, boxes[num_boxes - 1]);
        float stack_height = top_pos.y + box_size * 0.5f;
        
        /* Check stability */
        if (frame > 0 && fabsf(stack_height - prev_height) < 0.001f) {
            stable_frames++;
        } else {
            stable_frames = 0;
        }
        prev_height = stack_height;
        
        /* Print every 60 frames */
        if (frame % 60 == 0) {
            const char* status = "Settling";
            if (stable_frames > 60) status = "Stable";
            else if (frame < 60) status = "Initial";
            
            printf("%5d | %9.3f | %12.3f | %s\n",
                   frame, top_pos.y, stack_height, status);
        }
        
        /* Check if stack fell over */
        bool stack_fell = false;
        for (int i = 0; i < num_boxes; i++) {
            lg_vec3_t pos = lg_get_position(world, boxes[i]);
            if (fabsf(pos.x) > 5.0f || fabsf(pos.z) > 5.0f) {
                stack_fell = true;
                break;
            }
        }
        
        if (stack_fell) {
            printf("\n*** Stack collapsed at frame %d! ***\n", frame);
            break;
        }
    }
    
    /* Final statistics */
    printf("\n=== Final Stack State ===\n");
    float expected_height = num_boxes * box_size;
    lg_vec3_t final_top = lg_get_position(world, boxes[num_boxes - 1]);
    float final_height = final_top.y + box_size * 0.5f;
    
    printf("Expected height: %.2fm\n", expected_height);
    printf("Final height:    %.2fm\n", final_height);
    printf("Height error:    %.2f%%\n", 100.0f * fabsf(final_height - expected_height) / expected_height);
    
    /* Show each box position */
    printf("\nBox positions (from bottom to top):\n");
    for (int i = 0; i < num_boxes; i++) {
        lg_vec3_t pos = lg_get_position(world, boxes[i]);
        printf("  Box %2d: (%.3f, %.3f, %.3f)\n", i, pos.x, pos.y, pos.z);
    }
    
    /* Cleanup */
    lg_world_destroy(world);
    
    return 0;
}
