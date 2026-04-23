/*
 * Lagrange Example: Earth Orbit
 * Orbital mechanics demonstration using Velocity Verlet integration
 * 
 * Build: gcc -o orbit orbit.c -I../include -lm
 */

#define LAGRANGE_IMPLEMENTATION
#include "lagrange.h"

#include <stdio.h>
#include <math.h>

int main(void) {
    printf("=== Lagrange: Earth Orbit ===\n\n");
    
    /* Use smaller timestep for orbital stability */
    float dt = 0.1f;  /* 0.1 second per step */
    
    /* Create physics world */
    lg_world_config_t config = LG_WORLD_CONFIG_DEFAULT;
    config.time_step = dt;
    config.gravity[0] = 0.0f;
    config.gravity[1] = 0.0f;
    config.gravity[2] = 0.0f;
    config.enable_collision = false;
    
    lg_world_t* world = lg_world_create(&config);
    if (!world) {
        printf("Failed to create world!\n");
        return 1;
    }
    
    /* Use semi-implicit Euler (more stable for this simple implementation) */
    lg_world_set_integrator(world, LG_INTEGRATOR_SEMI_IMPLICIT_EULER);
    
    /* Create satellite in Low Earth Orbit (400 km altitude) */
    double altitude = 400e3;
    double orbit_radius = LG_RADIUS_EARTH + altitude;
    double orbital_vel = lg_orbital_velocity_circular(LG_MU_EARTH, orbit_radius);
    double period = lg_orbital_period(LG_MU_EARTH, orbit_radius);
    
    /* Calculate steps per orbit */
    int steps_per_orbit = (int)(period / dt);
    
    printf("Orbit Parameters:\n");
    printf("  Altitude: %.1f km\n", altitude / 1000.0);
    printf("  Orbit radius: %.1f km\n", orbit_radius / 1000.0);
    printf("  Orbital velocity: %.3f km/s\n", orbital_vel / 1000.0);
    printf("  Orbital period: %.1f minutes (%.2f hours)\n", 
           period / 60.0, period / 3600.0);
    printf("  Integration: Semi-implicit Euler with %.1fs timestep (%d steps/orbit)\n", 
           dt, steps_per_orbit);
    printf("\n");
    
    lg_entity_t satellite = lg_entity_create(world);
    lg_set_position(world, satellite, lg_vec3((float)orbit_radius, 0.0f, 0.0f));
    
    lg_body_t sat_body = lg_body(1000.0f);
    sat_body.velocity = lg_vec3(0.0f, (float)orbital_vel, 0.0f);
    sat_body.linear_damping = 0.0f;  /* No damping for orbital mechanics */
    sat_body.angular_damping = 0.0f;
    lg_set_body(world, satellite, &sat_body);
    
    /* Simulate for 3 orbits */
    int total_steps = steps_per_orbit * 3;
    
    printf("Simulating %d orbits (%d steps)...\n\n", 3, total_steps);
    printf("Orbit | Day | Hour | Distance (km) | Speed (km/s) | Error (m)\n");
    printf("------|-----|------|---------------|--------------|----------\n");
    
    double max_error = 0.0;
    double sum_error = 0.0;
    int samples = 0;
    
    for (int step = 0; step <= total_steps; step++) {
        /* Apply Earth's gravity */
        lg_vec3_t sat_pos = lg_get_position(world, satellite);
        lg_vec3_t earth_pos = lg_vec3_zero();
        
        lg_body_t* sat = lg_get_body(world, satellite);
        if (sat) {
            lg_vec3_t accel = lg_gravity_accel(sat_pos, earth_pos, LG_MU_EARTH);
            lg_vec3_t force = lg_vec3_scale(accel, sat->mass);
            lg_apply_force(world, satellite, force);
        }
        
        /* Step physics */
        lg_world_step(world);
        
        /* Track statistics */
        sat_pos = lg_get_position(world, satellite);
        double radius = lg_vec3_len(sat_pos);
        double error = fabs(radius - orbit_radius);
        
        max_error = fmax(max_error, error);
        sum_error += error;
        samples++;
        
        /* Print at each orbit completion */
        if (step % steps_per_orbit == 0) {
            lg_vec3_t sat_vel = lg_get_velocity(world, satellite);
            double speed = lg_vec3_len(sat_vel);
            double hours = step * dt / 3600.0;
            int days = (int)(hours / 24.0);
            int hour = (int)(hours) % 24;
            int orbit = step / steps_per_orbit;
            
            printf("  %2d  | %3d | %4d | %13.3f | %12.3f | %8.1f\n",
                   orbit, days, hour, radius / 1000.0, speed / 1000.0, error);
        }
    }
    
    printf("\n=== Results ===\n");
    printf("Expected orbital radius: %.3f km\n", orbit_radius / 1000.0);
    printf("Max radius deviation:    %.3f m (%.4f%%)\n", 
           max_error, (max_error / orbit_radius) * 100.0);
    printf("Average radius error:    %.3f m\n", sum_error / samples);
    
    if (max_error < 1000.0) {
        printf("\nOrbit is stable! ✓\n");
    } else {
        printf("\nOrbit drift detected - consider smaller timestep\n");
    }
    
    lg_world_destroy(world);
    return 0;
}
