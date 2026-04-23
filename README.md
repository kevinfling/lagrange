# Lagrange

A header-only C library for physics simulation, featuring rigid body dynamics and orbital mechanics.

## Features

- **Header-only**: Just include and go - no build step needed
- **Pure C11**: Clean C API with no external dependencies
- **Modular**: Include only what you need
- **Efficient**: Small types passed by value, optimized math operations
- **Two precision modes**: Float for games, double for orbital mechanics

### Physics Features

- Rigid body dynamics with forces, torques, and impulses
- Multiple integration methods (Euler, Symplectic, Verlet, RK4)
- Collision detection (sphere, box, plane)
- Sleeping system for performance
- Newtonian gravity with J2 perturbation
- Orbital mechanics calculations

## Quick Start

```c
#define LAGRANGE_IMPLEMENTATION
#include <lagrange.h>
#include <stdio.h>

int main(void) {
    // Create physics world
    lg_world_t* world = lg_world_create(NULL);
    
    // Create a ball
    lg_entity_t ball = lg_entity_create(world);
    lg_set_position(world, ball, lg_vec3(0.0f, 10.0f, 0.0f));
    
    lg_body_t body = lg_body(1.0f);  // 1 kg
    lg_set_body(world, ball, &body);
    
    lg_collider_t collider = lg_collider_sphere(0.5f);
    lg_set_collider(world, ball, &collider);
    
    // Simulate
    for (int i = 0; i < 600; i++) {
        lg_world_step(world);
        
        lg_vec3_t pos = lg_get_position(world, ball);
        printf("t=%.2f: y=%.2f\n", i * world->config.time_step, pos.y);
    }
    
    lg_world_destroy(world);
    return 0;
}
```

## Installation

### Header-only (Recommended)

Just copy the `include/` directory to your project:

```bash
cp -r include/lagrange /path/to/your/project/include/
cp include/lagrange.h /path/to/your/project/include/
```

### Using CMake

```bash
mkdir build && cd build
cmake ..
cmake --build .
sudo cmake --install .
```

Then in your CMakeLists.txt:
```cmake
find_package(lagrange REQUIRED)
target_link_library(your_app lagrange::lagrange m)
```

## Building Examples and Tests

```bash
mkdir build && cd build
cmake .. -DLAGRANGE_BUILD_EXAMPLES=ON -DLAGRANGE_BUILD_TESTS=ON
make

# Run tests
ctest

# Run examples
./bouncing_ball
./orbit
./stacking
```

## API Overview

### Math (`lagrange/math.h`)

```c
lg_vec3_t v = lg_vec3(1.0f, 2.0f, 3.0f);
lg_vec3_t w = lg_vec3_add(v, lg_vec3(4.0f, 5.0f, 6.0f));
float len = lg_vec3_len(v);
lg_vec3_t n = lg_vec3_norm(v);

lg_quat_t q = lg_quat_from_axis_angle(lg_vec3_up(), M_PI / 2.0f);
lg_vec3_t rotated = lg_quat_rotate(q, v);
```

### World (`lagrange/world.h`)

```c
lg_world_t* world = lg_world_create(NULL);  /* Default config */
lg_world_t* world = lg_world_create(&(lg_world_config_t){
    .time_step = 1.0f / 120.0f,
    .gravity = {0.0f, -9.81f, 0.0f}
});

lg_entity_t entity = lg_entity_create(world);
lg_entity_destroy(world, entity);
```

### Components

```c
/* Transform */
lg_set_position(world, entity, lg_vec3(0.0f, 10.0f, 0.0f));
lg_vec3_t pos = lg_get_position(world, entity);
lg_set_rotation(world, entity, lg_quat_identity());

/* Rigid Body */
lg_body_t body = lg_body(1.0f);  /* Dynamic body, 1 kg */
lg_body_t static_body = lg_body_static();
lg_set_body(world, entity, &body);

lg_apply_force(world, entity, lg_vec3(0.0f, 100.0f, 0.0f));
lg_apply_impulse(world, entity, lg_vec3(0.0f, 10.0f, 0.0f));

/* Collider */
lg_collider_t sphere = lg_collider_sphere(0.5f);
lg_collider_t box = lg_collider_box(1.0f, 2.0f, 3.0f);  /* Half extents */
lg_collider_t capsule = lg_collider_capsule(0.5f, 2.0f);
lg_set_collider(world, entity, &sphere);

/* Material */
lg_material_t mat = {
    .restitution = 0.5f,  /* Bounciness */
    .friction = 0.7f
};
lg_set_material(world, entity, &mat);
```

### Simulation

```c
/* Fixed timestep */
lg_world_step(world);

/* Variable timestep (interpolates) */
lg_world_update(world, delta_time);

/* Get interpolation factor for smooth rendering */
float alpha = lg_world_get_alpha(world);
```

### Orbital Mechanics (`lagrange/gravity.h`)

```c
/* Calculate orbital velocity for circular orbit */
double v = lg_orbital_velocity_circular(LG_MU_EARTH, 7000e3);

/* Calculate gravitational acceleration */
lg_vec3_t accel = lg_gravity_accel(satellite_pos, earth_pos, LG_MU_EARTH);

/* Apply to body */
lg_vec3_t force = lg_vec3_scale(accel, body_mass);
lg_apply_force(world, satellite, force);
```

## Naming Conventions

All public symbols use the `lg_` prefix:

| Type | Example |
|------|---------|
| Types | `lg_vec3_t`, `lg_body_t`, `lg_world_t` |
| Functions | `lg_vec3_add()`, `lg_world_step()` |
| Constants | `LG_VEC3_ZERO`, `LG_MU_EARTH` |

## Module Structure

```
include/lagrange.h          # Main umbrella header
include/lagrange/
  ├── types.h               # Entity handles, config
  ├── math.h                # vec3, vec3d, quat, mat4
  ├── transform.h           # Position, rotation, scale
  ├── body.h                # Rigid body dynamics
  ├── collider.h            # Collision shapes
  ├── integrator.h          # Integration methods
  ├── gravity.h             # Orbital mechanics
  ├── world.h               # World and entity management
  └── sim.h                 # Physics step and collision
```

Include individual headers for faster compilation, or use the umbrella header for convenience.

## Performance Tips

1. **Use stack allocation**: Small types (vec3, quat) are designed for pass-by-value
2. **Batch operations**: Process multiple entities in loops for cache efficiency
3. **Use sleeping**: Enable `allow_sleep` on bodies that come to rest
4. **Fixed timestep**: Use `lg_world_step()` for consistent physics
5. **Interpolation**: Use `lg_world_get_alpha()` for smooth visual rendering

## License

MIT License - see LICENSE file for details.

## Contributing

Contributions are welcome! Please follow the existing code style:
- Snake_case for all identifiers
- `lg_` prefix for all public symbols
- Types end with `_t`
- Static inline for header-only implementations
