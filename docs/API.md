# LaGrange API Reference

Complete API documentation for the LaGrange physics library.

## Table of Contents

1. [Overview](#overview)
2. [Core Modules](#core-modules)
3. [Simulation](#simulation)
4. [Orbital Mechanics](#orbital-mechanics)
5. [Advanced Modules](#advanced-modules)
6. [Mathematical Reference](#mathematical-reference)

---

## Overview

LaGrange is a header-only C11 physics library providing:

- **Rigid body dynamics** with forces, torques, and impulses
- **Multiple integrators**: Euler, Symplectic Euler, Velocity Verlet, RK4
- **Collision detection**: Sphere, Box, Capsule, Plane shapes
- **Orbital mechanics**: Newtonian gravity with perturbations
- **Advanced features**: N-body, wavelet analysis, fractional calculus, Koopman operators

### Architecture

```
Entity (lg_entity_t) - lightweight handle
    ├── Transform (position, rotation, scale)
    ├── Body (mass, velocity, forces)
    ├── Collider (shape, bounds)
    └── Material (friction, restitution)
```

### Precision Model

| Domain | Type | Use Case |
|--------|------|----------|
| Game Physics | `float` / `lg_vec3_t` | Real-time simulation |
| Orbital Mechanics | `double` / `lg_vec3d_t` | High-precision propagation |

---

## Core Modules

### math.h - Vector Mathematics

#### lg_vec3_t - 3D Vector (float)

```c
typedef struct {
    float x, y, z;
} lg_vec3_t;
```

**Construction**

| Function | Description |
|----------|-------------|
| `lg_vec3(x, y, z)` | Construct from components |
| `lg_vec3_zero()` | Zero vector (0, 0, 0) |
| `lg_vec3_one()` | One vector (1, 1, 1) |
| `lg_vec3_up/right/forward()` | Unit axis vectors |

**Basic Operations**

| Function | Formula | Description |
|----------|---------|-------------|
| `lg_vec3_add(a, b)` | $\mathbf{a} + \mathbf{b}$ | Component-wise addition |
| `lg_vec3_sub(a, b)` | $\mathbf{a} - \mathbf{b}$ | Component-wise subtraction |
| `lg_vec3_scale(v, s)` | $s \cdot \mathbf{v}$ | Scalar multiplication |
| `lg_vec3_mul(a, b)` | $\mathbf{a} \odot \mathbf{b}$ | Component-wise product |
| `lg_vec3_neg(v)` | $-\mathbf{v}$ | Negation |

**Vector Products**

| Function | Formula | Description |
|----------|---------|-------------|
| `lg_vec3_dot(a, b)` | $\mathbf{a} \cdot \mathbf{b} = \sum_i a_i b_i$ | Dot product |
| `lg_vec3_cross(a, b)` | $\mathbf{a} \times \mathbf{b}$ | Cross product |

**Cross Product Formula:**
$$\mathbf{a} \times \mathbf{b} = \begin{bmatrix} a_y b_z - a_z b_y \\ a_z b_x - a_x b_z \\ a_x b_y - a_y b_x \end{bmatrix}$$

**Length & Normalization**

| Function | Formula | Description |
|----------|---------|-------------|
| `lg_vec3_len_sq(v)` | $\|\mathbf{v}\|^2 = \mathbf{v} \cdot \mathbf{v}$ | Squared length |
| `lg_vec3_len(v)` | $\|\mathbf{v}\| = \sqrt{\mathbf{v} \cdot \mathbf{v}}$ | Length |
| `lg_vec3_norm(v)` | $\hat{\mathbf{v}} = \mathbf{v} / \|\mathbf{v}\|$ | Unit vector |
| `lg_vec3_dist(a, b)` | $\|\mathbf{a} - \mathbf{b}\|$ | Distance |

**Interpolation**

| Function | Formula | Description |
|----------|---------|-------------|
| `lg_vec3_lerp(a, b, t)` | $\mathbf{a} + t(\mathbf{b} - \mathbf{a})$ | Linear interpolation |

#### lg_vec3d_t - 3D Vector (double)

Double-precision variant for orbital mechanics. Same function names with `d` suffix:
- `lg_vec3d()`, `lg_vec3d_add()`, `lg_vec3d_len()`, etc.

#### lg_quat_t - Quaternion

```c
typedef struct {
    float x, y, z, w;
} lg_quat_t;
```

Quaternions represent rotations in 3D space: $q = x\mathbf{i} + y\mathbf{j} + z\mathbf{k} + w$

**Construction**

| Function | Description |
|----------|-------------|
| `lg_quat(x, y, z, w)` | From components |
| `lg_quat_identity()` | Identity (0, 0, 0, 1) |
| `lg_quat_from_axis_angle(axis, angle)` | Axis-angle to quaternion |
| `lg_quat_from_to(from, to)` | Rotation between vectors |

**Hamilton Product**

For $q_1 = (x_1, y_1, z_1, w_1)$ and $q_2 = (x_2, y_2, z_2, w_2)$:

$$q_1 \otimes q_2 = \begin{bmatrix}
w_1 x_2 + x_1 w_2 + y_1 z_2 - z_1 y_2 \\
w_1 y_2 - x_1 z_2 + y_1 w_2 + z_1 x_2 \\
w_1 z_2 + x_1 y_2 - y_1 x_2 + z_1 w_2 \\
w_1 w_2 - x_1 x_2 - y_1 y_2 - z_1 z_2
\end{bmatrix}$$

**Rotation**

To rotate vector $\mathbf{v}$ by quaternion $q$:
$$\mathbf{v}' = q \otimes (0, \mathbf{v}) \otimes q^*$$

where $q^* = (-x, -y, -z, w)$ is the conjugate.

**Spherical Linear Interpolation (SLERP)**

$$\text{slerp}(q_1, q_2, t) = \frac{\sin((1-t)\theta)}{\sin\theta} q_1 + \frac{\sin(t\theta)}{\sin\theta} q_2$$

where $\theta = \arccos(q_1 \cdot q_2)$

---

### types.h - Core Types

#### lg_entity_t

```c
typedef uint64_t lg_entity_t;
#define LG_ENTITY_INVALID ((lg_entity_t)0)
```

Opaque handle to entities in the physics world.

#### lg_world_config_t

```c
typedef struct {
    float gravity[3];        // Default: {0, -9.81, 0}
    float time_step;         // Default: 1/60 s
    int velocity_iterations; // Default: 6
    int position_iterations; // Default: 2
    bool enable_sleeping;    // Default: true
    float sleep_threshold;   // Default: 0.01
    float linear_damping;    // Default: 0.0
    float angular_damping;   // Default: 0.0
    bool enable_collision;   // Default: true
    int collision_iterations;// Default: 4
    size_t max_entities;     // Default: 10000
    size_t max_bodies;       // Default: 5000
    size_t max_colliders;    // Default: 5000
} lg_world_config_t;
```

---

### body.h - Rigid Body Dynamics

#### lg_body_type_t

```c
typedef enum {
    LG_BODY_DYNAMIC,    // Affected by forces
    LG_BODY_KINEMATIC,  // Moved by code, not forces
    LG_BODY_STATIC      // Immovable
} lg_body_type_t;
```

#### lg_body_t

```c
typedef struct {
    float mass;
    float inv_mass;
    lg_vec3_t inertia;      // Principal moments (diagonal)
    lg_vec3_t inv_inertia;
    lg_vec3_t velocity;
    lg_vec3_t force;
    lg_vec3_t angular_velocity;
    lg_vec3_t torque;
    float linear_damping;
    float angular_damping;
    lg_body_type_t type;
    bool allow_sleep;
    bool is_sleeping;
} lg_body_t;
```

**Newton's Second Law:**
$$\mathbf{F} = m\mathbf{a} \quad \Rightarrow \quad \mathbf{a} = \mathbf{F}/m = \mathbf{F} \cdot m_{inv}$$

**Angular Motion:**
$$\boldsymbol{\tau} = \mathbf{I}\boldsymbol{\alpha} \quad \Rightarrow \quad \boldsymbol{\alpha} = \mathbf{I}^{-1}\boldsymbol{\tau}$$

where $\mathbf{I}$ is the inertia tensor (currently diagonal only).

**Kinetic Energy:**
$$KE = \frac{1}{2}m\|\mathbf{v}\|^2 + \frac{1}{2}\boldsymbol{\omega}^T \mathbf{I} \boldsymbol{\omega}$$

**Construction**

| Function | Description |
|----------|-------------|
| `lg_body(mass)` | Dynamic body with given mass |
| `lg_body_static()` | Static body (infinite mass) |
| `lg_body_kinematic()` | Kinematic body |

**Force Application**

| Function | Formula | Description |
|----------|---------|-------------|
| `lg_body_apply_force(b, f)` | $\mathbf{F} \mathrel{+}= \mathbf{f}$ | Accumulate force |
| `lg_body_apply_force_at(b, f, p, c)` | $\mathbf{F} \mathrel{+}= \mathbf{f}$, $\boldsymbol{\tau} \mathrel{+}= (\mathbf{p} - \mathbf{c}) \times \mathbf{f}$ | Force at point |
| `lg_body_apply_torque(b, t)` | $\boldsymbol{\tau} \mathrel{+}= \mathbf{t}$ | Accumulate torque |
| `lg_body_apply_impulse(b, j)` | $\mathbf{v} \mathrel{+}= \mathbf{j}/m$ | Instantaneous velocity change |
| `lg_body_apply_impulse_at(b, j, p, c)` | $\mathbf{v} \mathrel{+}= \mathbf{j}/m$, $\boldsymbol{\omega} \mathrel{+}= \mathbf{I}^{-1}((\mathbf{p} - \mathbf{c}) \times \mathbf{j})$ | Impulse at point |

---

### collider.h - Collision Shapes

#### lg_shape_type_t

```c
typedef enum {
    LG_SHAPE_SPHERE,
    LG_SHAPE_BOX,
    LG_SHAPE_CAPSULE,
    LG_SHAPE_CYLINDER,
    LG_SHAPE_PLANE
} lg_shape_type_t;
```

#### Volume Formulas

| Shape | Volume Formula |
|-------|---------------|
| Sphere ($r$) | $V = \frac{4}{3}\pi r^3$ |
| Box ($h_x, h_y, h_z$ half-extents) | $V = 8h_x h_y h_z$ |
| Capsule ($r$, $h$ half-height) | $V = \frac{4}{3}\pi r^3 + \pi r^2 (2h)$ |
| Cylinder ($r$, $h$ half-height) | $V = \pi r^2 (2h)$ |
| Plane | $\infty$ |

#### Inertia Tensor Formulas

For mass $m$:

| Shape | Principal Moments |
|-------|-------------------|
| Sphere ($r$) | $I_{xx} = I_{yy} = I_{zz} = \frac{2}{5}mr^2$ |
| Box ($w, h, d$) | $I_{xx} = \frac{1}{12}m(h^2 + d^2)$, etc. |
| Capsule | Approximate as cylinder + point masses |
| Cylinder | $I_{xx} = I_{zz} = \frac{1}{12}m(3r^2 + 4h^2)$, $I_{yy} = \frac{1}{2}mr^2$ |

**Construction**

| Function | Parameters | Description |
|----------|------------|-------------|
| `lg_collider_sphere(r)` | radius | Sphere centered at origin |
| `lg_collider_box(hx, hy, hz)` | half-extents | Axis-aligned box |
| `lg_collider_capsule(r, h)` | radius, height | Capsule aligned with Y-axis |
| `lg_collider_cylinder(r, h)` | radius, height | Cylinder aligned with Y-axis |
| `lg_collider_plane(n, d)` | normal, distance | Infinite plane |

---

### transform.h - Spatial Transformations

#### lg_transform_t

```c
typedef struct {
    lg_vec3_t position;
    lg_quat_t rotation;
    lg_vec3_t scale;
} lg_transform_t;
```

**Transform Hierarchy**

Local to world transformation:
$$\mathbf{p}_{world} = \mathbf{t} + \mathbf{q} \otimes (\mathbf{s} \odot \mathbf{p}_{local}) \otimes \mathbf{q}^*$$

World to local (inverse):
$$\mathbf{p}_{local} = (\mathbf{q}^* \otimes (\mathbf{p}_{world} - \mathbf{t})) \oslash \mathbf{s}$$

**Functions**

| Function | Description |
|----------|-------------|
| `lg_transform()` | Identity transform |
| `lg_transform_at(pos)` | Translation only |
| `lg_transform_full(pos, rot, scale)` | All components |
| `lg_transform_point(t, local)` | Local → World |
| `lg_transform_inv_point(t, world)` | World → Local |
| `lg_transform_direction(t, local)` | Direction (no translation) |
| `lg_transform_forward/up/right(t)` | Local axes in world space |

---

## Simulation

### integrator.h - Numerical Integration

#### Integration Methods

**Explicit Euler (First Order)**

$$\mathbf{v}_{n+1} = \mathbf{v}_n + \mathbf{a}_n \Delta t$$
$$\mathbf{x}_{n+1} = \mathbf{x}_n + \mathbf{v}_n \Delta t$$

- Simple but unstable for oscillatory systems
- Energy grows unbounded

**Semi-Implicit (Symplectic) Euler**

$$\mathbf{v}_{n+1} = \mathbf{v}_n + \mathbf{a}_n \Delta t$$
$$\mathbf{x}_{n+1} = \mathbf{x}_n + \mathbf{v}_{n+1} \Delta t$$

- Updates position with new velocity
- Better energy conservation
- Default for game physics

**Velocity Verlet**

$$\mathbf{x}_{n+1} = \mathbf{x}_n + \mathbf{v}_n \Delta t + \frac{1}{2}\mathbf{a}_n \Delta t^2$$
$$\mathbf{v}_{n+1/2} = \mathbf{v}_n + \frac{1}{2}\mathbf{a}_n \Delta t$$
$$\mathbf{a}_{n+1} = \mathbf{F}(\mathbf{x}_{n+1})/m$$
$$\mathbf{v}_{n+1} = \mathbf{v}_{n+1/2} + \frac{1}{2}\mathbf{a}_{n+1} \Delta t$$

- Time-reversible
- Good energy conservation
- Requires storing previous acceleration

**Runge-Kutta 4 (RK4)**

$$k_1 = f(t_n, y_n)$$
$$k_2 = f(t_n + \frac{\Delta t}{2}, y_n + \frac{\Delta t}{2}k_1)$$
$$k_3 = f(t_n + \frac{\Delta t}{2}, y_n + \frac{\Delta t}{2}k_2)$$
$$k_4 = f(t_n + \Delta t, y_n + \Delta t \cdot k_3)$$
$$y_{n+1} = y_n + \frac{\Delta t}{6}(k_1 + 2k_2 + 2k_3 + k_4)$$

- Fourth-order accuracy
- More expensive (4 evaluations per step)
- Good for high-precision needs

#### Usage

```c
lg_world_t* world = lg_world_create(NULL);
lg_world_set_integrator(world, LG_INTEGRATOR_SEMI_IMPLICIT_EULER);
// or LG_INTEGRATOR_VELOCITY_VERLET, LG_INTEGRATOR_RK4
```

---

### world.h - World Management

#### lg_world_t

Container for all physics entities and simulation state.

**Lifecycle**

| Function | Description |
|----------|-------------|
| `lg_world_create(config)` | Create world (NULL config for defaults) |
| `lg_world_destroy(world)` | Free all resources |
| `lg_world_reset(world)` | Clear all entities, reset time |

**Entity Management**

| Function | Description |
|----------|-------------|
| `lg_entity_create(world)` | Create new entity |
| `lg_entity_destroy(world, entity)` | Remove entity |
| `lg_entity_valid(world, entity)` | Check if entity exists |

**Component Access**

| Function | Returns | Description |
|----------|---------|-------------|
| `lg_get_transform(world, e)` | `lg_transform_t*` | Transform component |
| `lg_get_body(world, e)` | `lg_body_t*` | Body component |
| `lg_get_collider(world, e)` | `lg_collider_t*` | Collider component |
| `lg_get_material(world, e)` | `lg_material_t*` | Material component |

**Shortcuts**

| Function | Description |
|----------|-------------|
| `lg_get/set_position(world, e)` | Position accessor |
| `lg_get/set_rotation(world, e)` | Rotation accessor |
| `lg_get/set_velocity(world, e)` | Velocity accessor |
| `lg_get/set_mass(world, e)` | Mass accessor |
| `lg_apply_force(world, e, f)` | Apply force to body |
| `lg_apply_impulse(world, e, j)` | Apply impulse to body |

**Iteration**

```c
void my_callback(lg_world_t* w, lg_entity_t e, void* data) {
    // Process entity e
}

lg_world_for_each(world, my_callback, user_data);
```

---

### sim.h - Simulation Step

#### lg_world_step()

Performs one fixed timestep simulation:

```c
lg_world_step(world);  // Advances by world->config.time_step
```

**Step Order:**

1. **Apply Gravity**: $\mathbf{F}_g = m\mathbf{g}$ for all dynamic bodies
2. **Integrate**: Update velocities and positions using selected integrator
3. **Collision Detection**: Broad phase → Narrow phase
4. **Collision Response**: Impulse-based resolution with position correction
5. **Clear Forces**: Reset force accumulators

#### lg_world_update()

Variable timestep with interpolation:

```c
lg_world_update(world, delta_time);  // Handles sub-stepping
float alpha = lg_world_get_alpha(world);  // Interpolation factor [0,1]
```

#### Collision Detection

**lg_contact_t**

```c
typedef struct {
    lg_entity_t entity_a, entity_b;
    lg_vec3_t point;        // Contact point in world space
    lg_vec3_t normal;       // Points from A to B
    float penetration;      // Positive = overlap
    float restitution;      // Combined bounciness
    float friction;         // Combined friction
} lg_contact_t;
```

**Supported Collisions:**

| Pair | Algorithm | Status |
|------|-----------|--------|
| Sphere-Sphere | Distance test | ✅ |
| Sphere-Plane | Signed distance | ✅ |
| Box-Plane | Vertex projection | ✅ |
| Box-Box | SAT | 🔴 Missing |
| Sphere-Box | Minkowski portal | 🔴 Missing |
| Sphere-Capsule | Line segment distance | 🔴 Missing |

---

## Orbital Mechanics

### gravity.h - Newtonian Gravity

#### Physical Constants

| Constant | Value | Units |
|----------|-------|-------|
| `LG_G` | $6.67430 \times 10^{-11}$ | m³/(kg·s²) |
| `LG_MU_EARTH` | $3.986004418 \times 10^{14}$ | m³/s² |
| `LG_MU_SUN` | $1.32712440018 \times 10^{20}$ | m³/s² |
| `LG_RADIUS_EARTH` | $6.371 \times 10^6$ | m |

#### Newton's Law of Gravitation

Force between two masses:
$$\mathbf{F} = -\frac{G m_1 m_2}{r^2} \hat{\mathbf{r}} = -\frac{\mu m_2}{r^2} \hat{\mathbf{r}}$$

where $\mu = G m_1$ is the standard gravitational parameter.

**Functions**

| Function | Formula | Description |
|----------|---------|-------------|
| `lg_gravity_point()` | $\mathbf{F} = -\frac{\mu m}{r^2}\hat{\mathbf{r}}$ | Gravitational force |
| `lg_gravity_accel()` | $\mathbf{a} = -\frac{\mu}{r^2}\hat{\mathbf{r}}$ | Acceleration (force/mass) |

#### Orbital Velocity

Circular orbit velocity:
$$v_{circ} = \sqrt{\frac{\mu}{r}}$$

Escape velocity:
$$v_{esc} = \sqrt{\frac{2\mu}{r}} = \sqrt{2} \cdot v_{circ}$$

Orbital period (Kepler's 3rd Law):
$$T = 2\pi\sqrt{\frac{a^3}{\mu}}$$

where $a$ is the semi-major axis.

#### J2 Perturbation

Oblateness correction for Earth's equatorial bulge:
$$\mathbf{a}_{J2} = \frac{3}{2}\frac{J_2 \mu R^2}{r^4}\left[\frac{x}{r}\left(5\frac{z^2}{r^2}-1\right), \frac{y}{r}\left(5\frac{z^2}{r^2}-1\right), \frac{z}{r}\left(5\frac{z^2}{r^2}-3\right)\right]$$

where $J_2 \approx 1.08263 \times 10^{-3}$ for Earth.

#### Orbital Elements

Classical orbital elements:

| Element | Symbol | Description |
|---------|--------|-------------|
| Semi-major axis | $a$ | Orbit size |
| Eccentricity | $e$ | Orbit shape (0=circle) |
| Inclination | $i$ | Tilt from reference plane |
| RAAN | $\Omega$ | Longitude of ascending node |
| Arg. of periapsis | $\omega$ | Orientation of ellipse |
| True anomaly | $\nu$ | Position on orbit |

**State Vector ↔ Elements Conversion**

```c
lg_orbital_elements_t elem = lg_state_to_elements(pos, vel, mu);
lg_elements_to_state(&elem, mu, &pos, &vel);
```

---

## Advanced Modules

### particle.h - N-Body Simulation

#### Barnes-Hut Tree

Approximates N-body gravity in $O(N \log N)$ using spatial subdivision:

1. Build octree containing all particles
2. Compute center of mass for each node
3. For each particle, traverse tree:
   - If $\frac{s}{d} < \theta$ (MAC criterion), use node approximation
   - Otherwise, traverse children

where $s$ is node size, $d$ is distance to COM, $\theta \approx 0.5$.

**Status:** 🔴 Mostly stubbed (see TODO.md)

#### Individual Time Steps

Particles use different $\Delta t$ based on local dynamical time:

$$\Delta t_i \sim \eta \frac{|\mathbf{a}_i|}{|\dot{\mathbf{a}}_i|}$$

where $\eta \approx 0.1$ is accuracy parameter.

### wavelet.h - Multi-Resolution Analysis

#### Lifting Scheme

In-place wavelet transform using prediction and update steps:

**Haar Wavelet (Simplest):**
- Split: even/odd samples
- Predict: $d = o - e$ (detail = difference)
- Update: $s = e + d/2$ (smooth = average)

**Multi-Level Decomposition:**
$$f(t) = \sum_k c_{J,k} \phi_{J,k}(t) + \sum_{j=1}^J \sum_k d_{j,k} \psi_{j,k}(t)$$

where $\phi$ is scaling function, $\psi$ is wavelet.

**Status:** ⚠️ Partial (Haar only fully implemented)

### fractional.h - Fractional Calculus

#### Fractional Integral (Riemann-Liouville)

$$I^\alpha[f](t) = \frac{1}{\Gamma(\alpha)} \int_0^t (t-s)^{\alpha-1} f(s) \, ds$$

for $0 < \alpha < 1$.

**Discrete Approximation:**
$$I^\alpha_n \approx \sum_{k=0}^{n-1} b_{n-k} f_k$$

where $b_j = \frac{\Delta t^\alpha}{\Gamma(1+\alpha)}[(j+1)^\alpha - j^\alpha]$

#### Applications

- **Anomalous diffusion**: Non-Brownian particle motion
- **Tidal dissipation**: Viscoelastic memory effects
- **Yarkovsky effect**: Thermal inertia on asteroids

### koopman.h - Koopman Operator Theory

#### Extended Dynamic Mode Decomposition (EDMD)

Lifts nonlinear dynamics to linear operator in observable space:

$$\boldsymbol{\psi}(\mathbf{x}_{n+1}) = \mathbf{K} \boldsymbol{\psi}(\mathbf{x}_n)$$

where $\boldsymbol{\psi}: \mathbb{R}^n \to \mathbb{C}^D$ is observable function and $\mathbf{K}$ is the Koopman operator.

**DCT Observables:**
$$\psi_k(\mathbf{x}) = \exp(i \omega_k \cdot \mathbf{x})$$

**Learning:**
$$\mathbf{K} = \boldsymbol{\Psi}_Y \boldsymbol{\Psi}_X^+$$

where $^+$ denotes pseudoinverse.

**Status:** 🔴 Mostly stubbed (see TODO.md)

### singularity.h - Black Hole Physics

#### Schwarzschild Metric

Spacetime interval:
$$ds^2 = -\left(1-\frac{r_s}{r}\right)c^2 dt^2 + \left(1-\frac{r_s}{r}\right)^{-1} dr^2 + r^2 d\Omega^2$$

where $r_s = \frac{2GM}{c^2}$ is Schwarzschild radius.

**Key Radii:**

| Radius | Location | Significance |
|--------|----------|--------------|
| $r_s$ | Event horizon | Point of no return |
| $3r_s/2$ | Photon sphere | Unstable light orbit |
| $3r_s$ | ISCO | Innermost stable circular orbit |

**Gravitational Redshift:**
$$z = \frac{1}{\sqrt{1 - r_s/r}} - 1$$

**Status:** 🔴 Geodesic integration stubbed (see TODO.md)

### exoplanet.h - Procedural Generation

#### Gaussian Copula Sampling

Generate correlated planetary parameters:

1. Sample $\mathbf{z} \sim \mathcal{N}(0, \mathbf{I})$
2. Transform: $\mathbf{x} = \mathbf{L}\mathbf{z}$ (Cholesky of correlation matrix)
3. CDF transform: $\mathbf{u} = \Phi(\mathbf{x})$
4. Inverse CDF: $\mathbf{y} = F^{-1}(\mathbf{u})$

**Status:** ⚠️ Partial (simplified distributions)

---

## Mathematical Reference

### Common Formulas

#### Vector Identities

Triple product:
$$\mathbf{a} \cdot (\mathbf{b} \times \mathbf{c}) = \mathbf{b} \cdot (\mathbf{c} \times \mathbf{a}) = \mathbf{c} \cdot (\mathbf{a} \times \mathbf{b})$$

Vector triple product:
$$\mathbf{a} \times (\mathbf{b} \times \mathbf{c}) = \mathbf{b}(\mathbf{a} \cdot \mathbf{c}) - \mathbf{c}(\mathbf{a} \cdot \mathbf{b})$$

#### Quaternion Identities

Rotation composition:
$$q_{AB} \otimes q_{BC} = q_{AC}$$

Inverse rotation:
$$q^{-1} = q^* / |q|^2$$

Rotation of vector:
$$\mathbf{v}' = q \otimes (0, \mathbf{v}) \otimes q^*$$

#### Collision Response

**Impulse magnitude** (with restitution $e$):
$$j = \frac{-(1+e)\mathbf{v}_{rel} \cdot \mathbf{n}}{m_1^{-1} + m_2^{-1}}$$

**Position correction** (to prevent sinking):
$$\mathbf{x}_1 \mathrel{{-}{=}} \frac{m_1^{-1}}{m_1^{-1} + m_2^{-1}} \cdot \text{percent} \cdot d \cdot \mathbf{n}$$

where $d$ is penetration depth.

### Constants Reference

| Name | Value | Description |
|------|-------|-------------|
| $\pi$ | 3.14159265 | Circle ratio |
| $G$ | $6.67430 \times 10^{-11}$ | Gravitational constant |
| $c$ | $2.998 \times 10^8$ | Speed of light |
| $M_\odot$ | $1.989 \times 10^{30}$ | Solar mass (kg) |
| $M_\oplus$ | $5.972 \times 10^{24}$ | Earth mass (kg) |
| AU | $1.496 \times 10^{11}$ | Astronomical unit (m) |

---

## Usage Examples

### Basic Rigid Body Simulation

```c
#define LAGRANGE_IMPLEMENTATION
#include "lagrange.h"

lg_world_t* world = lg_world_create(NULL);

// Create a falling ball
lg_entity_t ball = lg_entity_create(world);
lg_set_position(world, ball, lg_vec3(0, 10, 0));

lg_body_t body = lg_body(1.0f);  // 1 kg
lg_set_body(world, ball, &body);

lg_collider_t collider = lg_collider_sphere(0.5f);
lg_set_collider(world, ball, &collider);

// Ground plane
lg_entity_t ground = lg_entity_create(world);
lg_body_t ground_body = lg_body_static();
lg_set_body(world, ground, &ground_body);
lg_collider_t ground_collider = lg_collider_plane(lg_vec3_up(), 0);
lg_set_collider(world, ground, &ground_collider);

// Simulate
for (int i = 0; i < 600; i++) {
    lg_world_step(world);
    lg_vec3_t pos = lg_get_position(world, ball);
    printf("t=%.2f: y=%.3f\n", i * world->config.time_step, pos.y);
}

lg_world_destroy(world);
```

### Orbital Mechanics

```c
// Satellite in Low Earth Orbit
lg_vec3_t pos = lg_vec3(7e6f, 0, 0);  // 7000 km from center
lg_vec3_t center = lg_vec3_zero();

// Compute circular orbit velocity
double r = 7e6;
double v_circ = lg_orbital_velocity_circular(LG_MU_EARTH, r);
lg_vec3_t vel = lg_vec3(0, 0, (float)v_circ);

// Apply gravity
lg_vec3_t accel = lg_gravity_accel(pos, center, LG_MU_EARTH);
lg_apply_force(world, satellite, lg_vec3_scale(accel, mass));
```

---

## See Also

- [TODO.md](../TODO.md) - Development status and missing features
- [README.md](../README.md) - Quick start and overview
