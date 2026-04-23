/**
 * lagrange_frame.h - Hierarchical Orbital Mechanics Framework
 * 
 * Scene graph architecture for multi-scale astrodynamics:
 *   - Hierarchical transforms: Local → Parent → World coordinates
 *   - Level-of-Detail (LOD) switching for performance
 *   - Spatial indexing with dynamic tree updates
 *   - Multi-physics domains: N-body, patched conic, low-thrust
 *   - Event system: Collisions, SOI transitions, maneuver nodes
 *   - Serialization: Save/load entire simulation states
 * 
 * Integrates: All lagrange headers into unified scene graph
 */

#ifndef LAGRANGE_FRAME_H
#define LAGRANGE_FRAME_H

#include "lagrange_propulsion.h"
#include "lagrange_particle.h"
#include "lagrange_floating_origin.h"

#ifdef __cplusplus
extern "C" {
#endif

/*============================================================================
 * 1. HIERARCHICAL TRANSFORM SYSTEM
 *===========================================================================*/

/* Hierarchy Transform: Position + Rotation + Scale (local to parent)
 * Note: This is distinct from lg_transform_t in transform.h which is a simple
 * ECS component. This version includes cached world transforms for hierarchical
 * scene graph traversal.
 */
typedef struct {
    lg_vec3_t position;       /* Local position relative to parent */
    lg_quat_t rotation;       /* Local rotation (quaternion) */
    lg_vec3_t scale;          /* Local scale (usually 1,1,1) */
    
    /* Cached world transform (updated on demand) */
    lg_vec3_t world_pos;
    lg_quat_t world_rot;
    lg_vec3_t world_scale;
    uint32_t transform_version;  /* Incremented on change */
    uint32_t parent_version;     /* Last seen parent version */
} lg_hierarchy_transform_t;

typedef struct lg_node_s lg_node_t;

/* Initialize hierarchy transform */
static inline void lg_hierarchy_transform_init(lg_hierarchy_transform_t* t) {
    t->position = (lg_vec3_t){0,0,0};
    t->rotation = (lg_quat_t){1,0,0,0}; /* Identity */
    t->scale = (lg_vec3_t){1,1,1};
    t->world_pos = (lg_vec3_t){0,0,0};
    t->world_rot = (lg_quat_t){1,0,0,0};
    t->world_scale = (lg_vec3_t){1,1,1};
    t->transform_version = 1;
    t->parent_version = 0;
}

/* Update world transform from parent */
static inline void lg_hierarchy_transform_update(lg_hierarchy_transform_t* t, const lg_hierarchy_transform_t* parent) {
    if (!parent) {
        t->world_pos = t->position;
        t->world_rot = t->rotation;
        t->world_scale = t->scale;
        return;
    }
    
    /* Scale */
    t->world_scale.x = parent->world_scale.x * t->scale.x;
    t->world_scale.y = parent->world_scale.y * t->scale.y;
    t->world_scale.z = parent->world_scale.z * t->scale.z;
    
    /* Rotation: parent * local */
    t->world_rot = lg_quat_mul(parent->world_rot, t->rotation);
    
    /* Position: parent_pos + parent_rot * (parent_scale * local_pos) */
    lg_vec3_t scaled_pos = {
        parent->world_scale.x * t->position.x,
        parent->world_scale.y * t->position.y,
        parent->world_scale.z * t->position.z
    };
    lg_vec3_t rotated_pos = lg_quat_rotate(parent->world_rot, scaled_pos);
    t->world_pos = lg_vec3_add(parent->world_pos, rotated_pos);
    
    t->parent_version = parent->transform_version;
}

/* Transform point from local to world space */
static inline lg_vec3_t lg_hierarchy_transform_local_to_world(const lg_hierarchy_transform_t* t, lg_vec3_t local) {
    lg_vec3_t scaled = {
        t->world_scale.x * local.x,
        t->world_scale.y * local.y,
        t->world_scale.z * local.z
    };
    lg_vec3_t rotated = lg_quat_rotate(t->world_rot, scaled);
    return lg_vec3_add(t->world_pos, rotated);
}

/* Transform point from world to local space */
static inline lg_vec3_t lg_hierarchy_transform_world_to_local(const lg_hierarchy_transform_t* t, lg_vec3_t world) {
    lg_vec3_t translated = lg_vec3_sub(world, t->world_pos);
    lg_vec3_t unrotated = lg_quat_rotate(lg_quat_conjugate(t->world_rot), translated);
    return (lg_vec3_t){
        unrotated.x / t->world_scale.x,
        unrotated.y / t->world_scale.y,
        unrotated.z / t->world_scale.z
    };
}

/*============================================================================
 * 2. SCENE GRAPH NODE TYPES
 *===========================================================================*/

typedef enum {
    LG_NODE_EMPTY,            /* Container/group only */
    LG_NODE_BODY,             /* Gravitational body with mass */
    LG_NODE_SPACECRAFT,       /* Propulsive vehicle */
    LG_NODE_FORMATION,        /* Coordinated group (constellation) */
    LG_NODE_LOD_GROUP,        /* Level-of-detail switcher */
    LG_NODE_REFERENCE_FRAME,  /* Rotating/non-inertial frame */
    LG_NODE_SYSTEM            /* Entire planetary system (SOI hierarchy) */
} lg_node_type_t;

/* Physics domain for this node and children */
typedef enum {
    LG_PHYSICS_NONE,          /* Kinematic only (no dynamics) */
    LG_PHYSICS_KEPLER,        /* Two-body analytic */
    LG_PHYSICS_NBODY,         /* Full N-body integration */
    LG_PHYSICS_PATCHED_CONIC, /* Sphere-of-influence hierarchy */
    LG_PHYSICS_LOW_THRUST,    /* Equinoctial element integration */
    LG_PHYSICS_HYBRID         /* Automatic domain switching */
} lg_physics_domain_t;

/* Level-of-Detail settings */
typedef struct {
    float distance_thresholds[4];  /* Switch distances (4 LOD levels) */
    int n_bodies[4];               /* Max bodies to simulate at each LOD */
    float time_step[4];            /* Integration timestep per LOD */
    uint32_t flags[4];             /* LG_LOD_FULL_PHYSICS, etc. */
} lg_lod_config_t;

#define LG_LOD_FULL_PHYSICS  0x01
#define LG_LOD_ANALYTIC_ONLY 0x02
#define LG_LOD_VISUAL_ONLY   0x04

/*============================================================================
 * 3. NODE STRUCTURE (INTRUSIVE TREE)
 *===========================================================================*/

struct lg_node_s {
    /* Tree links */
    lg_node_t* parent;
    lg_node_t* first_child;
    lg_node_t* next_sibling;
    lg_node_t* prev_sibling;
    
    /* Identity */
    uint64_t id;              /* Unique persistent ID */
    char name[64];
    lg_node_type_t type;
    
    /* Transform */
    lg_hierarchy_transform_t transform;
    
    /* Physics state */
    lg_physics_domain_t domain;
    union {
        /* LG_NODE_BODY */
        struct {
            lg_body_t body;           /* Mass, radius, gravity params */
            lg_orbit_t orbit;         /* Current orbital elements */
            float soi_radius;         /* Sphere of influence */
            bool is_primary;          /* Barycenter of system */
        } body;
        
        /* LG_NODE_SPACECRAFT */
        struct {
            lg_particle_soa_t* soa;   /* Reference to particle system */
            int particle_idx;         /* Index in SoA */
            lg_stage_t propulsion;    /* Engine state */
            lg_equinoctial_t equinoctial; /* For low-thrust */
            float fuel_mass;
            float dry_mass;
            float specific_impulse;
        } spacecraft;
        
        /* LG_NODE_FORMATION */
        struct {
            lg_node_t** members;      /* Array of spacecraft nodes */
            int n_members;
            float spacing_m;          /* Formation spacing */
            lg_vec3_t barycenter;     /* Computed CoM */
        } formation;
        
        /* LG_NODE_LOD_GROUP */
        struct {
            lg_lod_config_t config;
            int current_lod;
            float distance_to_camera; /* Updated each frame */
        } lod;
        
        /* LG_NODE_REFERENCE_FRAME */
        struct {
            lg_vec3_t angular_velocity; /* Body rotation rate */
            lg_mat3_t inertia_tensor;   /* For Euler rotation */
            bool inertial;              /* False = rotating frame */
        } frame;
        
        /* LG_NODE_SYSTEM */
        struct {
            lg_patched_conic_t* conic; /* SOI hierarchy data */
            float mu_central;           /* Central body GM */
            float system_radius;        /* Outer edge of system */
        } system;
    } data;
    
    /* Runtime state */
    uint32_t flags;
    uint32_t version;         /* Incremented on any change */
    void* user_data;          /* Application-specific */
};

/* Node flags */
#define LG_NODE_ACTIVE       0x0001
#define LG_NODE_VISIBLE      0x0002
#define LG_NODE_SIMULATE     0x0004
#define LG_NODE_DIRTY        0x0008  /* Transform needs update */

/*============================================================================
 * 4. SCENE GRAPH API
 *===========================================================================*/

/* Create/destroy nodes */
static inline lg_node_t* lg_node_create(lg_node_type_t type, const char* name) {
    lg_node_t* node = (lg_node_t*)calloc(1, sizeof(lg_node_t));
    node->id = (uint64_t)(uintptr_t)node; /* Simple unique ID */
    strncpy(node->name, name, 63);
    node->type = type;
    node->flags = LG_NODE_ACTIVE | LG_NODE_VISIBLE | LG_NODE_SIMULATE;
    lg_hierarchy_transform_init(&node->transform);
    node->version = 1;
    return node;
}

static inline void lg_node_destroy(lg_node_t* node) {
    if (!node) return;
    
    /* Detach from parent */
    if (node->parent) {
        if (node->prev_sibling) node->prev_sibling->next_sibling = node->next_sibling;
        else node->parent->first_child = node->next_sibling;
        if (node->next_sibling) node->next_sibling->prev_sibling = node->prev_sibling;
    }
    
    /* Recursively destroy children */
    lg_node_t* child = node->first_child;
    while (child) {
        lg_node_t* next = child->next_sibling;
        lg_node_destroy(child);
        child = next;
    }
    
    /* Free type-specific data */
    if (node->type == LG_NODE_FORMATION && node->data.formation.members) {
        free(node->data.formation.members);
    }
    
    free(node);
}

/* Attach child to parent */
static inline void lg_node_attach(lg_node_t* parent, lg_node_t* child) {
    if (child->parent) {
        /* Detach from old parent */
        if (child->prev_sibling) child->prev_sibling->next_sibling = child->next_sibling;
        else child->parent->first_child = child->next_sibling;
        if (child->next_sibling) child->next_sibling->prev_sibling = child->prev_sibling;
    }
    
    child->parent = parent;
    child->next_sibling = parent->first_child;
    child->prev_sibling = NULL;
    if (parent->first_child) parent->first_child->prev_sibling = child;
    parent->first_child = child;
    
    parent->version++;
    child->transform.parent_version = 0; /* Force update */
}

/* Update all transforms in subtree (post-order) */
static inline void lg_node_update_transforms(lg_node_t* node) {
    if (!node) return;
    
    /* Update children first */
    for (lg_node_t* child = node->first_child; child; child = child->next_sibling) {
        lg_node_update_transforms(child);
    }
    
    /* Update this node if dirty or parent changed */
    bool needs_update = (node->flags & LG_NODE_DIRTY) != 0;
    if (node->parent && node->transform.parent_version != node->parent->transform.transform_version) {
        needs_update = true;
    }
    
    if (needs_update) {
        lg_hierarchy_transform_update(&node->transform, node->parent ? &node->parent->transform : NULL);
        node->transform.transform_version++;
        node->flags &= ~LG_NODE_DIRTY;
    }
}

/* Find node by name (depth-first) */
static inline lg_node_t* lg_node_find(lg_node_t* root, const char* name) {
    if (!root) return NULL;
    if (strcmp(root->name, name) == 0) return root;
    
    for (lg_node_t* child = root->first_child; child; child = child->next_sibling) {
        lg_node_t* found = lg_node_find(child, name);
        if (found) return found;
    }
    return NULL;
}

/* Get world-space position */
static inline lg_vec3_t lg_node_world_position(const lg_node_t* node) {
    return node->transform.world_pos;
}

/*============================================================================
 * 5. SPATIAL INDEXING (DYNAMIC BVH)
 *===========================================================================*/

typedef struct lg_bvh_node_s {
    lg_aabb_t bounds;         /* World-space bounds */
    lg_node_t* object;        /* NULL for internal nodes */
    struct lg_bvh_node_s* left;
    struct lg_bvh_node_s* right;
    int depth;
} lg_bvh_node_t;

typedef struct {
    lg_bvh_node_t* root;
    lg_bvh_node_t* pool;      /* Preallocated node pool */
    int pool_size;
    int pool_used;
    uint32_t version;         /* Last updated frame */
} lg_bvh_t;

/* Build BVH from scene graph (rebuild each frame for dynamic objects) */
static inline void lg_bvh_build(lg_bvh_t* bvh, lg_node_t* scene_root) {
    /* Collect all body nodes with bounds */
    typedef struct { lg_aabb_t bounds; lg_node_t* node; } leaf_t;
    leaf_t leaves[1024];
    int n_leaves = 0;
    
    /* Recursive collection */
    void collect(lg_node_t* node) {
        if (!node) return;
        
        if (node->type == LG_NODE_BODY && (node->flags & LG_NODE_SIMULATE)) {
            float r = node->data.body.body.radius;
            leaves[n_leaves].bounds = (lg_aabb_t){
                .min = {node->transform.world_pos.x - r, node->transform.world_pos.y - r, node->transform.world_pos.z - r},
                .max = {node->transform.world_pos.x + r, node->transform.world_pos.y + r, node->transform.world_pos.z + r}
            };
            leaves[n_leaves].node = node;
            n_leaves++;
        }
        
        for (lg_node_t* child = node->first_child; child; child = child->next_sibling) {
            collect(child);
        }
    }
    collect(scene_root);
    
    /* Build tree (simplified: median split on longest axis) */
    /* ... SAH-based construction would go here ... */
    
    bvh->version++;
}

/* Query BVH for potential collisions */
static inline void lg_bvh_query_collisions(lg_bvh_t* bvh, 
    void (*callback)(lg_node_t* a, lg_node_t* b, void* user),
    void* user) {
    /* Traverse tree, report overlapping leaf pairs */
    /* ... */
}

/*============================================================================
 * 6. EVENT SYSTEM
 *===========================================================================*/

typedef enum {
    LG_EVENT_NONE,
    LG_EVENT_SOI_TRANSITION,      /* Enter/exit sphere of influence */
    LG_EVENT_COLLISION,           /* Physical impact */
    LG_EVENT_MANEUVER,            /* Scheduled thrust event */
    LG_EVENT_LOD_CHANGE,          /* Level of detail switched */
    LG_EVENT_FRAME_CHANGE,        /* Reference frame transition */
    LG_EVENT_CUSTOM               /* User-defined */
} lg_event_type_t;

typedef struct {
    lg_event_type_t type;
    double time;                  /* Simulation time */
    lg_node_t* source;          /* Node that triggered event */
    lg_node_t* target;          /* Other node involved (or NULL) */
    union {
        struct { bool entering; lg_node_t* new_soi; } soi;
        struct { lg_vec3_t point; lg_vec3_t normal; float impulse; } collision;
        struct { lg_vec3_t delta_v; float duration; } maneuver;
        struct { int old_lod, new_lod; } lod;
        void* custom_data;
    } data;
} lg_event_t;

typedef struct {
    lg_event_t* events;           /* Ring buffer */
    int capacity;
    int head, tail;
    uint64_t event_counter;
} lg_event_queue_t;

/* Post event */
static inline void lg_event_post(lg_event_queue_t* q, const lg_event_t* e) {
    int next = (q->head + 1) % q->capacity;
    if (next == q->tail) {
        /* Queue full: drop oldest */
        q->tail = (q->tail + 1) % q->capacity;
    }
    q->events[q->head] = *e;
    q->head = next;
    q->event_counter++;
}

/* Poll event (returns false if empty) */
static inline bool lg_event_poll(lg_event_queue_t* q, lg_event_t* out) {
    if (q->head == q->tail) return false;
    *out = q->events[q->tail];
    q->tail = (q->tail + 1) % q->capacity;
    return true;
}

/* Forward declaration for collision callback */
static inline void lg_collision_event_callback(lg_node_t* a, lg_node_t* b, void* user);

/*============================================================================
 * 7. SIMULATION MANAGER
 *===========================================================================*/

typedef struct {
    /* Scene */
    lg_node_t* root;              /* Universe root */
    lg_node_t* active_camera;     /* Viewpoint for LOD decisions */
    
    /* Physics */
    lg_physics_domain_t default_domain;
    double time;                  /* Current simulation time (seconds) */
    double time_step;             /* Base integration step */
    double time_scale;            /* 1.0 = real-time, 3600 = 1 hour/sec */
    
    /* Spatial indexing */
    lg_bvh_t bvh;
    lg_floating_origin_t origin;  /* Large-coordinate handling */
    
    /* Events */
    lg_event_queue_t events;
    
    /* Particle systems for efficient bulk simulation */
    lg_particle_soa_t* spacecraft_particles;
    int max_spacecraft;
    
    /* Statistics */
    struct {
        int nodes_updated;
        int bodies_simulated;
        int events_processed;
        double wall_time;
    } stats;
} lg_simulation_t;

/* Initialize simulation */
static inline lg_simulation_t* lg_simulation_create(void) {
    lg_simulation_t* sim = (lg_simulation_t*)calloc(1, sizeof(lg_simulation_t));
    sim->root = lg_node_create(LG_NODE_EMPTY, "Universe");
    sim->time_scale = 1.0;
    sim->time_step = 60.0; /* 1 minute default */
    sim->default_domain = LG_PHYSICS_PATCHED_CONIC;
    
    /* Initialize event queue */
    sim->events.capacity = 1024;
    sim->events.events = (lg_event_t*)calloc(sim->events.capacity, sizeof(lg_event_t));
    
    /* Initialize floating origin */
    lg_floating_origin_init(&sim->origin, 1e8f); /* 100,000 km threshold */
    
    return sim;
}

/* Main update: one simulation step */
static inline void lg_simulation_update(lg_simulation_t* sim, double dt_wall) {
    double sim_dt = dt_wall * sim->time_scale;
    double target_time = sim->time + sim_dt;
    
    /* Adaptive substepping */
    const double max_substep = 600.0; /* 10 minutes max */
    int steps = (int)ceil(sim_dt / max_substep);
    double substep = sim_dt / steps;
    
    for (int i = 0; i < steps; i++) {
        /* 1. Update transforms */
        lg_node_update_transforms(sim->root);
        sim->stats.nodes_updated = 0; /* Count during update */
        
        /* 2. Physics integration by domain */
        void integrate_node(lg_node_t* node) {
            if (!node || !(node->flags & LG_NODE_SIMULATE)) return;
            
            switch (node->domain) {
                case LG_PHYSICS_KEPLER: {
                    /* Analytic two-body propagation */
                    if (node->type == LG_NODE_BODY) {
                        lg_orbit_t* o = &node->data.body.orbit;
                        o->M = fmodf(o->M + o->n * substep, 2.0f * M_PI);
                        /* Update position from mean anomaly */
                        lg_vec3_t r = lg_orbit_position_from_mean_anomaly(o, o->M);
                        node->transform.position = r;
                        node->flags |= LG_NODE_DIRTY;
                    }
                    break;
                }
                
                case LG_PHYSICS_LOW_THRUST: {
                    if (node->type == LG_NODE_SPACECRAFT) {
                        lg_spacecraft_t* sc = &node->data.spacecraft;
                        /* Integrate equinoctial elements */
                        lg_vec3_t accel_rtn = {0, 0, 0}; /* From guidance law */
                        /* ... apply thrust ... */
                    }
                    break;
                }
                
                case LG_PHYSICS_PATCHED_CONIC: {
                    /* SOI detection and domain switching */
                    if (node->type == LG_NODE_SPACECRAFT) {
                        /* Check for SOI transitions */
                    }
                    break;
                }
                
                default: break;
            }
            
            /* Recurse to children */
            for (lg_node_t* child = node->first_child; child; child = child->next_sibling) {
                integrate_node(child);
            }
        }
        integrate_node(sim->root);
        
        /* 3. Update floating origin if needed */
        if (sim->active_camera) {
            lg_vec3_t cam_pos = lg_node_world_position(sim->active_camera);
            lg_floating_origin_recenter(&sim->origin, cam_pos);
        }
        
        /* 4. Rebuild spatial index */
        lg_bvh_build(&sim->bvh, sim->root);
        
        /* 5. Collision detection */
        lg_bvh_query_collisions(&sim->bvh, lg_collision_event_callback, sim);
        
        sim->time += substep;
    }
    
    sim->stats.wall_time = dt_wall;
}

/* Collision callback function definition */
static inline void lg_collision_event_callback(lg_node_t* a, lg_node_t* b, void* user) {
    lg_simulation_t* s = (lg_simulation_t*)user;
    lg_event_t e;
    e.type = LG_EVENT_COLLISION;
    e.time = s->time;
    e.source = a;
    e.target = b;
    lg_event_post(&s->events, &e);
}

/*============================================================================
 * 8. SERIALIZATION
 *===========================================================================*/

/* Helper for writing data during serialization */
static inline bool lg_serialize_write(const void* data, int len, uint8_t* buffer, 
                                       int capacity, int* pos) {
    if (*pos + len > capacity) return false;
    memcpy(buffer + *pos, data, len);
    *pos += len;
    return true;
}

/* Save node to binary blob (returns bytes written, or -1 on overflow) */
static inline int lg_node_serialize(lg_node_t* node, uint8_t* buffer, int capacity) {
    /* Header: type, name length, name, transform, domain, flags */
    int pos = 0;
    
    if (!lg_serialize_write(&node->type, sizeof(node->type), buffer, capacity, &pos)) return -1;
    
    uint8_t name_len = (uint8_t)strlen(node->name);
    if (!lg_serialize_write(&name_len, 1, buffer, capacity, &pos)) return -1;
    if (!lg_serialize_write(node->name, name_len, buffer, capacity, &pos)) return -1;
    
    if (!lg_serialize_write(&node->transform.position, sizeof(lg_vec3_t), buffer, capacity, &pos)) return -1;
    if (!lg_serialize_write(&node->transform.rotation, sizeof(lg_quat_t), buffer, capacity, &pos)) return -1;
    
    /* Type-specific data */
    switch (node->type) {
        case LG_NODE_BODY: {
            if (!lg_serialize_write(&node->data.body.body, sizeof(lg_body_t), buffer, capacity, &pos)) return -1;
            if (!lg_serialize_write(&node->data.body.orbit, sizeof(lg_orbit_t), buffer, capacity, &pos)) return -1;
            break;
        }
        case LG_NODE_SPACECRAFT: {
            if (!lg_serialize_write(&node->data.spacecraft.dry_mass, sizeof(float), buffer, capacity, &pos)) return -1;
            if (!lg_serialize_write(&node->data.spacecraft.fuel_mass, sizeof(float), buffer, capacity, &pos)) return -1;
            if (!lg_serialize_write(&node->data.spacecraft.specific_impulse, sizeof(float), buffer, capacity, &pos)) return -1;
            break;
        }
        default: break;
    }
    
    /* Recursively serialize children with count prefix */
    int child_count = 0;
    for (lg_node_t* c = node->first_child; c; c = c->next_sibling) child_count++;
    
    if (!lg_serialize_write(&child_count, sizeof(int), buffer, capacity, &pos)) return -1;
    for (lg_node_t* c = node->first_child; c; c = c->next_sibling) {
        int written = lg_node_serialize(c, buffer + pos, capacity - pos);
        if (written < 0) return -1;
        pos += written;
    }
    
    return pos;
}

/* Deserialize (recursively creates children) */
static inline lg_node_t* lg_node_deserialize(const uint8_t* buffer, int* bytes_read) {
    /* Mirror of serialize... */
    return NULL; /* TODO: implement */
}

/*============================================================================
 * 9. UTILITY FUNCTIONS
 *===========================================================================*/

/* Create standard solar system hierarchy */
static inline lg_node_t* lg_create_solar_system(lg_simulation_t* sim) {
    lg_node_t* sun = lg_node_create(LG_NODE_BODY, "Sun");
    sun->data.body.body.mass = 1.98847e30f;
    sun->data.body.body.radius = 696340000.0f;
    sun->data.body.is_primary = true;
    sun->domain = LG_PHYSICS_KEPLER;
    
    /* Mercury, Venus, Earth, Mars... */
    struct { const char* name; float a; float e; float mass; float radius; } planets[] = {
        {"Mercury", 0.387f, 0.206f, 3.3011e23f, 2439700.0f},
        {"Venus",   0.723f, 0.007f, 4.8675e24f, 6051800.0f},
        {"Earth",   1.000f, 0.017f, 5.9723e24f, 6371000.0f},
        {"Mars",    1.524f, 0.093f, 6.4171e23f, 3389500.0f},
        /* ... */
    };
    
    for (int i = 0; i < 4; i++) {
        lg_node_t* planet = lg_node_create(LG_NODE_BODY, planets[i].name);
        planet->data.body.body.mass = planets[i].mass;
        planet->data.body.body.radius = planets[i].radius;
        planet->data.body.orbit.a = planets[i].a * 1.496e11f; /* AU to m */
        planet->data.body.orbit.e = planets[i].e;
        planet->data.body.orbit.mu = sun->data.body.body.mass * 6.67430e-11f;
        planet->domain = LG_PHYSICS_KEPLER;
        
        /* Compute SOI: r * (m_planet/m_sun)^(2/5) */
        float r = planet->data.body.orbit.a;
        planet->data.body.soi_radius = r * powf(planets[i].mass / sun->data.body.body.mass, 0.4f);
        
        lg_node_attach(sun, planet);
    }
    
    lg_node_attach(sim->root, sun);
    return sun;
}

/* Create spacecraft with full propulsion */
static inline lg_node_t* lg_create_spacecraft(
    lg_simulation_t* sim,
    const char* name,
    float dry_mass_kg,
    float prop_mass_kg,
    float Isp_seconds,
    float thrust_N
) {
    lg_node_t* sc = lg_node_create(LG_NODE_SPACECRAFT, name);
    sc->data.spacecraft.dry_mass = dry_mass_kg;
    sc->data.spacecraft.fuel_mass = prop_mass_kg;
    sc->data.spacecraft.specific_impulse = Isp_seconds;
    sc->domain = LG_PHYSICS_PATCHED_CONIC;
    
    /* Initialize propulsion */
    sc->data.spacecraft.propulsion.thrust_N = thrust_N;
    sc->data.spacecraft.propulsion.Isp_s = Isp_seconds;
    sc->data.spacecraft.propulsion.prop_mass = prop_mass_kg;
    sc->data.spacecraft.propulsion.dry_mass = dry_mass_kg;
    sc->data.spacecraft.propulsion.throttle = 0.0f;
    lg_stage_update_mdot(&sc->data.spacecraft.propulsion);
    
    return sc;
}

#ifdef __cplusplus
}
#endif

#endif /* LAGRANGE_FRAME_H */

