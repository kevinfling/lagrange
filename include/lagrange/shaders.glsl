/**
 * lagrange_shaders.glsl - Complete Space Visualization Shader Pack
 * 
 * Shader modules:
 *   - Skybox: Cubemap and equirectangular rendering with tone mapping
 *   - Nebula: Raymarched volumetric fog with procedural noise
 *   - Starfield: Multi-layer parallax with proper motion blur
 *   - Engine trails: GPU particle ribbons with fade decay
 *   - Post-process: Bloom, lens dirt, chromatic warp, film grain
 *   - UI overlay: Orbit paths, SOI spheres, maneuver nodes
 * 
 * Designed for OpenGL 3.3+ / GLSL 330 core
 */

#ifndef LAGRANGE_SHADERS_GLSL
#define LAGRANGE_SHADERS_GLSL

/*============================================================================
 * 1. UNIFORM BLOCKS & SHARED DEFINITIONS
 *===========================================================================*/

layout(std140) uniform CameraBlock {
    mat4 view;
    mat4 projection;
    mat4 view_proj;
    mat4 inv_view;
    mat4 inv_proj;
    vec3 camera_pos;
    float camera_fov;
    vec3 camera_dir;
    float time;
    vec2 resolution;
    float near_plane;
    float far_plane;
};

layout(std140) uniform SceneBlock {
    vec3 sun_direction;
    float sun_intensity;
    vec3 ambient_color;
    float exposure;
    float fog_density;
    float fog_height_falloff;
    vec2 padding;
};

/* Tone mapping operators */
vec3 tonemap_aces(vec3 x) {
    // ACES filmic approximation
    float a = 2.51, b = 0.03, c = 2.43, d = 0.59, e = 0.14;
    return clamp((x * (a * x + b)) / (x * (c * x + d) + e), 0.0, 1.0);
}

vec3 tonemap_reinhard(vec3 x) {
    return x / (1.0 + x);
}

/* Random/hash functions */
float hash(vec2 p) {
    return fract(1e4 * sin(17.0 * p.x + p.y * 0.1) * 
                 (0.1 + abs(sin(p.y * 13.0 + p.x))));
}

float hash3(vec3 p) {
    return fract(sin(dot(p, vec3(12.9898, 78.233, 45.164))) * 43758.5453);
}

/* Value noise */
float noise(vec3 p) {
    vec3 i = floor(p);
    vec3 f = fract(p);
    f = f * f * (3.0 - 2.0 * f); // Smoothstep
    
    float n = i.x + i.y * 57.0 + 113.0 * i.z;
    return mix(mix(mix(hash3(i + vec3(0,0,0)), hash3(i + vec3(1,0,0)), f.x),
                   mix(hash3(i + vec3(0,1,0)), hash3(i + vec3(1,1,0)), f.x), f.y),
               mix(mix(hash3(i + vec3(0,0,1)), hash3(i + vec3(1,0,1)), f.x),
                   mix(hash3(i + vec3(0,1,1)), hash3(i + vec3(1,1,1)), f.x), f.y), f.z);
}

/* FBM for nebula clouds */
float fbm(vec3 p, int octaves) {
    float v = 0.0;
    float a = 0.5;
    vec3 shift = vec3(100.0);
    for (int i = 0; i < octaves; ++i) {
        v += a * noise(p);
        p = p * 2.0 + shift;
        a *= 0.5;
    }
    return v;
}

/* Rotation matrix from axis-angle */
mat3 rotate_axis(vec3 axis, float angle) {
    float c = cos(angle), s = sin(angle);
    float t = 1.0 - c;
    axis = normalize(axis);
    return mat3(
        t * axis.x * axis.x + c,        t * axis.x * axis.y - s * axis.z, t * axis.x * axis.z + s * axis.y,
        t * axis.x * axis.y + s * axis.z, t * axis.y * axis.y + c,          t * axis.y * axis.z - s * axis.x,
        t * axis.x * axis.z - s * axis.y, t * axis.y * axis.z + s * axis.x, t * axis.z * axis.z + c
    );
}

/*============================================================================
 * 2. SKYBOX VERTEX SHADER
 *===========================================================================*/

#ifdef LAGRANGE_SHADER_SKYBOX_VERT
layout(location = 0) in vec3 a_position; // Cube vertices [-1,1]

out vec3 v_direction;

void main() {
    // Remove translation from view matrix
    mat4 rot_view = mat4(mat3(view));
    vec4 clip_pos = projection * rot_view * vec4(a_position, 1.0);
    
    // Force depth to far plane
    gl_Position = clip_pos.xyww;
    v_direction = a_position;
}
#endif

/*============================================================================
 * 3. SKYBOX FRAGMENT SHADER (Equirectangular + Procedural Stars)
 *===========================================================================*/

#ifdef LAGRANGE_SHADER_SKYBOX_FRAG
in vec3 v_direction;
out vec4 frag_color;

uniform sampler2D u_skybox_equirect; // Optional: 2:1 equirectangular map
uniform float u_use_procedural;      // 0 = texture, 1 = procedural

/* Procedural star field */
vec3 render_procedural_stars(vec3 dir) {
    vec3 color = vec3(0.0);
    
    // Multiple octaves of point stars
    for (int layer = 0; layer < 4; layer++) {
        float scale = pow(2.0, float(layer)) * 100.0;
        vec3 p = dir * scale;
        
        // Grid cells
        vec3 cell = floor(p);
        vec3 frac = fract(p);
        
        // Random star in each cell
        float star = hash3(cell);
        if (star > 0.97) {
            // Star position within cell
            vec3 star_pos = cell + vec3(
                hash(cell.xy),
                hash(cell.yz),
                hash(cell.zx)
            );
            
            float dist = length(p - star_pos);
            float brightness = smoothstep(0.15, 0.0, dist);
            
            // Color temperature based on hash
            float temp = hash3(cell + vec3(1.0));
            vec3 star_color;
            if (temp < 0.3) star_color = vec3(0.8, 0.9, 1.0); // Blue-white
            else if (temp < 0.7) star_color = vec3(1.0, 0.95, 0.8); // Yellow-white
            else star_color = vec3(1.0, 0.7, 0.5); // Orange-red
            
            // Twinkle
            float twinkle = 0.8 + 0.2 * sin(time * 3.0 + hash3(cell) * 10.0);
            
            color += star_color * brightness * twinkle * (0.5 / float(layer + 1));
        }
    }
    
    // Milky way band (simplified)
    float galactic_plane = abs(dot(normalize(dir), vec3(0.0, 0.1, 1.0)));
    float milky_way = pow(1.0 - galactic_plane, 3.0) * 0.3;
    color += vec3(0.9, 0.85, 0.95) * milky_way * fbm(dir * 5.0, 4);
    
    // Nebula background (very subtle)
    float neb = fbm(dir * 3.0 + vec3(time * 0.01), 3);
    color += vec3(0.1, 0.15, 0.3) * neb * 0.1;
    
    return color;
}

/* Equirectangular texture lookup */
vec3 sample_equirect(vec3 dir) {
    vec2 uv;
    uv.x = atan(dir.z, dir.x) / (2.0 * 3.14159265) + 0.5;
    uv.y = asin(clamp(dir.y, -1.0, 1.0)) / 3.14159265 + 0.5;
    return texture(u_skybox_equirect, uv).rgb;
}

void main() {
    vec3 dir = normalize(v_direction);
    
    vec3 color;
    if (u_use_procedural > 0.5) {
        color = render_procedural_stars(dir);
    } else {
        color = sample_equirect(dir);
    }
    
    // Tone mapping and gamma correction
    color *= exposure;
    color = tonemap_aces(color);
    color = pow(color, vec3(1.0 / 2.2)); // Gamma
    
    frag_color = vec4(color, 1.0);
}
#endif

/*============================================================================
 * 4. NEBULA RAYMARCHING VERTEX
 *===========================================================================*/

#ifdef LAGRANGE_SHADER_NEBULA_VERT
layout(location = 0) in vec3 a_position;

out vec3 v_world_pos;
out vec3 v_local_pos;

uniform mat4 u_model;

void main() {
    v_local_pos = a_position;
    v_world_pos = (u_model * vec4(a_position, 1.0)).xyz;
    gl_Position = view_proj * vec4(v_world_pos, 1.0);
}
#endif

/*============================================================================
 * 5. NEBULA RAYMARCHING FRAGMENT (Volumetric Fog)
 *===========================================================================*/

#ifdef LAGRANGE_SHADER_NEBULA_FRAG
in vec3 v_world_pos;
in vec3 v_local_pos;
out vec4 frag_color;

uniform vec3 u_nebula_center;
uniform float u_nebula_radius;
uniform vec3 u_nebula_color;
uniform float u_density_scale;

/* SDF for nebula bounds */
float sd_sphere(vec3 p, float r) {
    return length(p) - r;
}

/* Procedural nebula density */
float nebula_density(vec3 p) {
    // Warp space for organic shapes
    vec3 q = p * 0.5;
    float warp = fbm(q + fbm(q * 2.0, 3), 3);
    
    // Base cloud shape
    float d = sd_sphere(p - u_nebula_center, u_nebula_radius);
    
    // Add detail
    float detail = fbm(p * 3.0 + warp, 5) * 0.3;
    
    // Density falls off at edges
    float density = smoothstep(0.5, -0.2, d + detail);
    
    // Internal structure
    density *= (0.5 + 0.5 * fbm(p * 8.0, 4));
    
    return max(0.0, density * u_density_scale);
}

/* Raymarch through volume */
vec4 raymarch_nebula(vec3 ro, vec3 rd, float tmax) {
    const int MAX_STEPS = 64;
    const float STEP_SIZE = u_nebula_radius / 20.0;
    
    vec4 accum = vec4(0.0);
    float t = 0.0;
    
    // Random offset to reduce banding
    t += hash(gl_FragCoord.xy) * STEP_SIZE;
    
    for (int i = 0; i < MAX_STEPS && t < tmax && accum.a < 0.95; i++) {
        vec3 p = ro + rd * t;
        
        float density = nebula_density(p);
        
        if (density > 0.01) {
            // Lighting from sun direction
            float light = max(0.0, dot(normalize(p - u_nebula_center), -sun_direction));
            light = 0.2 + 0.8 * pow(light, 2.0);
            
            // Emission + absorption
            vec3 emission = u_nebula_color * density * light;
            float absorption = density * 0.5;
            
            // Accumulate (back-to-front for simplicity, or front-to-back with proper alpha)
            float alpha = (1.0 - accum.a) * absorption;
            accum.rgb += emission * alpha;
            accum.a += alpha;
        }
        
        // Adaptive step: smaller in dense regions
        float step = STEP_SIZE * (1.0 + density * 2.0);
        t += step;
    }
    
    return accum;
}

void main() {
    vec3 ro = camera_pos;
    vec3 rd = normalize(v_world_pos - camera_pos);
    
    // Ray-sphere intersection for entry/exit
    vec3 oc = ro - u_nebula_center;
    float a = dot(rd, rd);
    float b = 2.0 * dot(oc, rd);
    float c = dot(oc, oc) - u_nebula_radius * u_nebula_radius;
    float discriminant = b * b - 4.0 * a * c;
    
    vec4 color = vec4(0.0);
    
    if (discriminant > 0.0) {
        float t0 = (-b - sqrt(discriminant)) / (2.0 * a);
        float t1 = (-b + sqrt(discriminant)) / (2.0 * a);
        
        float t_enter = max(0.0, t0);
        float t_exit = min(t1, far_plane);
        
        if (t_exit > t_enter) {
            color = raymarch_nebula(ro + rd * t_enter, rd, t_exit - t_enter);
        }
    }
    
    // Blend with existing scene (assume premultiplied or proper blend mode)
    frag_color = color;
}
#endif

/*============================================================================
 * 6. STARFIELD/PARTICLE VERTEX (Parallax Layers)
 *===========================================================================*/

#ifdef LAGRANGE_SHADER_STARFIELD_VERT
layout(location = 0) in vec3 a_position;      // Vertex
layout(location = 1) in vec4 a_instance_pos;  // xyz = world pos, w = size
layout(location = 2) in vec4 a_instance_vel; // xyz = velocity, w = fade
layout(location = 3) in vec4 a_instance_color;

out vec4 v_color;
out float v_size;
out float v_fade;

uniform float u_parallax_depth; // 0 = distant, 1 = close

void main() {
    // Parallax: shift based on camera movement
    vec3 world_pos = a_instance_pos.xyz;
    
    // Deeper layers move slower (inverse parallax)
    float parallax_scale = 1.0 - u_parallax_depth * 0.9;
    vec3 camera_relative = camera_pos * parallax_scale;
    
    // Add velocity streak
    vec3 streak = a_instance_vel.xyz * 0.1; // Trail length
    vec3 view_pos = world_pos - camera_relative + streak * a_position.x;
    
    // Billboard
    vec2 billboard = a_position.xy * a_instance_pos.w;
    view_pos += (inverse(mat3(view)) * vec3(billboard, 0.0)).xyz;
    
    gl_Position = projection * view * vec4(view_pos, 1.0);
    
    v_color = a_instance_color;
    v_size = a_instance_pos.w;
    v_fade = a_instance_vel.w;
    
    // Size attenuation
    float dist = length(view_pos);
    gl_PointSize = a_instance_pos.w * (1000.0 / dist);
}
#endif

/*============================================================================
 * 7. STARFIELD/PARTICLE FRAGMENT
 *===========================================================================*/

#ifdef LAGRANGE_SHADER_STARFIELD_FRAG
in vec4 v_color;
in float v_size;
in float v_fade;
out vec4 frag_color;

void main() {
    // Circular particle with soft edge
    vec2 uv = gl_PointCoord * 2.0 - 1.0;
    float dist = length(uv);
    
    // Core + glow
    float core = 1.0 - smoothstep(0.0, 0.5, dist);
    float glow = exp(-dist * 3.0) * 0.5;
    
    float alpha = (core + glow) * v_fade;
    
    // Motion blur streak (anisotropic)
    // Would need velocity direction passed through
    
    frag_color = vec4(v_color.rgb, v_color.a * alpha);
}
#endif

/*============================================================================
 * 8. ENGINE TRAIL VERTEX (Ribbon/Line)
 *===========================================================================*/

#ifdef LAGRANGE_SHADER_TRAIL_VERT
layout(location = 0) in vec3 a_position;  // Along trail: x = segment, y = width, z = fade
layout(location = 1) in vec3 a_normal;

out vec2 v_uv;
out float v_fade;
out vec3 v_world_pos;

uniform mat4 u_model;
uniform vec3 u_trail_points[64]; // World-space ribbon points
uniform int u_trail_count;
uniform float u_trail_width;
uniform float u_trail_fade_start;

void main() {
    int segment = int(a_position.x);
    float t = fract(a_position.x);
    
    // Catmull-Rom or linear interpolation between points
    vec3 p0 = u_trail_points[max(0, segment - 1)];
    vec3 p1 = u_trail_points[segment];
    vec3 p2 = u_trail_points[min(u_trail_count - 1, segment + 1)];
    vec3 p3 = u_trail_points[min(u_trail_count - 1, segment + 2)];
    
    // Cubic interpolation
    vec3 pos = 0.5 * (
        2.0 * p1 +
        (p2 - p0) * t +
        (2.0 * p0 - 5.0 * p1 + 4.0 * p2 - p3) * t * t +
        (3.0 * p1 - p0 - 3.0 * p2 + p3) * t * t * t
    );
    
    // Width expansion from velocity/heat
    vec3 tangent = normalize(p2 - p1);
    vec3 bitangent = normalize(cross(tangent, camera_pos - pos));
    
    pos += bitangent * a_position.y * u_trail_width * (1.0 + t * 2.0);
    
    v_world_pos = pos;
    v_uv = vec2(t, a_position.y * 0.5 + 0.5);
    v_fade = 1.0 - (float(segment) + t) / float(u_trail_count);
    v_fade = smoothstep(0.0, u_trail_fade_start, v_fade);
    
    gl_Position = view_proj * vec4(pos, 1.0);
}
#endif

/*============================================================================
 * 9. ENGINE TRAIL FRAGMENT (Heat distortion + Glow)
 *===========================================================================*/

#ifdef LAGRANGE_SHADER_TRAIL_FRAG
in vec2 v_uv;
in float v_fade;
in vec3 v_world_pos;
out vec4 frag_color;

uniform vec3 u_trail_color;
uniform float u_heat_intensity;
uniform float u_time;

void main() {
    // Core heat (white-blue)
    vec3 core_color = vec3(0.9, 0.95, 1.0);
    // Outer glow (engine color)
    vec3 outer_color = u_trail_color;
    
    // Radial gradient
    float center = 1.0 - abs(v_uv.y * 2.0 - 1.0);
    center = pow(center, 0.5);
    
    // Flicker noise
    float flicker = 0.9 + 0.1 * sin(u_time * 20.0 + v_uv.x * 50.0);
    
    // Heat distortion shimmer (would be displacement in vertex shader)
    float shimmer = sin(v_uv.x * 30.0 - u_time * 5.0) * 0.5 + 0.5;
    
    vec3 color = mix(outer_color, core_color, center * u_heat_intensity);
    color *= flicker * (1.0 + shimmer * 0.2);
    
    float alpha = center * v_fade * u_heat_intensity;
    
    // Additive blending expected
    frag_color = vec4(color * alpha, alpha);
}
#endif

/*============================================================================
 * 10. POST-PROCESS VERTEX (Full-screen triangle)
 *===========================================================================*/

#ifdef LAGRANGE_SHADER_POST_VERT
const vec2 positions[3] = vec2[](
    vec2(-1.0, -1.0),
    vec2( 3.0, -1.0),
    vec2(-1.0,  3.0)
);

out vec2 v_uv;

void main() {
    gl_Position = vec4(positions[gl_VertexID], 0.0, 1.0);
    v_uv = positions[gl_VertexID] * 0.5 + 0.5;
}
#endif

/*============================================================================
 * 11. POST-PROCESS FRAGMENT (Bloom, Warp, Grain)
 *===========================================================================*/

#ifdef LAGRANGE_SHADER_POST_FRAG
in vec2 v_uv;
out vec4 frag_color;

uniform sampler2D u_scene_color;      // HDR scene
uniform sampler2D u_scene_depth;      // For depth effects
uniform sampler2D u_bloom_blur;       // Blurred highlights
uniform sampler2D u_lens_dirt;        // Lens dirt overlay

uniform float u_bloom_intensity;
uniform float u_warp_speed;           // 0 = normal, >0 = speed effect
uniform float u_chromatic_aberration;
uniform float u_vignette;
uniform float u_film_grain;

/* Depth-based fog */
vec3 apply_fog(vec3 color, float depth, vec3 fog_color) {
    float fog_factor = 1.0 - exp(-fog_density * depth * fog_height_falloff);
    return mix(color, fog_color, clamp(fog_factor, 0.0, 1.0));
}

/* Speed warp: radial blur + FOV stretch */
vec2 warp_uv(vec2 uv, float strength) {
    vec2 center = uv - 0.5;
    float dist = length(center);
    float angle = atan(center.y, center.x);
    
    // Radial stretch
    float stretch = 1.0 + strength * dist * 2.0;
    center *= stretch;
    
    // Rotation blur
    float blur = strength * dist * 0.1;
    angle += blur * sin(time * 10.0);
    
    return vec2(
        cos(angle) * length(center),
        sin(angle) * length(center)
    ) + 0.5;
}

/* Chromatic aberration */
vec3 chromatic_sample(sampler2D tex, vec2 uv, float amount) {
    vec2 dir = normalize(uv - 0.5);
    float r = texture(tex, uv + dir * amount).r;
    float g = texture(tex, uv).g;
    float b = texture(tex, uv - dir * amount).b;
    return vec3(r, g, b);
}

/* Film grain */
float grain(vec2 uv) {
    return hash(uv + fract(time));
}

void main() {
    vec2 uv = v_uv;
    
    // Speed warp effect
    if (u_warp_speed > 0.01) {
        uv = warp_uv(uv, u_warp_speed);
    }
    
    // Scene sampling with chromatic aberration
    vec3 scene = chromatic_sample(u_scene_color, uv, u_chromatic_aberration * u_warp_speed);
    
    // Depth for fog
    float depth = texture(u_scene_depth, uv).r;
    depth = 2.0 * near_plane * far_plane / (far_plane + near_plane - (2.0 * depth - 1.0) * (far_plane - near_plane));
    
    scene = apply_fog(scene, depth, ambient_color * 0.5);
    
    // Bloom addition
    vec3 bloom = texture(u_bloom_blur, uv).rgb * u_bloom_intensity;
    scene += bloom;
    
    // Lens dirt (vignette + dust)
    float vignette = 1.0 - dot(uv - 0.5, uv - 0.5) * u_vignette;
    vec3 dirt = texture(u_lens_dirt, uv).rgb;
    scene *= vignette * (1.0 - dirt * 0.3);
    
    // Film grain
    scene *= (1.0 - u_film_grain * 0.5) + grain(uv * resolution) * u_film_grain;
    
    // Final tone mapping
    scene = tonemap_aces(scene * exposure);
    scene = pow(scene, vec3(1.0 / 2.2));
    
    frag_color = vec4(scene, 1.0);
}
#endif

/*============================================================================
 * 12. ORBIT PATH VISUALIZATION (Line rendering)
 *===========================================================================*/

#ifdef LAGRANGE_SHADER_ORBIT_VERT
layout(location = 0) in vec3 a_position; // Point on orbit
layout(location = 1) in float a_param;   // 0 = start, 1 = end for fade

out float v_param;
out vec3 v_world_pos;

uniform mat4 u_model; // Usually identity, or planet-local

void main() {
    v_world_pos = (u_model * vec4(a_position, 1.0)).xyz;
    v_param = a_param;
    gl_Position = view_proj * vec4(v_world_pos, 1.0);
}
#endif

#ifdef LAGRANGE_SHADER_ORBIT_FRAG
in float v_param;
in vec3 v_world_pos;
out vec4 frag_color;

uniform vec3 u_orbit_color;
uniform float u_orbit_brightness;
uniform float u_fade_distance;

void main() {
    // Distance fade
    float dist = length(v_world_pos - camera_pos);
    float fade = exp(-dist / u_fade_distance);
    
    // Parametric fade at ends
    float end_fade = smoothstep(0.0, 0.1, v_param) * smoothstep(1.0, 0.9, v_param);
    
    // Pulsing highlight
    float pulse = 0.8 + 0.2 * sin(time * 2.0);
    
    vec3 color = u_orbit_color * u_orbit_brightness * pulse;
    float alpha = fade * end_fade;
    
    frag_color = vec4(color, alpha);
}
#endif

/*============================================================================
 * 13. SOI SPHERE VISUALIZATION
 *===========================================================================*/

#ifdef LAGRANGE_SHADER_SOI_VERT
layout(location = 0) in vec3 a_position;

out vec3 v_normal;
out vec3 v_world_pos;

uniform mat4 u_model;
uniform float u_radius;

void main() {
    v_normal = a_position;
    v_world_pos = (u_model * vec4(a_position * u_radius, 1.0)).xyz;
    gl_Position = view_proj * vec4(v_world_pos, 1.0);
}
#endif

#ifdef LAGRANGE_SHADER_SOI_FRAG
in vec3 v_normal;
in vec3 v_world_pos;
out vec4 frag_color;

uniform vec3 u_color;
uniform float u_opacity;
uniform float u_pulse;

void main() {
    // Fresnel effect
    vec3 view_dir = normalize(camera_pos - v_world_pos);
    float fresnel = pow(1.0 - abs(dot(view_dir, normalize(v_normal))), 2.0);
    
    // Grid pattern on sphere
    vec3 n = normalize(v_normal);
    float grid = step(0.98, fract(n.x * 10.0)) + step(0.98, fract(n.y * 10.0)) + step(0.98, fract(n.z * 10.0));
    
    // Pulsing boundary
    float pulse = 0.5 + 0.5 * sin(time * u_pulse);
    
    vec3 color = mix(u_color, vec3(1.0), grid * 0.5);
    float alpha = fresnel * u_opacity * (0.5 + 0.5 * pulse);
    
    frag_color = vec4(color, alpha);
}
#endif

/*============================================================================
 * 14. C INTERFACE (for C99 integration)
 *===========================================================================*/

#ifdef __cplusplus
extern "C" {
#endif

/* Shader program handles */
typedef struct {
    uint32_t skybox_prog;
    uint32_t nebula_prog;
    uint32_t starfield_prog;
    uint32_t trail_prog;
    uint32_t post_prog;
    uint32_t orbit_prog;
    uint32_t soi_prog;
    
    /* Uniform locations cached */
    struct {
        int exposure;
        int time;
        int camera_pos;
        int warp_speed;
    } post_locs;
} lg_shader_pack_t;

/* Compile all shaders */
lg_shader_pack_t* lg_shader_pack_compile(void);
void lg_shader_pack_destroy(lg_shader_pack_t* pack);

/* Update per-frame uniforms */
void lg_shader_pack_update_camera(lg_shader_pack_t* pack, const lg_camera_t* cam);
void lg_shader_pack_update_time(lg_shader_pack_t* pack, float time, float exposure);

/* Effect triggers */
void lg_shader_pack_set_warp(lg_shader_pack_t* pack, float speed_factor); // 0-1
void lg_shader_pack_set_bloom(lg_shader_pack_t* pack, float intensity);

#ifdef __cplusplus
}
#endif

#endif /* LAGRANGE_SHADERS_GLSL */

