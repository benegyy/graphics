#include <iostream>
#include "parser.h"
#include "ppm.h"
#include <cmath>
#include <math.h>
#include <limits>
//Ray-Tracer Hw
//Authors : Alpay Avşar - 2171254
//          Bengisu Karaca - 2448538

typedef unsigned char RGB[3];

using namespace parser;

enum INFO
{
    MISS = 0,
    HIT
};
struct Ray
{
    Vec3f origin, direction;
};
const RGB BAR_COLOR[8] =
    {
        {255, 255, 255}, // 100% White
        {255, 255, 0},   // Yellow
        {0, 255, 255},   // Cyan
        {0, 255, 0},     // Green
        {255, 0, 255},   // Magenta
        {255, 0, 0},     // Red
        {0, 0, 255},     // Blue
        {0, 0, 0},       // Black
};

Scene myscene;
Vec3f colorize(Ray ray, const Camera &cam, int recursion_depth);

struct HitInfo
{
    INFO info;
    float t;
    Vec3f position;
    Vec3f normal;
    int material_id;
    int touched_object_id;
};
float dotProduct(Vec3f vec1, Vec3f vec2)
{
    return vec1.x * vec2.x + vec1.y * vec2.y + vec1.z * vec2.z;
}
Vec3f scalar_product(float s, Vec3f vec)
{
    Vec3f c{};
    c.x = vec.x * s;
    c.y = vec.y * s;
    c.z = vec.z * s;
    return c;
}
Vec3f negate(Vec3f vec)
{
    return scalar_product(-1, vec);
}
Vec3f add(Vec3f a, Vec3f b)
{
    Vec3f result;
    result.x = a.x + b.x;
    result.y = a.y + b.y;
    result.z = a.z + b.z;
    return result;
}
Vec3f get_otd(Ray r, float t)
{
    return add(r.origin, scalar_product(t, r.direction));
}
float sqr_vec(Vec3f vec)
{
    return dotProduct(vec, vec);
}
Vec3f vector_scale(const Vec3f &lhs, const Vec3f &rhs)
{
    Vec3f result;
    result.x = lhs.x * rhs.x;
    result.y = lhs.y * rhs.y;
    result.z = lhs.z * rhs.z;
    return result;
}
Vec3f multScaler(Vec3f a, float s)
{
    Vec3f result;
    result.x = a.x * s;
    result.y = a.y * s;
    result.z = a.z * s;
    return result;
}

Vec3f subtract(Vec3f a, Vec3f b)
{
    Vec3f result;
    result.x = a.x - b.x;
    result.y = a.y - b.y;
    result.z = a.z - b.z;
    return result;
}
Vec3f cross_product(Vec3f vec1, Vec3f vec2)
{
    Vec3f result;
    result.x = (vec1.y * vec2.z) - (vec1.z * vec2.y);
    result.y = (vec1.z * vec2.x) - (vec1.x * vec2.z);
    result.z = (vec1.x * vec2.y) - (vec1.y * vec2.x);
    return result;
}

Vec3f normalize(Vec3f vec)
{
    return multScaler(vec, 1.0 / sqrt(dotProduct(vec, vec)));
}

float vectorLength(Vec3f vec)
{
    return sqrtf(powf(vec.x, 2) + powf(vec.y, 2) + powf(vec.z, 2));
}
float pointDistance(Vec3f vec1, Vec3f vec2)
{
    return vectorLength(add(vec1, multScaler(vec2, -1)));
}
float Determinant(const Vec3f &a, const Vec3f &b, const Vec3f &c)
{
    float result = a.x * (b.y * c.z - c.y * b.z) + a.y * (c.x * b.z - b.x * c.z) + a.z * (b.x * c.y - b.y * c.x);
    return result;
}

float matrix_determinant(float m[3][3])
{
    float result;
    result = (m[0][0] * (m[1][1] * m[2][2] - m[2][1] * m[1][2])) - (m[1][0] * (m[2][2] * m[0][1] - m[0][2] * m[2][1])) + (m[2][0] * (m[0][1] * m[1][2] - m[0][2] * m[1][1]));
    return result;
}

unsigned char *InitializeImage(int width, int height)
{
    auto *image = new unsigned char[width * height * 3];
    return image;
}

Ray generateRay(int i, int j, const Camera &cam)
{
    Ray result{};
    float su, sv;
    Vec3f m, q, s, u, v, w, e;
    float nx = cam.image_width;
    float ny = cam.image_height;
    float top, bottom, left, right;
    left = cam.near_plane.x;
    right = cam.near_plane.y;
    top = cam.near_plane.w;
    bottom = cam.near_plane.z;

    su = (i + 0.5) * (right - left) / nx;
    sv = (j + 0.5) * (top - bottom) / ny;

    e = cam.position;
    v = cam.up;
    w = multScaler(cam.gaze, -1);
    u = cross_product(v, w);

    m = add(e, (multScaler(cam.gaze, cam.near_distance)));
    q = add(m, add(multScaler(u, left), multScaler(v, top)));
    s = add(q, add(multScaler(u, su), multScaler(v, -sv)));
    result.origin = e;
    result.direction = add(s, multScaler(e, -1));
    return result;
}

HitInfo sphereIntersection(Ray ray, const Vec3f &center, const Sphere &sphere)

{
    HitInfo hitInfo = {MISS, -1, {0, 0, 0}, {0, 0, 0}, -1, -1};

    float A, B, C;

    float delta;

    float t;

    C = (ray.origin.x - center.x) * (ray.origin.x - center.x) + (ray.origin.y - center.y) * (ray.origin.y - center.y) +
        (ray.origin.z - center.z) * (ray.origin.z - center.z) - (sphere.radius * sphere.radius);

    B = 2 * ray.direction.x * (ray.origin.x - center.x) + 2 * ray.direction.y * (ray.origin.y - center.y) + 2 * ray.direction.z * (ray.origin.z - center.z);

    A = ray.direction.x * ray.direction.x + ray.direction.y * ray.direction.y + ray.direction.z * ray.direction.z;

    delta = B * B - 4 * A * C;

    if (delta < 0)
        return hitInfo; // no intersection
    else
    {
        delta = sqrt(delta);
        t = (-B - delta) / (2 * A); // closest intersection t1
        hitInfo.info = HIT;
        hitInfo.t = t;
        hitInfo.position = add(ray.origin, multScaler(ray.direction, t));
        Vec3f surface_normal = normalize(subtract(hitInfo.position, center));
        hitInfo.normal = surface_normal;
        hitInfo.material_id = sphere.material_id;
        return hitInfo;
    }
}
HitInfo triangleIntersection(Ray ray, Scene &scene, parser::Triangle triangle, const int material_id)
{
    HitInfo hitInfo = {MISS, -1, {0, 0, 0}, {0, 0, 0}, -1, -1};
    Vec3f face1 = scene.vertex_data[triangle.indices.v0_id - 1];
    Vec3f face2 = scene.vertex_data[triangle.indices.v1_id - 1];
    Vec3f face3 = scene.vertex_data[triangle.indices.v2_id - 1];

    float a_matrix[3][3] = {
        {face1.x - face2.x, face1.x - face3.x, ray.direction.x},
        {face1.y - face2.y, face1.y - face3.y, ray.direction.y},
        {face1.z - face2.z, face1.z - face3.z, ray.direction.z}};

    float beta_matrix[3][3] = {
        {face1.x - ray.origin.x, face1.x - face3.x, ray.direction.x},
        {face1.y - ray.origin.y, face1.y - face3.y, ray.direction.y},
        {face1.z - ray.origin.z, face1.z - face3.z, ray.direction.z}

    };

    float gama_matrix[3][3] = {
        {face1.x - face2.x, face1.x - ray.origin.x, ray.direction.x},
        {face1.y - face2.y, face1.y - ray.origin.y, ray.direction.y},
        {face1.z - face2.z, face1.z - ray.origin.z, ray.direction.z}

    };
    float t_matrix[3][3] = {
        {face1.x - face2.x, face1.x - face3.x, face1.x - ray.origin.x},
        {face1.y - face2.y, face1.y - face3.y, face1.y - ray.origin.y},
        {face1.z - face2.z, face1.z - face3.z, face1.z - ray.origin.z}

    };

    float a_det = matrix_determinant(a_matrix);
    if (a_det == 0)
        return hitInfo;
    float gama_det = matrix_determinant(gama_matrix);
    float beta_det = matrix_determinant(beta_matrix);
    float t_det = matrix_determinant(t_matrix);
    float beta = beta_det / a_det;
    float gama = gama_det / a_det;
    float t = t_det / a_det;

    if (t > 0 && beta > 0 && gama > 0 && beta + gama < 1)
    {
        hitInfo.info = HIT;
        hitInfo.t = t;
        hitInfo.position = add(ray.origin, multScaler(ray.direction, t));
        hitInfo.normal = normalize(cross_product(subtract(face2, face1), subtract(face3, face1)));
        hitInfo.material_id = material_id;
        return hitInfo;
    }
    else
    {
        return hitInfo;
    }
}

HitInfo triangleIntersection2(Ray ray, const Vec3f &a, const Vec3f &b, const Vec3f &c, const int material_id)
{
    HitInfo hitInfo = {MISS, -1, {0, 0, 0}, {0, 0, 0}, -1, -1};

    Vec3f a_b = subtract(a, b);
    Vec3f a_c = subtract(a, c);
    Vec3f a_o = subtract(a, ray.origin);

    float det_A = Determinant(a_b, a_c, ray.direction);
    if (det_A == 0.0)
    {
        return hitInfo;
    }
    float t = (Determinant(a_b, a_c, a_o)) / det_A;
    if (t <= 0.0)
    {
        return hitInfo;
    }
    float gamma = Determinant(a_b, a_o, ray.direction) / det_A;
    if (gamma < 0 || gamma > 1)
    {
        return hitInfo;
    }
    float beta = Determinant(a_o, a_c, ray.direction) / det_A;
    if (beta < 0 || beta > 1 || beta + gamma > 1)
    {
        return hitInfo;
    }
    hitInfo.material_id = material_id;
    hitInfo.t = t;
    hitInfo.position = add(ray.origin, multScaler(ray.direction, t));
    hitInfo.normal = normalize(cross_product(subtract(b, a), subtract(c, a)));
    hitInfo.info = HIT;
    return hitInfo;
}

HitInfo mesh_intersect(Ray ray, Scene &scene, parser::Mesh mesh)
{
    HitInfo hitInfo = {MISS, -1, {0, 0, 0}, {0, 0, 0}, -1, -1};
    hitInfo.t = std::numeric_limits<float>::max();
    int face_count = mesh.faces.size();
    for (int face_no = 0; face_no < face_count; face_no++)
    {
        HitInfo current_attempt;
        Vec3f v0 = scene.vertex_data[mesh.faces[face_no].v0_id - 1];
        Vec3f v1 = scene.vertex_data[mesh.faces[face_no].v1_id - 1];
        Vec3f v2 = scene.vertex_data[mesh.faces[face_no].v2_id - 1];

        ///
        current_attempt = triangleIntersection2(ray, v0, v1, v2, mesh.material_id);
        ///

        if (current_attempt.info == HIT && current_attempt.t < hitInfo.t)
        {
            hitInfo = current_attempt;
            hitInfo.info = HIT;
        }
        else
        {
            continue;
        }
    }
    return hitInfo;
}

HitInfo FindClosestIntersection(const Ray &ray)
{
    HitInfo hitInfo = {MISS, -1, {0, 0, 0}, {0, 0, 0}, -1, -1};
    float min_t = std::numeric_limits<float>::max();

    // Check spheres.
    for (const auto &sphere : myscene.spheres)
    {
        Vec3f sphere_center = myscene.vertex_data.at(sphere.center_vertex_id - 1);
        HitInfo current = sphereIntersection(ray, sphere_center, sphere);
        float curr_t = current.t;
        if (curr_t > 0 && curr_t < min_t)
        {
            min_t = curr_t;
            hitInfo = current;
        }
    }

    // Check triangles.
    for (const auto &triangle : myscene.triangles)
    {
        HitInfo current = triangleIntersection(ray, myscene, triangle, triangle.material_id);
        float curr_t = current.t;
        if (curr_t > 0 && curr_t < min_t)
        {
            min_t = curr_t;
            hitInfo = current;
        }
    }
    // Check meshes.
    for (const auto &mesh : myscene.meshes)
    {

        HitInfo current = mesh_intersect(ray, myscene, mesh);
        float current_t = current.t;
        float curr_t = current.t;
        if (curr_t > 0 && curr_t < min_t)
        {
            min_t = curr_t;
            hitInfo = current;
        }
    }
    return hitInfo;
}

Vec3f ambient_shading(const Vec3f &coefficient, const Vec3f &radiance)
{
    return vector_scale(coefficient, radiance);
}

bool isNotShadow(const Ray &ray, const Vec3f &position, const Scene &scene)
{
    float t = vectorLength(subtract(position, ray.origin)); // Assumed ray.d is unit vector

    for (std::size_t sid = 0; sid < scene.spheres.size(); sid++)
    {
        Vec3f c = scene.vertex_data[scene.spheres[sid].center_vertex_id - 1];
        float r = scene.spheres[sid].radius;
        float t_intersect = sphereIntersection(ray, c, scene.spheres[sid]).t;

        if ((t_intersect >= 0) && (t_intersect < t))
            return false;
    }

    for (std::size_t tid = 0; tid < scene.triangles.size(); tid++)
    {
        Face indice = scene.triangles[tid].indices;
        int material_id = scene.triangles[tid].material_id;
        Vec3f a = scene.vertex_data[indice.v0_id - 1];
        Vec3f b = scene.vertex_data[indice.v1_id - 1];
        Vec3f c = scene.vertex_data[indice.v2_id - 1];

        float t_intersect = triangleIntersection2(ray, a, b, c, material_id).t;

        if ((t_intersect >= 0) && (t_intersect < t))
            return false;
    }

    for (std::size_t meid = 0; meid < scene.meshes.size(); meid++)
    {
        int material_id = scene.meshes[meid].material_id;
        for (std::size_t fid = 0; fid < scene.meshes[meid].faces.size(); fid++)
        {
            Face face = scene.meshes[meid].faces[fid];

            Vec3f a = scene.vertex_data[face.v0_id - 1];
            Vec3f b = scene.vertex_data[face.v1_id - 1];
            Vec3f c = scene.vertex_data[face.v2_id - 1];

            float t_intersect = triangleIntersection2(ray, a, b, c, material_id).t;

            if ((t_intersect >= 0) && (t_intersect < t))
                return false;
        }
    }
    return true;
}

Vec3f diffuse_shading(HitInfo hitinfo, const Vec3f &coefficient, const Scene &scene)
{
    Vec3f result{0, 0, 0};

    for (const PointLight &pointLight : scene.point_lights)
    {

        Vec3f wi = subtract(pointLight.position, hitinfo.position);
        float r = vectorLength(wi);
        wi.x = wi.x / r;
        wi.y = wi.y / r;
        wi.z = wi.z / r;

        Ray shadow_ray;
        shadow_ray.origin = add(hitinfo.position, scalar_product(scene.shadow_ray_epsilon, wi));
        shadow_ray.direction = wi;

        // Check if the point is in shadow
        if (isNotShadow(shadow_ray, pointLight.position, scene))
        {
            float cos_theta = std::max((float)0, dotProduct(wi, hitinfo.normal));

            result = add(result, scalar_product((cos_theta / (r * r)), vector_scale(coefficient, pointLight.intensity)));
        }
    }

    return result;
}

Vec3f specular_shading(const Ray &ray, HitInfo hitInfo, const Vec3f &coefficient, const float phong, const Scene &scene)
{
    Vec3f result{0, 0, 0};
    Vec3f wo = normalize(negate(ray.direction));

    for (const PointLight &pointLight : scene.point_lights)
    {
        // Check if the point is in pointLight

        Vec3f light_vector = subtract(pointLight.position, hitInfo.position);
        float light_distance = sqrt(dotProduct(light_vector, light_vector));
        Vec3f light_contribution;
        Vec3f tmp = subtract(hitInfo.position, pointLight.position);
        light_contribution.x = pointLight.intensity.x / std::pow(std::sqrt(dotProduct(tmp, tmp)), 2);
        light_contribution.y = pointLight.intensity.y / std::pow(std::sqrt(dotProduct(tmp, tmp)), 2);
        light_contribution.z = pointLight.intensity.z / std::pow(std::sqrt(dotProduct(tmp, tmp)), 2);
        Vec3f wi = subtract(pointLight.position, hitInfo.position);
        Ray shadow_ray;
        shadow_ray.origin = add(hitInfo.position, scalar_product(scene.shadow_ray_epsilon, wi));
        shadow_ray.direction = wi;

        // Check if the point is in shadow
        if (isNotShadow(shadow_ray, pointLight.position, scene))
        {
            Vec3f h = normalize(add(normalize(scalar_product(-1, ray.direction)), normalize(light_vector)));

            // Calculate the specular contribution using Blinn-Phong model
            float cos_alpha = dotProduct(hitInfo.normal, h);
            if (cos_alpha > 0)
            {
                float specular_term = std::pow(cos_alpha, phong);
                result = result + scalar_product(specular_term, vector_scale(coefficient, light_contribution));
            }
        }
    }

    return result;
}

Ray reflection_ray(const Vec3f &o, const Vec3f &d, const Scene &scene)
{
    Ray rr;
    rr.origin = o;
    rr.direction = d;
    HitInfo HitInfo;
    HitInfo.t = INFINITY;

    for (auto &current_sphere : scene.spheres)
    {
        Vec3f c = scene.vertex_data[current_sphere.center_vertex_id - 1];
        float r = current_sphere.radius;
        float t_intersect = sphereIntersection(rr, c, current_sphere).t;

        if ((t_intersect >= 0) && (t_intersect < HitInfo.t))
        {
            HitInfo.t = t_intersect;
            HitInfo.material_id = current_sphere.material_id;
            HitInfo.normal = subtract(add(rr.origin, scalar_product(t_intersect, rr.direction)), c);
            HitInfo.normal = normalize(HitInfo.normal);
        }
    }

    for (auto &current_triangle : scene.triangles)
    {
        Face indice = current_triangle.indices;
        Vec3f a = scene.vertex_data[indice.v0_id - 1];
        Vec3f b = scene.vertex_data[indice.v1_id - 1];
        Vec3f c = scene.vertex_data[indice.v2_id - 1];

        float t_intersect = triangleIntersection2(rr, a, b, c, current_triangle.material_id).t;

        if ((t_intersect >= 0) && (t_intersect < HitInfo.t))
        {
            HitInfo.t = t_intersect;
            HitInfo.material_id = current_triangle.material_id;
            HitInfo.normal = cross_product(b - a, c - b);
            HitInfo.normal = normalize(HitInfo.normal);
        }
    }

    for (const auto &mesh : myscene.meshes)
    {
        for (std::size_t fid = 0; fid < mesh.faces.size(); fid++)
        {
            Face face = mesh.faces[fid];

            Vec3f a = scene.vertex_data[face.v0_id - 1];
            Vec3f b = scene.vertex_data[face.v1_id - 1];
            Vec3f c = scene.vertex_data[face.v2_id - 1];

            float t_intersect = triangleIntersection2(rr, a, b, c, mesh.material_id).t;

            if ((t_intersect >= 0) && (t_intersect < HitInfo.t))
            {
                HitInfo.t = t_intersect;
                HitInfo.material_id = mesh.material_id;
                HitInfo.normal = cross_product(b - a, c - b);
                HitInfo.normal = normalize(HitInfo.normal);
            }
        }
    }

    return rr;
}

Vec3f specular_reflection(const Ray &ray, HitInfo hitinfo, const Vec3f &coefficient, const Scene &scene, const Camera &cam, int recursion_depth)
{
    if (!coefficient.x && !coefficient.y && !coefficient.z)
        return Vec3f{0, 0, 0};
    else if (recursion_depth > scene.max_recursion_depth) 
        return Vec3f{0, 0, 0};
    else
    {
        Vec3f wo = normalize(-ray.direction);
        Vec3f wr = -wo + 2 * hitinfo.normal * dotProduct(hitinfo.normal, wo);

        Ray rr = reflection_ray(hitinfo.position + wr * scene.shadow_ray_epsilon, wr, scene);
        return vector_scale(coefficient, colorize(rr, cam, recursion_depth));
    }
}

Vec3f colorize(Ray ray, const Camera &cam, int recursion_depth)
{

    HitInfo HitInfo = FindClosestIntersection(ray);
    Vec3f color = {0, 0, 0};

    if (HitInfo.info == HIT)
    {
        Material touched_mat = myscene.materials[HitInfo.material_id - 1];

        // implement reflection
        // ambient
        // diffuse
        // specular
        
        Vec3f mirror_reflection = specular_reflection(ray, HitInfo, touched_mat.mirror, myscene, cam, recursion_depth + 1);
        Vec3f ambient_component = ambient_shading(touched_mat.ambient, myscene.ambient_light);
        Vec3f diffuse_component = diffuse_shading(HitInfo, touched_mat.diffuse, myscene);
        Vec3f specular_component = specular_shading(ray, HitInfo, touched_mat.specular, touched_mat.phong_exponent, myscene);

        color = color + ambient_component + diffuse_component + specular_component + mirror_reflection;
        color.x = color.x > 255 ? 255 : color.x;
        color.y = color.y > 255 ? 255 : color.y;
        color.z = color.z > 255 ? 255 : color.z;
        return color;
    }
    else if (HitInfo.info == MISS || HitInfo.t < 0)
    {
        color.x = (float)myscene.background_color.x;
        color.y = (float)myscene.background_color.y;
        color.z = (float)myscene.background_color.z;
        return color;
    }
    return color;

}

unsigned char *Render(const Camera &cam, const int width, const int height)
{
    unsigned char *image = InitializeImage(width, height);
    unsigned int index = 0;

    for (int row = 0; row < cam.image_height; row++)
    {
        for (int column = 0; column < cam.image_width; column++)
        {
            Vec3f pixel_color{};
            Ray ray = generateRay(column, row, cam);
            pixel_color = colorize(ray, cam, 0);
            RGB clamped_pixel_color;
            clamped_pixel_color[0] = pixel_color.x > 255 ? 255 : (int)pixel_color.x + 0.5;
            clamped_pixel_color[1] = pixel_color.y > 255 ? 255 : (int)pixel_color.y + 0.5;
            clamped_pixel_color[2] = pixel_color.z > 255 ? 255 : (int)pixel_color.z + 0.5;
            image[index++] = clamped_pixel_color[0];
            image[index++] = clamped_pixel_color[1];
            image[index++] = clamped_pixel_color[2];
        }
    }
    return image;
}

int main(int argc, char *argv[])
{
    // Sample usage for reading an XML scene file
    myscene.loadFromXml(argv[1]);
    // The code below creates a test pattern and writes
    // it to a PPM file to demonstrate the usage of the
    // ppm_write function.
    //
    // Normally, you would be running your ray tracing
    // code here to produce the desired image.

    for (const auto &camera : myscene.cameras)
    {
        auto width = camera.image_width;
        auto height = camera.image_height;
        unsigned char *image = Render(camera, width, height);
        write_ppm(camera.image_name.c_str(), image, width, height);
    }
}
