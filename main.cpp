#include <iostream>
#include <cstdint>
#include <string>
#include <vector>
#include <unordered_map>
#include <limits>
#include "Bitmap.h"
#include "Geometry.h"

constexpr unsigned width = 1024;
constexpr unsigned height = 1024;
constexpr float fov = M_PI/2.0f;
constexpr size_t RECURSION_DEPTH = 2;
constexpr float EPSILON = 0.0001;
constexpr float MAX_FLOAT = std::numeric_limits<float>::max();
constexpr unsigned HORIZON = 1000;
const Vec3f BACKGROUND_COLOR(0.2, 0.7, 0.8);

//возвращает пересечение со сценой(сферы, текстура, примитив)
//изменяет hit - т.пересечения, N - нормаль, material
bool scene_intersect(
        const Vec3f &orig,
        const Vec3f &dir,
        const std::vector<Sphere> &spheres,
        const std::vector<Triangle> &triangles,
        const Plane &plane,
        Vec3f &hit, Vec3f &N,
        Material &material) {

    float spheres_dist = MAX_FLOAT;
    for (auto const& sphere: spheres) {
        float dist_i;
        if (sphere.intersect(orig, dir, N, hit, dist_i) && dist_i < spheres_dist) {
            spheres_dist = dist_i;
            material = sphere.material;
        }
    }

    for (auto const& tr: triangles) {
        float dist_i;
        if (tr.intersect(orig, dir, N, hit, dist_i) && dist_i < spheres_dist) {
            spheres_dist = dist_i;
            material = tr.material;
        }
    }

    float plane_dist = MAX_FLOAT;
    Vec3f nn;
    Vec3f point;
    float t = MAX_FLOAT;
    if (plane.intersect(orig, dir, nn, point, t) && fabs(point.x) < 30 && point.z < -10
        &&  t < spheres_dist) {
        N = nn;
        hit = point;
        plane_dist = t;
        //шахматная раскраска
        if (!(int(0.4 * hit.x + 100) + int(0.4 * hit.z + 100) & 1)) {
            material = plane.m2;
        } else {
            material = plane.m1;
        }
    }
    return std::min(spheres_dist, plane_dist) < HORIZON;
}

//возвращает материал/цвет
Vec3f cast_ray(
        const Vec3f &orig,
        const Vec3f &dir,
        const std::vector<Sphere> &spheres,
        const std::vector<Triangle> &triangles,
        const std::vector<Light> &lights,
        const Plane &plane,
        size_t depth = 0) {

    Vec3f point, N;
    Material material;


    float diffuse_light_intensity = 0;
    float specular_light_intensity = 0;

    if (depth > RECURSION_DEPTH ||
        !scene_intersect(orig, dir, spheres, triangles, plane, point, N, material)) {//здесь point, N, material, изменяются
        return BACKGROUND_COLOR;
    }

    Vec3f refr_dir = refract(dir, N, material.refractive_index).normalize();//вектор преломления
    Vec3f refr_orig = refr_dir*N < 0 ? point - N*EPSILON : point + N*EPSILON;//точка пересечения с объектом для преломления со смещением
    Vec3f refr_color = cast_ray(refr_orig, refr_dir, spheres, triangles, lights, plane, depth + 1);//пускаем преломленный луч

    Vec3f refl_dir = reflect(dir, N).normalize();//вектор отражения
    Vec3f refl_orig = refl_dir*N < 0 ? point - N*EPSILON : point + N*EPSILON; //точка пересечения с объектом для отражения со смещением
    Vec3f refl_color = cast_ray(refl_orig, refl_dir, spheres, triangles, lights, plane, depth + 1);//пускаем отраженный луч

    for (auto const &light: lights) {
        Vec3f light_dir      = (light.position - point);//вектор света к точке пересечения

        float light_dist = light_dir.norm();//расстояние до точки
        light_dir = light_dir.normalize();

        Vec3f shadow_orig = light_dir*N < 0 ? point - N*EPSILON : point + N*EPSILON; //проверка нахождения точки в тени
        Vec3f shadow_pt, shadow_N;

        Material tmpmat;
        if (scene_intersect(shadow_orig, light_dir, spheres, triangles, plane,shadow_pt, shadow_N, tmpmat)
            && (shadow_pt - shadow_orig).norm() < light_dist)//если луч пересекает объекты сцены на своем пути, то пропускаем
            continue;
        //Модель Фонга
        specular_light_intensity += powf(std::max(0.f, reflect(light_dir, N)*dir),
                material.specular_exponent)*light.intensity;
        diffuse_light_intensity  += light.intensity * std::max(0.f, light_dir*N);
    }

    return material.diffuse_color * diffuse_light_intensity * material.albedo[0] +
        Vec3f(1., 1., 1.)*specular_light_intensity * material.albedo[1] +
        refl_color*material.albedo[2] + refr_color*material.albedo[3];
}

int toRGB(float R, float G, float B) {
    unsigned r = 255*(R < 0 ? 0 : R);
    unsigned g = 255*(G < 0 ? 0 : G);
    unsigned b = 255*(B < 0 ? 0 : B);
    return (b<<16) | (g<<8) | r;
}

int main(int argc, const char** argv)
{
    std::unordered_map<std::string, std::string> cmdLineParams;
    for(int i = 0; i < argc; ++i) {
        std::string key(argv[i]);
        if(key.size() > 0 && key[0]=='-') {
            if(i != argc-1) // not last argument {
                cmdLineParams[key] = argv[i+1];
                ++i;
            }
        else {
            cmdLineParams[key] = "";
        }
    }

    std::string outFilePath = "out.bmp";
    if(cmdLineParams.find("-out") != cmdLineParams.end())
    outFilePath = cmdLineParams["-out"];

    int sceneId = 1;
    if(cmdLineParams.find("-scene") != cmdLineParams.end())
        sceneId = atoi(cmdLineParams["-scene"].c_str());
    if (sceneId != 1) {
        return 0;
    }

    size_t num_thr = 0;
    if(cmdLineParams.find("-threads") != cmdLineParams.end())
        num_thr = atoi(cmdLineParams["-threads"].c_str());
//-------------------------------------------------------------/

    std::vector<Vec3f> image(width*height);

    const Material glass(1.5, Vec4f(0.0,  0.5, 0.1, 0.8), Vec3f(0.6, 0.7, 0.8),  125.);
    const Material yellow(1.0, Vec4f(0.6,  0.3, 0.1, 0.0), Vec3f(1.0, 0.77, 0.45),   50.);
    const Material mirror(1.0, Vec4f(0.0, 10.0, 0.8, 0.0), Vec3f(1.0, 1.0, 1.0), 1425.);
    const Material for_pl1(1.0, Vec4f(0.6,  0.3, 0.1, 0.0), Vec3f(0.5,0.5, 0.5),   50.);
    const Material for_pl2(1.0, Vec4f(0.6,  0.3, 0.1, 0.0), Vec3f(0.069, 0.234, 0.213),   50.);
    const Material pink(1.0, Vec4f(0.6,  0.3, 0.1, 0.0), Vec3f(1.0, 0.478, 0.45),   50.);


    Plane plane(Vec3f(0,1,0), Vec3f(0, -6, 0), for_pl1, for_pl2);

    std::vector<Sphere> spheres;
    spheres.emplace_back(Sphere(Vec3f( 6, -2, -21), 2.5, pink));
    spheres.emplace_back(Sphere(Vec3f( 5,    5,   -20), 3,      mirror));
    spheres.emplace_back(Sphere(Vec3f(-8, -2, -24), 2.5,      glass));
    spheres.emplace_back(Sphere(Vec3f(-2,    -4,   -14), 2,      yellow));


    std::vector<Triangle> triangles;
    Triangle tr1(Vec3f(-3, -2, -26),
                 Vec3f( 3, -2, -26),
                 Vec3f( 0,  3, -26),
                 glass);
    triangles.emplace_back(tr1);

    std::vector<Light>  lights;
    lights.emplace_back(Light(Vec3f(-30, 25,  25), 1.4));
    lights.emplace_back(Vec3f( 30, 50, -35), 1.4);
    lights.emplace_back(Light(Vec3f( 30, 25,  35), 1.4));

    //для устранения ступенчености
    std::vector <std::pair<float, float>> bias = {
            {-0.125,  0.375},
            {0.375,  0.125},
            {0.375, 0.125},
            {0.125, -0.375},
    };
#pragma omp parallel for num_threads(num_thr)
    for (size_t j = 0; j < height; ++j) {
        for (size_t i = 0; i < width; ++i) {
            for (auto _bias: bias) {
                float x = i - float(width)/2 + _bias.first;
                float y = j  - float(height)/2 + _bias.second;
                float z = - width/ (tan(fov));

                Vec3f dir = Vec3f(x, y, 4*z).normalize();
                image[i+j*width] = image[i+j*width] + cast_ray(Vec3f(0,0,0), dir, spheres, triangles, lights, plane);
            }
            image[i+j*width] = image[i+j*width] * 0.25;
        }
    }
    std::vector<uint32_t> buf(width*height);
#pragma omp parallel for num_threads(num_thr)
    for (size_t i = 0; i < image.size(); ++i) {
        Vec3f &c = image[i];
        float max = std::max(c[0], std::max(c[1], c[2]));
        if (max > 1) {
            buf[i] = toRGB(image[i][0] / max, image[i][1] / max, image[i][2] / max);
        } else {
            buf[i] = toRGB(image[i][0], image[i][1], image[i][2]);
        }
    }

    SaveBMP(outFilePath.c_str(), buf.data(), width, height);

    std::cout << "end." << std::endl;
    return 0;
}