//
// Created by alex on 15.03.2020.
//

#ifndef RT_TRACER_H
#define RT_TRACER_H

#include <cmath>
#include <iostream>
#include <vector>
#include <cassert>
#include <array>

template <size_t N, typename T>
class vec {
    std::array <T, N> data_;
public:
    vec() {
        for (auto elem: data_) {
            elem = T();
        }
    }
    T& operator[](const size_t i) {
        assert(i < N);
        return data_[i];
    }
    const T& operator[](const size_t i) const {
        assert(i < N);
        return data_[i];
    }
};

typedef vec<3, float> Vec3f;
typedef vec<4, float> Vec4f;

template <typename T> struct vec<2,T> {
    vec() : x(T()), y(T()) {}
    vec(T X, T Y) : x(X), y(Y) {}
    template <class U> vec<2,T>(const vec<2,U> &v);
    T& operator[](const size_t i)       { assert(i<2); return i<=0 ? x : y; }
    const T& operator[](const size_t i) const { assert(i<2); return i<=0 ? x : y; }
    T x,y;
};

template <typename T>
struct vec<3,T> {
    T x,y,z;

    vec() : x(T()), y(T()), z(T()) {}
    vec(T X, T Y, T Z) : x(X), y(Y), z(Z) {}

    T& operator[](const size_t i) {
        assert(i<3);
        return i<=0 ? x : (1==i ? y : z);
    }
    const T& operator[](const size_t i) const {
        assert(i<3);
        return i<=0 ? x : (1==i ? y : z);
    }
    float norm() const {
        return std::sqrt(x*x+y*y+z*z);
    }
    vec<3,T> & normalize(T l=1) {
        *this = (*this)*(l/norm());
        return *this;
    }
};

template <typename T>
struct vec<4,T> {
    T x,y,z,w;

    vec() : x(T()), y(T()), z(T()), w(T()) {}
    vec(T X, T Y, T Z, T W) : x(X), y(Y), z(Z), w(W) {}

    T& operator[](const size_t i) {
        assert(i<4);
        return i<=0 ? x : (1==i ? y : (2==i ? z : w));
    }
    const T& operator[](const size_t i) const {
        assert(i<4);
        return i<=0 ? x : (1==i ? y : (2==i ? z : w));
    }
};

template<size_t N,typename T>
T operator*(const vec<N,T>& left, const vec<N,T>& right) {
    T result = T();
    for (size_t i = 0; i < N; ++i) {
        result += left[i] * right[i];
    }
    return result;
}

template<size_t N,typename T>
vec<N,T> operator+(const vec<N,T> &left, const vec<N,T> &right) {
    vec<N,T> res;
    for (size_t i = 0; i < N; ++i) {
        res[i] = left[i] + right[i];
    }
    return res;
}

template<size_t N,typename T>
vec<N,T> operator-(const vec<N,T> &left, const vec<N,T> &right) {
    vec<N,T> res;
    for (size_t i = 0; i < N; ++i) {
        res[i] = left[i] - right[i];
    }
    return res;
}

template<size_t N,typename T,typename U>
vec<N,T> operator*(const vec<N,T> &left, const U& c) {
    vec<N,T> res;
    for (size_t i=N; i--; );
    for (size_t i = 0; i < N; ++i) {
        res[i] = left[i] * c;
    }
    return res;
}

template<size_t N,typename T>
vec<N,T> operator-(const vec<N,T> &v) {
    return v * T(-1);
}

template <typename T>
vec<3,T> cross(vec<3,T> v1, vec<3,T> v2) {
    return vec<3,T>(v1.y*v2.z - v1.z*v2.y, v1.z*v2.x - v1.x*v2.z, v1.x*v2.y - v1.y*v2.x);
}

template <size_t N, typename T>
std::ostream& operator<<(std::ostream& out, const vec<N,T>& v) {
    for(unsigned int i = 0; i < N; ++i) {
        out << v[i] << " ";
    }
    return out ;
}

//вектор отражения
Vec3f reflect(const Vec3f &I, const Vec3f &N) {
    return I - N*(I*N)*2;
}

//вектор преломления по закону Снелла
Vec3f refract(const Vec3f &I, const Vec3f &N, const float eta_t, const float eta_i=1.f) {
    float acos = - std::max(-1.f, std::min(1.f, I*N));
    if (acos < 0) {
        return refract(I, -N, eta_i, eta_t); // if the ray comes from the inside the object, swap the air and the media
    }
    float eta = eta_i / eta_t;
    float k = 1 - eta*eta*(1 - acos*acos);
    return k<0 ? Vec3f(1,0,0) : I*eta + N*(eta*acos - sqrtf(k)); // k<0 = total reflection, no ray to refract. I refract it anyways, this has no physical meaning
}

Vec3f crossProd(const Vec3f &a, const Vec3f &b)
{
    return Vec3f(
            a.y * b.z - a.z * b.y,
            a.z * b.x - a.x * b.z,
            a.x * b.y - a.y * b.x
    );
}

struct Light {
    Vec3f position;
    float intensity;

    Light(const Vec3f &p, float i): position(p), intensity(i){}
};

struct Material {
    float refractive_index;
    Vec4f albedo;
    Vec3f  diffuse_color;
    float specular_exponent;
    Material(const float &r, const Vec4f &a, const Vec3f &color, const float &spec):
        refractive_index(r), albedo(a), diffuse_color(color), specular_exponent(spec) {}
    Material():
        refractive_index(1), albedo(1,0,0, 0), diffuse_color(), specular_exponent() {}
};

class Object
{
public:
    Object() {}
    virtual ~Object() {}
    // Method to compute the intersection of the object with a ray
    // Returns true if an intersection was found, false otherwise
    // See method implementation in children class for details
    virtual bool intersect(const Vec3f &, const Vec3f &, Vec3f &, Vec3f &, float &) const = 0;
};

struct Sphere: public Object {
    Material material;
    Vec3f center;
    float radius;

    Sphere(const Vec3f &c, const float &r, const Material &m): center(c), radius(r), material(m) {}
    bool ray_intersect1(const Vec3f &orig, const Vec3f &dir, float &d)  const {
        Vec3f vpc = center - orig; // from orig to center
        Vec3f pc = dir * (vpc * dir);//
        float dist = vpc * vpc - pc * pc;
        return dist <= radius * radius;
    }

    bool intersect(const Vec3f &orig, const Vec3f &dir, Vec3f &N, Vec3f &hit,float &t) const override
    {
        float t0, t1; // solutions for t if the ray intersects
        // geometric solution
        Vec3f L = center - orig;
        float tca = L * dir;
        if (tca < 0) {
            return false;
        }
        float d2 = L * L - tca * tca;
        if (d2 > radius * radius) {
            return false;
        }
        float thc = sqrtf(radius * radius - d2);
        t0 = tca - thc;
        t1 = tca + thc;
        if (t0 > t1) {
            std::swap(t0, t1);
        }
        if (t0 < 0) {
            t0 = t1; // if t0 is negative, let's use t1 instead
            if (t0 < 0) return false; // both t0 and t1 are negative
        }
        t = t0;
        hit = hit = orig + dir*t;
        N = (hit - center).normalize();//нормаль в точке пересечения
        return true;
    }
};

struct Triangle: public Object {
    Vec3f v0;
    Vec3f v1;
    Vec3f v2;
    Material material;
    Triangle(const Vec3f &a, const Vec3f &b, const Vec3f &c, const Material &m):
        v0(a), v1(b), v2(c), material(m) {}
    bool intersect(const Vec3f &orig, const Vec3f &dir, Vec3f &tr_N, Vec3f &hit,float &t) const override {
        Vec3f N = crossProd(v1 - v0, v2 - v0).normalize();//нормаль плоскости
        float denom = N * dir;
        if (fabs(denom) < 0.0001) {//плоскость треугольника параллельна лучу
            return false;
        }
        float dist =  (N * orig + N * v0) / denom;//расстояние до точки пересечения
        if (dist < 0) {
            return false;//треугольник позади луча
        }
        Vec3f point = orig + dir * dist;//точка пересечения
        //далее проверяем на наличие в треуголнике относительно каждого ребра
        Vec3f edge0 = v1 - v0;
        Vec3f C0 = point - v0;
        if (N * crossProd(edge0, C0) < 0) {
            return false;
        }
        Vec3f edge1 = v2 - v1;
        Vec3f C1 = point - v1;
        if (N * crossProd(edge1, C1) < 0) {
            return false;
        }

        Vec3f edge2 = v0 - v2;
        Vec3f C2 = point - v2;
        if (N * crossProd(edge2, C2) < 0) {
            return false;
        }
        tr_N = N;
        hit = orig + dir * dist;
        t  = dist;
        return true;
    }
};

struct Plane: public Object {
    Vec3f norm;
    Vec3f p0;//точка плоскости
    Material m1;
    Material m2;
    Plane(const Vec3f &n, const Vec3f &p, const Material &_m1, const Material &_m2):
        norm(n), p0(p), m1(_m1), m2(_m2) {}
    bool intersect(const Vec3f &orig, const Vec3f &dir, Vec3f &N, Vec3f &hit, float &t) const override {
        float denom = norm * dir;
        if (fabs(denom) > 1e-6) {
            t = (p0 - orig)* norm / denom;
            hit = orig + dir * t;
            N = norm;
            return t>0;
        }
        return false;
    }
};

#endif //RT_TRACER_H
