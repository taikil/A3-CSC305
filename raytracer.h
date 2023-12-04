#pragma once

#include <vector>
#include <string>

// Define a basic Vector class to represent 3D points and vectors
struct Vector
{
    double x, y, z;

    Vector(double x, double y, double z) : x(x), y(y), z(z) {}

    Vector operator+(const Vector &other) const
    {
        return Vector(x + other.x, y + other.y, z + other.z);
    }

    Vector operator-(const Vector &other) const
    {
        return Vector(x - other.x, y - other.y, z - other.z);
    }

    Vector operator*(const Vector &other) const
    {
        return Vector(x * other.x, y * other.y, z * other.z);
    }

    Vector operator/(const Vector &other) const
    {
        return Vector(x / other.x, y / other.y, z / other.z);
    }

    Vector operator*(double scalar) const
    {
        return Vector(x * scalar, y * scalar, z * scalar);
    }

    Vector &operator+=(const Vector &other)
    {
        x += other.x;
        y += other.y;
        z += other.z;
        return *this;
    }

    Vector &operator=(const Vector &other)
    {
        if (this != &other)
        {
            x = other.x;
            y = other.y;
            z = other.z;
        }
        return *this;
    }

    Vector normalize() const
    {
        double length = std::sqrt(x * x + y * y + z * z);
        return Vector(x / length, y / length, z / length);
    }

    double dot(const Vector &other) const
    {
        return x * other.x + y * other.y + z * other.z;
    }

    Vector clamp() const
    {
        Vector output = Vector(x, y, z);
        if (x > 1)
        {
            output = Vector(1, y, z);
        }
        if (y > 1)
        {
            output = Vector(x, 1, z);
        }
        if (z > 1)
        {
            output = Vector(x, y, 1);
        }
        return output;
    }

    double length() const
    {
        return std::sqrt(x * x + y * y + z * z);
    }

    bool operator==(const Vector &other) const
    {
        return x == other.x && y == other.y && z == other.z;
    }

    bool operator!=(const Vector &other) const
    {
        return x != other.x || y != other.y || z != other.z;
    }

    friend std::ostream &operator<<(std::ostream &os, const Vector &vector);
    Vector() : x(0), y(0), z(0) {}
};

std::ostream &operator<<(std::ostream &os, const Vector &vector)
{
    os << "{ x: " << vector.x << ", y: " << vector.y << ", z: " << vector.z << " }";
    return os;
}

Vector operator*(double scalar, const Vector &v)
{
    return Vector(scalar * v.x, scalar * v.y, scalar * v.z);
}
// Ray class
struct Ray
{
    Vector origin, direction;
    int depth;

    Ray(const Vector &origin, const Vector &direction, int depth) : origin(origin), direction(direction), depth(depth) {}
    friend std::ostream &operator<<(std::ostream &os, const Ray &ray);
};

std::ostream &operator<<(std::ostream &os, const Ray &ray)
{
    os << "{ Origin: " << ray.origin << ", Direction: " << ray.direction << ", z: " << ray.depth << " }";
    return os;
}

struct Sphere
{
    std::string name;
    Vector position, scale, color;
    double ka, kd, ks, kr, n;

    Sphere(const std::string &name, const Vector &position, const Vector &scale,
           const Vector &color, double ka, double kd, double ks, double kr, double n)
        : name(name), position(position), scale(scale), color(color), ka(ka), kd(kd), ks(ks), kr(kr), n(n) {}
    Sphere() : name(""), position(0, 0, 0), scale(0, 0, 0), color(0, 0, 0), ka(0), kd(0), ks(0), kr(0), n(0) {}
};

struct Light
{
    std::string name;
    Vector position, color;

    Light(const std::string &name, const Vector &position, const Vector &color)
        : name(name), position(position), color(color) {}
};

// Scene Class
struct Scene
{
    double near, left, right, bottom, top, width, height;
    std::vector<Sphere> spheres;
    std::vector<Light> lights;
    Vector backgroundColor, ambientColor;
    std::string outputFilename;

    Scene() : near(0), left(0), right(0), bottom(0), top(0) {}

    // Constructor with parameters
    Scene(double n, double l, double r, double b, double t, double w, double h,
          const std::vector<Sphere> &sph, const std::vector<Light> &lts,
          const Vector &bg, const Vector &ac, const std::string &output)
        : near(n), left(l), right(r), bottom(b), top(t), width(w), height(h), spheres(sph), lights(lts),
          backgroundColor(bg), ambientColor(ac), outputFilename(output) {}
};

bool intersect(const Ray &ray, const Sphere &sphere, double &t);
bool closest_intersection(const Ray &ray, const std::vector<Sphere> &spheres, Sphere &closestSphere, double &t, bool &inside);
Vector trace(const Ray &ray, const Scene &scene);
std::vector<std::string> split(const std::string &str, char delimiter);
Scene readInputFile(const std::string &filename);
void writePPM(const std::string &filename, const std::vector<std::vector<Vector>> &pixelValues);
void printScene(const Scene &scene);
