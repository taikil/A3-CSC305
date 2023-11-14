#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cmath>

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

    Vector operator*(double scalar) const
    {
        return Vector(x * scalar, y * scalar, z * scalar);
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

    Vector() : x(0), y(0), z(0) {}
};

// Ray class
struct Ray
{
    Vector origin, direction;

    Ray(const Vector &origin, const Vector &direction) : origin(origin), direction(direction) {}
};

struct Sphere
{
    std::string name;
    Vector position, scale, color;
    double ka, kd, ks, kr, n;

    Sphere(const std::string &name, const Vector &position, const Vector &scale,
           const Vector &color, double ka, double kd, double ks, double kr, double n)
        : name(name), position(position), scale(scale), color(color), ka(ka), kd(kd), ks(ks), kr(kr), n(n) {}
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
    double near, left, right, bottom, top;
    std::vector<Sphere> spheres;
    std::vector<Light> lights;
    Vector backgroundColor, ambientColor;
    std::string outputFilename;

    Scene() : near(0), left(0), right(0), bottom(0), top(0) {}

    // Constructor with parameters
    Scene(double n, double l, double r, double b, double t,
          const std::vector<Sphere> &sph, const std::vector<Light> &lts,
          const Vector &bg, const Vector &ac, const std::string &output)
        : near(n), left(l), right(r), bottom(b), top(t), spheres(sph), lights(lts),
          backgroundColor(bg), ambientColor(ac), outputFilename(output) {}
};

// Function to check if a ray intersects with a sphere
bool intersect(const Ray &ray, const Sphere &sphere, double &t)
{
    Vector oc = ray.origin - sphere.position;
    double a = ray.direction.dot(ray.direction);
    double b = 2.0 * oc.dot(ray.direction);
    double c = oc.dot(oc) - sphere.scale.dot(sphere.scale);
    double discriminant = b * b - 4 * a * c;

    if (discriminant < 0)
    {
        return false;
    }
    else
    {
        double t1 = (-b - std::sqrt(discriminant)) / (2.0 * a);
        double t2 = (-b + std::sqrt(discriminant)) / (2.0 * a);

        t = (t1 < t2) ? t1 : t2;
        return true;
    }
}

// Function to trace rays and calculate pixel color
Vector trace(const Ray &ray, const std::vector<Sphere> &spheres, const std::vector<Light> &lights, const Vector &ambient)
{
    Vector color(0, 0, 0);
    double t;
    const double epsilon = 0.001;

    for (const auto &sphere : spheres)
    {
        if (intersect(ray, sphere, t) && t > epsilon)
        {
            Vector hitPoint = ray.origin + ray.direction * t;
            Vector normal = (hitPoint - sphere.position).normalize();
            Vector lightDir(0, 0, 0);

            for (const auto &light : lights)
            {
                lightDir = (light.position - hitPoint).normalize();
                double lambertian = std::max(0.0, normal.dot(lightDir));
                color = color + (sphere.color * lambertian) * light.color;
            }

            color = color + ambient * sphere.color;
        }
    }

    return color;
}

// Function to split a string into tokens
std::vector<std::string> split(const std::string &str, char delimiter)
{
    std::vector<std::string> tokens;
    std::istringstream tokenStream(str);
    std::string token;
    while (std::getline(tokenStream, token, delimiter))
    {
        tokens.push_back(token);
    }
    return tokens;
}

Scene readInputFile(const std::string &filename)
{
    std::ifstream inputFile(filename);

    if (!inputFile.is_open())
    {
        std::cerr << "Error opening the file: " << filename << std::endl;
        exit(1);
    }

    Scene scene;
    std::string line;

    while (std::getline(inputFile, line))
    {
        std::vector<std::string> tokens = split(line, ' ');

        if (tokens.size() >= 2)
        {
            if (tokens[0] == "NEAR")
            {
                scene.near = std::stod(tokens[1]);
            }
            else if (tokens[0] == "LEFT")
            {
                scene.left = std::stod(tokens[1]);
            }
            else if (tokens[0] == "RIGHT")
            {
                scene.right = std::stod(tokens[1]);
            }
            else if (tokens[0] == "BOTTOM")
            {
                scene.bottom = std::stod(tokens[1]);
            }
            else if (tokens[0] == "TOP")
            {
                scene.top = std::stod(tokens[1]);
            }
            else if (tokens[0] == "SPHERE")
            {
                // Parse sphere parameters and create a Sphere object
                Sphere sphere(tokens[1],
                              Vector(std::stod(tokens[2]), std::stod(tokens[3]), std::stod(tokens[4])),
                              Vector(std::stod(tokens[5]), std::stod(tokens[6]), std::stod(tokens[7])),
                              Vector(std::stod(tokens[8]), std::stod(tokens[9]), std::stod(tokens[10])),
                              std::stod(tokens[11]), std::stod(tokens[12]), std::stod(tokens[13]),
                              std::stod(tokens[14]), std::stod(tokens[15]));
                scene.spheres.push_back(sphere);
            }
            else if (tokens[0] == "LIGHT")
            {
                // Parse light parameters and create a Light object
                Light light(tokens[1],
                            Vector(std::stod(tokens[2]), std::stod(tokens[3]), std::stod(tokens[4])),
                            Vector(std::stod(tokens[5]), std::stod(tokens[6]), std::stod(tokens[7])));
                scene.lights.push_back(light);
            }
            else if (tokens[0] == "BACK")
            {
                scene.backgroundColor = Vector(std::stod(tokens[1]), std::stod(tokens[2]), std::stod(tokens[3]));
            }
            else if (tokens[0] == "AMBIENT")
            {
                scene.ambientColor = Vector(std::stod(tokens[1]), std::stod(tokens[2]), std::stod(tokens[3]));
            }
            else if (tokens[0] == "OUTPUT")
            {
                scene.outputFilename = tokens[1];
            }
        }
    }

    inputFile.close();
    return scene;
}

// print input file info (for testing purposes)
void printScene(const Scene &scene)
{
    std::cout << "NEAR: " << scene.near << std::endl;
    std::cout << "LEFT: " << scene.left << std::endl;
    std::cout << "RIGHT: " << scene.right << std::endl;
    std::cout << "BOTTOM: " << scene.bottom << std::endl;
    std::cout << "TOP: " << scene.top << std::endl;

    std::cout << "SPHERES:" << std::endl;
    for (const auto &sphere : scene.spheres)
    {
        std::cout << "  Name: " << sphere.name << std::endl;
        std::cout << "  Position: (" << sphere.position.x << ", " << sphere.position.y << ", " << sphere.position.z << ")" << std::endl;
        std::cout << "  Scale: (" << sphere.scale.x << ", " << sphere.scale.y << ", " << sphere.scale.z << ")" << std::endl;
        std::cout << "  Color: (" << sphere.color.x << ", " << sphere.color.y << ", " << sphere.color.z << ")" << std::endl;
        std::cout << "  Ka: " << sphere.ka << std::endl;
        std::cout << "  Kd: " << sphere.kd << std::endl;
        std::cout << "  Ks: " << sphere.ks << std::endl;
        std::cout << "  Kr: " << sphere.kr << std::endl;
        std::cout << "  N: " << sphere.n << std::endl;
    }
    std::cout << "LIGHTS:" << std::endl;
    for (const auto &light : scene.lights)
    {
        std::cout << "  Name: " << light.name << std::endl;
        std::cout << "  Position: (" << light.position.x << ", " << light.position.y << ", " << light.position.z << ")" << std::endl;
        std::cout << "  Color: (" << light.color.x << ", " << light.color.y << ", " << light.color.z << ")" << std::endl;
    }

    std::cout << "BACKGROUND COLOR: (" << scene.backgroundColor.x << ", " << scene.backgroundColor.y << ", " << scene.backgroundColor.z << ")" << std::endl;
    std::cout << "AMBIENT COLOR: (" << scene.ambientColor.x << ", " << scene.ambientColor.y << ", " << scene.ambientColor.z << ")" << std::endl;
    std::cout << "OUTPUT FILENAME: " << scene.outputFilename << std::endl;
}

using namespace std;

int main(int argc, char *argv[])
{
    Scene scene = readInputFile(argv[1]);
    printScene(scene);

    // Scene setup

    Vector ambient(0.2, 0.2, 0.2);
    Vector backgroundColor(1, 1, 1);

    int width = 600;
    int height = 600;

    // Image creation
    std::ofstream outputFile("testSample.ppm");
    outputFile << "P3\n"
               << width << " " << height << "\n255\n";

    for (int y = 0; y < height; ++y)
    {
        for (int x = 0; x < width; ++x)
        {
            double u = double(x) / width * 2 - 1;
            double v = 1 - double(y) / height * 2;

            Ray ray(Vector(0, 0, 0), Vector(u, v, -1).normalize());
        }
        outputFile << "\n";
    }

    outputFile.close();

    return 0;
}
