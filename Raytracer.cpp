#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cmath>
#include "raytracer.h"
#include "invert.cpp"

int invert(const Sphere &sphere, double (&inverseT)[4][4])
{
    double T[4][4]{};
    T[0][0] = sphere.scale.x;
    T[1][1] = sphere.scale.y;
    T[2][2] = sphere.scale.z;
    T[0][3] = sphere.position.x;
    T[1][3] = sphere.position.y;
    T[2][3] = sphere.position.z;
    T[3][3] = 1;

    invert_matrix(T, inverseT);
    return 0;
}

// Check if a ray intersects with a sphere
bool intersect(const Ray &ray, const Sphere &sphere, double &t)
{
    double rayLength = 0.0;
    if (ray.origin == Vector(0, 0, 0))
    {
        rayLength = ray.direction.length();
    }

    double inverseT[4][4]{};
    invert(sphere, inverseT);

    Vector invTrans = Vector(inverseT[0][3], inverseT[1][3], inverseT[2][3]);
    Vector invScale = Vector(inverseT[0][0], inverseT[1][1], inverseT[2][2]);
    Vector invOrigin = (ray.origin + invTrans);
    Vector invDirection = ray.direction * invScale;

    double a = invDirection.dot(invDirection);     // |c|^2
    double b = 2.0f * invOrigin.dot(invDirection); // 2 (S * c)
    double c = invOrigin.dot(invOrigin) - 1;       // S^2 - 1
    double discriminant = b * b - 4 * a * c;

    if (discriminant > 0)
    {
        // Find the nearest intersection point
        float t1 = (-b - std::sqrt(discriminant)) / (2.0f * a);
        float t2 = (-b + std::sqrt(discriminant)) / (2.0f * a);

        if (t1 > rayLength && t1 > 0 && (t1 < t2 || (t2 < 0 && t1 >= 0)))
        {
            t = t1;
            return true;
        }

        if (t2 > rayLength && t2 > 0 && (t2 < t1 || (t1 < 0 && t2 >= 0)))
        {
            t = t2;
            return true;
        }
    }

    return false;
}

// Find the closest intersection among spheres
bool closest_intersection(const Ray &ray, const std::vector<Sphere> &spheres, Sphere &closestSphere, double &t)
{
    double minT = std::numeric_limits<double>::infinity();

    for (const auto &sphere : spheres)
    {
        double currentT;

        if (intersect(ray, sphere, currentT))
        {
            if (currentT < minT)
            {
                minT = currentT;
                closestSphere = sphere;
            }
        }
    }

    t = minT;

    return (t < std::numeric_limits<double>::infinity());
}

// Function to perform ADS calculations for a given intersection point
Vector ADS(const Ray &ray, const Scene &scene, const Sphere &closestSphere, const Vector &localIntersection, const Vector &normal)
{
    // Ka*Ia[c]*O[c]
    Vector ambient = scene.ambientColor * closestSphere.ka * closestSphere.color;

    // Base for diffuse and specular
    Vector diffuse(0, 0, 0);
    Vector specular(0, 0, 0);

    for (const auto &light : scene.lights)
    {
        Vector lightDirection = (light.position - localIntersection).normalize();
        float NdotL = std::max(0.0, normal.dot(lightDirection));
        // Shadow ray calculation
        Sphere shadowClosestSphere;
        double shadowT;
        Ray shadowRay = Ray(localIntersection + lightDirection * 0.00001, lightDirection.normalize(), 0);
        // std::cout << "Shadow Ray: " << shadowRay << "\n";

        if (!closest_intersection(shadowRay, scene.spheres, shadowClosestSphere, shadowT) || shadowT > (light.position - localIntersection).length())
        {
            // No intersection with shadow ray or intersection beyond the light source
            // Calculate diffuse and specular terms

            // L[c] * kd * NdotL * O[c]
            diffuse += (light.color * closestSphere.kd * std::max(0.0f, NdotL) * closestSphere.color);

            Vector rDirection = (2.0 * normal.dot(lightDirection.normalize()) * normal - lightDirection.normalize()).normalize();
            float RdotV = std::max(0.0, rDirection.dot((ray.origin - localIntersection).normalize()));

            // L[c] * ks * RdotV^n
            specular += light.color * closestSphere.ks * std::pow(RdotV, closestSphere.n);
        }
        // If there is an intersection with the shadow ray, the point is in shadow, so no contribution.
    }

    Vector totalColor = ambient + diffuse + specular;
    return totalColor;
}

// Trace rays and calculate pixel color
Vector trace(const Ray &ray, const Scene &scene)
{
    int MAX_DEPTH = 3;
    if (ray.depth > MAX_DEPTH)
    {
        return Vector(0, 0, 0); // Return black
    }

    Sphere closestSphere = scene.spheres[0];

    double t;

    if (!closest_intersection(ray, scene.spheres, closestSphere, t) && ray.origin == Vector(0, 0, 0))
    {
        return scene.backgroundColor;
    }
    else if (!closest_intersection(ray, scene.spheres, closestSphere, t) && ray.origin == Vector(0, 0, 0))
    {
        return Vector(0, 0, 0);
    }

    // S + ct_h
    Vector localIntersection = ray.origin + (ray.direction * t);
    // Normal * Inverse Transpose of T
    Vector squared = Vector(closestSphere.scale.x * closestSphere.scale.x, closestSphere.scale.y * closestSphere.scale.y, closestSphere.scale.z * closestSphere.scale.z);
    Vector normal = ((localIntersection - closestSphere.position) / squared).normalize();

    // Calculate ADS components
    Vector adsColor = ADS(ray, scene, closestSphere, localIntersection, normal);

    // v = −2(N⋅c)⋅N+c.
    Vector reflectionDirection = (-2.0 * normal.dot(ray.direction) * normal + ray.direction).normalize();

    // Create the reflected ray
    Ray reflectedRay = Ray(localIntersection + reflectionDirection * 0.00001, reflectionDirection, ray.depth + 1);

    // Trace reflections
    // Vector reflectionColor = trace(reflectedRay, scene);
    Vector reflectionColor = Vector(0, 0, 0);

    // Combine ADS and reflections
    Vector totalColor = adsColor + (closestSphere.kr * reflectionColor);
    return totalColor;
}

// // Split file info into tokens
std::vector<std::string> split(const std::string &str, const std::string &delimiters)
{
    std::vector<std::string> tokens;
    std::size_t start = 0, end;

    while ((end = str.find_first_of(delimiters, start)) != std::string::npos)
    {
        // Extract token and remove leading/trailing whitespaces
        std::string token = str.substr(start, end - start);
        token.erase(token.find_last_not_of(" \t") + 1);
        token.erase(0, token.find_first_not_of(" \t"));

        // Add non-empty tokens to the result
        if (!token.empty())
            tokens.push_back(token);

        start = end + 1; // Move to the next character after the delimiter
    }

    // Add the last token if any
    std::string lastToken = str.substr(start);
    lastToken.erase(lastToken.find_last_not_of(" \t") + 1);
    lastToken.erase(0, lastToken.find_first_not_of(" \t"));

    if (!lastToken.empty())
        tokens.push_back(lastToken);

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
        std::vector<std::string> tokens = split(line, " \t");

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
            else if (tokens[0] == "RES")
            {
                scene.width = std::stod(tokens[1]);
                scene.height = std::stod(tokens[2]);
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

void writePPM(const std::string &filename, const std::vector<std::vector<Vector>> &pixelValues)
{
    std::ofstream ppmFile(filename);

    if (ppmFile.is_open())
    {
        // Write PPM header
        int width = pixelValues[0].size();
        int height = pixelValues.size();
        ppmFile << "P3\n"
                << width << " " << height << "\n255\n";

        // Write pixel values
        for (auto it = pixelValues.rbegin(); it != pixelValues.rend(); ++it)
        {
            const auto &row = *it;
            for (const auto &pixel : row)
            {
                int r = static_cast<int>(std::round(pixel.x * 255));
                int g = static_cast<int>(std::round(pixel.y * 255));
                int b = static_cast<int>(std::round(pixel.z * 255));

                // std::cout << "R " << r << " G " << g << " B" << b << std::endl;

                ppmFile << r << " " << g << " " << b << " ";
            }
            ppmFile << "\n";
        }

        ppmFile.close();
        std::cout << "Image generated successfully: " << filename << std::endl;
    }
    else
    {
        std::cerr << "Error: Unable to open file for writing." << std::endl;
    }
}

// print input file info (for testing purposes)
void printScene(const Scene &scene)
{
    std::cout << "NEAR: " << scene.near << std::endl;
    std::cout << "LEFT: " << scene.left << std::endl;
    std::cout << "RIGHT: " << scene.right << std::endl;
    std::cout << "BOTTOM: " << scene.bottom << std::endl;
    std::cout << "TOP: " << scene.top << std::endl;
    std::cout << "RES: " << scene.width << " x " << scene.height << std::endl;

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

int main(int argc, char *argv[])
{
    if (argc != 2)
    {
        std::cout << "Error, missing filename, call the command as \"./Raytracer [filename.txt]\"\n";
        exit(1);
    }
    Scene scene = readInputFile(argv[1]);
    printScene(scene);
    Vector eye(0, 0, 0);
    std::vector<std::vector<Vector>> pixelValues(scene.height, std::vector<Vector>(scene.width));

    // Loop through each pixel
    for (int r = 0; r < scene.height; ++r)
    {
        // -W + W2c/nCols
        float v = -scene.right + (2.0 * scene.right * r) / scene.height;
        for (int c = 0; c < scene.width; ++c)
        {

            // -H + H2r/nRows
            float u = -scene.top + (2.0 * scene.top * c) / scene.width;
            // Calculate the ray direction
            Vector pixelDirection = Vector(u, v, -scene.near);
            Ray ray(eye, pixelDirection, 0);
            // std::cout << ray << std::endl;
            // std::cout << "Pixel (x,y): (" << r << ", " << c << ")\n";
            Vector color = trace(ray, scene);
            pixelValues[r][c] = color;
        }
    }

    writePPM(scene.outputFilename, pixelValues);

    return 0;
}
