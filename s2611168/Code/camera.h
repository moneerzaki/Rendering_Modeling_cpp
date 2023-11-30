// camera.h
#ifndef CAMERA_H
#define CAMERA_H

#include <iostream>
#include <fstream>
#include <nlohmann/json.hpp>
#include "json_read.h"
#include "vec3.h"
#include "ray.h"
#include "hittable.h"
#include "color.h"
#include "basics.h"


using json = nlohmann::json;

class PinholeCamera {
public:
    // Constructor that takes a JSON object and initializes the camera
    PinholeCamera(const json& cameraData, const json& backgroundColor_json) :
        width(cameraData.at("width").get<int>()),
        height(cameraData.at("height").get<int>()),
        position(Vec3(cameraData.at("position")[0], cameraData.at("position")[1], cameraData.at("position")[2])),
        lookAt(Vec3(cameraData.at("lookAt")[0], cameraData.at("lookAt")[1], cameraData.at("lookAt")[2])),
        upVector(Vec3(cameraData.at("upVector")[0], cameraData.at("upVector")[1], cameraData.at("upVector")[2])),
        fov(cameraData.at("fov").get<double>()),
        exposure(cameraData.at("exposure").get<double>())
        // Initialize other camera parameters...

        // backgroundColor(backgroundColor[0], backgroundColor[1], backgroundColor[2])
        // backgroundColor(backgroundColor[0], backgroundColor[1], backgroundColor[2])
    {

        // overwrite width and height 
        width = 10000;
        height = 10000;
        // Resize the pixels vector to match the image resolution
        pixels.resize(height, std::vector<Vec3>(width, Vec3(0.0, 0.0, 0.0)));

        aspectRatio = static_cast<double>(width) / height;

        // Example (simplified):
        Vec3 w = (position - lookAt).normalize();
        Vec3 u = upVector.cross(w).normalize();
        Vec3 v = w.cross(u);

        double halfHeight = tan(degreesToRadians(fov) / 2);
        double halfWidth = aspectRatio * halfHeight;

        // upperLeftCorner = origin + halfWidth * u + halfHeight * v + w;
        // horizontal = -2 * halfWidth * u;
        // vertical = -2 * halfHeight * v;

        upperLeftCorner = position - halfWidth * u - halfHeight * v - w;
        horizontal = 2 * halfWidth * u;
        vertical = 2 * halfHeight * v;

        // Other initialization...
        // backgroundColor = Color(backgroundColor[0].get<double>(), backgroundColor[1].get<double>(), backgroundColor[2].get<double>());
        backgroundColor = Color(backgroundColor_json);

        // overwrite for width and heigth 
        width = 10000;
        height = 10000;
    }

    // binary rendering the scene using the camera
    void binaryRendering(const std::vector<Hittable*>& objects, const std::string& outputFileName) {
        std::ofstream ppmFile(outputFileName);

        if (!ppmFile.is_open()) {
            std::cerr << "Error: Unable to open output file." << std::endl;
            return;
        }

        ppmFile << "P3\n" << width << " " << height << "\n255\n";  // PPM header

        for (int y = height - 1; y >= 0; --y) {
            // std::cout << "\rScanlines remaining: " << y << ' ' << std::flush;
            for (int x = 0; x < width; ++x) {

                double u = double(x) / (width - 1);
                double v = double(y) / (height - 1);

                Ray ray = generateRay(u, v, 0.0);
                // Ray ray = get_ray(x, y, 1, 1);
                Color color = binaryRayColor(ray, objects); 

                // Color color = binaryRayColor(ray, scene);
                int ir = static_cast<int>(255.99 * color.red());
                int ig = static_cast<int>(255.99 * color.green());
                int ib = static_cast<int>(255.99 * color.blue());

                ppmFile << ir << " " << ig << " " << ib << "\n";
            }
        }

        std::cout << "\nDone." << std::endl;
        // ppmFile.flush();
        ppmFile.close();
    }

    // first method of ray generation
    // Ray generateRay(double x, double y, double time) const {
    //     // Generate a ray from the camera through the specified (x, y) coordinates on the image plane
    //     Vec3 direction = (upperLeftCorner + x * horizontal + y * vertical - origin).normalize();
    //     return Ray(origin, direction, time);
    // }


    // another method of generating rays 
    Ray generateRay(double x, double y, double time) const {
        // Convert (x, y) from pixel coordinates to normalized device coordinates
        double ndcX = (2.0 * x - width) / width;
        double ndcY = (height - 2.0 * y) / height;

        // Calculate the direction vector in camera space
        Vec3 directionInCameraSpace(ndcX, ndcY, -1.0);

        // Transform the direction to world space
        Vec3 direction = transformDirection(directionInCameraSpace);

        return Ray(position, direction.normalize(), time);
        // return Ray(origin, direction.normalize(), time);
    }

    // a thrid method of ray generation.
    // Ray generateRay(double x, double y, double time) const {
    //     // Convert (x, y) from pixel coordinates to normalized device coordinates
    //     double ndcX = (2.0 * x - width) / width;
    //     double ndcY = (height - 2.0 * y) / height;

    //     // Calculate the direction vector in camera space
    //     Vec3 directionInCameraSpace(ndcX, ndcY, -1.0);

    //     // Transform the direction to world space
    //     Vec3 direction = transformDirection(directionInCameraSpace);

    //     return Ray(position, direction.normalize(), time);
    // }

    Vec3 transformDirection(const Vec3& directionInCameraSpace) const {
        // Transform the direction from camera space to world space
        Vec3 w = (position - lookAt).normalize();
        Vec3 u = upVector.cross(w).normalize();
        Vec3 v = w.cross(u);

        // Calculate the world space direction
        Vec3 worldSpaceDirection = (u * directionInCameraSpace.x) +
                                (v * directionInCameraSpace.y) +
                                (w * directionInCameraSpace.z);

        return worldSpaceDirection;
    }



    void print() const {
        std::cout << "********PinholeCamera********" << std::endl;
        std::cout << "Width: " << width << std::endl;
        std::cout << "Height: " << height << std::endl;
        std::cout << "Position: " << position << std::endl;
        std::cout << "LookAt: " << lookAt << std::endl;
        std::cout << "UpVector: " << upVector << std::endl;
        std::cout << "FOV: " << fov << std::endl;
        std::cout << "Exposure: " << exposure << std::endl;
        std::cout << "Background Color: "; 
        backgroundColor.print();
        std::cout << std::endl;
        std::cout << "*****************************" << std::endl;
    }

private:
    // Private member variables
    double aspectRatio;
    int width;
    int height;
    Vec3 position;
    Vec3 lookAt;
    Vec3 upVector;
    double fov;
    double exposure;
    Color backgroundColor;
    std::vector<std::vector<Vec3>> pixels;

    // Add more private fields as needed
    Vec3 origin;
    Vec3 upperLeftCorner;
    Vec3 horizontal;
    Vec3 vertical;

    // Helper function to convert degrees to radians
    double degreesToRadians(double degrees) const {
        return degrees * (M_PI / 180.0);
    }
    
    // Update the binaryRayColor function
    Color binaryRayColor(const Ray& ray, const std::vector<Hittable*>& objects) {
        HitRecord record;
        for (const auto& object : objects) {
            if (object->hit(ray, 0.001, std::numeric_limits<double>::infinity(), record)) {
                // If the ray hits an object, return its binary color
                return object->binaryColor();
            }
        }
        // If no intersection, return the default background color
        return backgroundColor;
    }

    bool hitAnyObjects(const Ray& ray, const std::vector<Hittable*>& objects, HitRecord& record) const {
        bool hitAnything = false;
        double closestSoFar = DBL_MAX;

        for (const auto& object : objects) {
            if (object->hit(ray, 0.001, closestSoFar, record)) {
                hitAnything = true;
                closestSoFar = record.t;
            }
        }

        return hitAnything;
    }

    // Other camera parameters...

    // Add more private fields as needed
    // ...
};


#endif // CAMERA_H
