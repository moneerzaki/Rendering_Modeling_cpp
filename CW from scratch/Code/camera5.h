// pinhole_camera.h
// this camera file is to improve the already working phong rendering process
// by adding the reflection and refraction function to the phong rendering function
// and also to add the binary rendering function to the camera class

#ifndef PINHOLE_CAMERA_H
#define PINHOLE_CAMERA_H

#include <iostream>
#include <fstream>
#include <nlohmann/json.hpp>
#include "vec3.h"
#include "ray.h"
#include "color.h"
#include "hittable.h"

using json = nlohmann::json;
// std::mutex fileMutex;  // Declare a mutex for file access
std::mutex outputMutex;


class PinholeCamera {
public:
    PinholeCamera(const json& cameraData, const json& backgroundColorData, int nbounces)
    {
        width = cameraData.at("width").get<int>();
        height = cameraData.at("height").get<int>();
        position = Vec3(cameraData.at("position")[0], cameraData.at("position")[1], cameraData.at("position")[2]);
        lookAt = Vec3(cameraData.at("lookAt")[0], cameraData.at("lookAt")[1], cameraData.at("lookAt")[2]);
        upVector = Vec3(cameraData.at("upVector")[0], cameraData.at("upVector")[1], cameraData.at("upVector")[2]);
        fov = cameraData.at("fov").get<double>();
        exposure = cameraData.at("exposure").get<double>();
        backgroundColor = Color(backgroundColorData[0].get<double>(),
                                backgroundColorData[1].get<double>(),
                                backgroundColorData[2].get<double>());
        maxDepth = nbounces;

        initialize();
    }

    // Dummy functions to simulate rendering tasks
    void renderPart(const std::vector<Hittable*>& scene, int startRow, int endRow, std::ofstream& ppmFile) {
        for (int y = endRow; y >= startRow; --y) {
            for (int x = width - 1; x >= 0; --x) {
                double u = double(x) / (width - 1);
                double v = double(y) / (height - 1);

                Ray ray = generateRay(u, v);
                Color color = binaryRayColor(ray, scene);

                int ir = static_cast<int>(255.99 * color.red());
                int ig = static_cast<int>(255.99 * color.green());
                int ib = static_cast<int>(255.99 * color.blue());

                ppmFile << ir << " " << ig << " " << ib << "\n";
            }
        }
    }

    void binaryRenderingMT(const std::vector<Hittable*>& scene, std::string& outputFileName, int numThreads) {
        outputFileName += ".ppm";
        std::ofstream ppmFile(outputFileName);

        if (!ppmFile.is_open()) {
            std::cerr << "Error: Unable to open output file." << std::endl;
            return;
        }

        ppmFile << "P3\n" << width << " " << height << "\n255\n";  // PPM header

        std::vector<std::thread> threads;
        int rowsPerThread = height / numThreads;
        std::atomic<int> currentRow(0);

        for (int i = 0; i < numThreads; ++i) {
            int startRow = currentRow.load();
            int endRow = (i == numThreads - 1) ? 0 : startRow - rowsPerThread;

            threads.emplace_back([this, &scene, startRow, endRow, &ppmFile]() {
                renderPart(scene, startRow, endRow, ppmFile);
            });

            currentRow.store(endRow);
        }

        for (auto& thread : threads) {
            thread.join();
        }

        std::cout << "\nDone." << std::endl;
        ppmFile.close();
    }

    void binaryRendering(const std::vector<Hittable*>& scene, std::string& outputFileName){
        outputFileName += ".ppm";
        std::ofstream ppmFile(outputFileName);

        if (!ppmFile.is_open()) {
            std::cerr << "Error: Unable to open output file." << std::endl;
            return;
        }

        ppmFile << "P3\n" << width << " " << height << "\n255\n";  // PPM header

        for (int y = height - 1; y >= 0; --y) {
            // std::clog << "\rScanlines remaining: " << (y) << ' ' << std::flush;
            for (int x = width-1; x >= 0; --x) {
                double u = double(x) / (width - 1);
                double v = double(y) / (height - 1);

                Ray ray = generateRay(u, v);
                Color color = binaryRayColor(ray, scene);

                int ir = static_cast<int>(255.99 * color.red());
                int ig = static_cast<int>(255.99 * color.green());
                int ib = static_cast<int>(255.99 * color.blue());

                ppmFile << ir << " " << ig << " " << ib << "\n";
            }
        }

        std::cout << "\nDone." << std::endl;
        ppmFile.close();
    }

    // the main Blinn-Phong rendering function (working well but slow)
    void phongRendering(const std::vector<Hittable*>& hittableObjects,
                        const std::vector<Light*>& lights,
                        std::string& binaryOutput) {
        
        binaryOutput += ".ppm";
        initialize();
        std::ofstream ppmFile(binaryOutput);

        if (!ppmFile.is_open()) {
            std::cerr << "Error: Unable to open output file." << std::endl;
            return;
        }

        ppmFile << "P3\n" << width << " " << height << "\n255\n";  // PPM header

        for (int y = height - 1; y >= 0; --y) {
        // for (int y = 0; y < height; ++y) {
            // std::clog << "\rScanlines remaining: " << (y) << ' ' << std::flush;
                for (int x = width-1; x >= 0; --x) {
            //     for (int x = 0; x < width; ++x) {
                    double u = double(x) / (width - 1);
                    double v = double(y) / (height - 1);

                    Ray ray = generateRay(u, v);
                    Color color = phongRayColor(ray, hittableObjects, lights);

                    int ir = static_cast<int>(255.99 * color.red());
                    int ig = static_cast<int>(255.99 * color.green());
                    int ib = static_cast<int>(255.99 * color.blue());

                    // std::lock_guard<std::mutex> lock(outputMutex);
                    ppmFile << ir << " " << ig << " " << ib << "\n";
                }
        }

        std::cout << "\nDone." << std::endl;
        ppmFile.close();
    }

    // a try to make the phong rendering function multithreaded 
    void phongRenderingMT(const std::vector<Hittable*>& hittableObjects,
                        const std::vector<Light*>& lights,
                        std::string& binaryOutput) {
        binaryOutput += ".ppm";
        initialize();
        std::ofstream ppmFile(binaryOutput);

        if (!ppmFile.is_open()) {
            std::cerr << "Error: Unable to open output file." << std::endl;
            return;
        }

        ppmFile << "P3\n" << width << " " << height << "\n255\n";  // PPM header

        std::vector<std::future<void>> futures;

        for (int y = height - 1; y >= 0; --y) {
            futures.push_back(std::async(std::launch::async, [this, y, &hittableObjects, &lights, &ppmFile]() {
                for (int x = width - 1; x >= 0; --x) {
                    double u = double(x) / (width - 1);
                    double v = double(y) / (height - 1);

                    Ray ray = generateRay(u, v);
                    Color color = phongRayColor(ray, hittableObjects, lights);

                    int ir = static_cast<int>(255.99 * color.red());
                    int ig = static_cast<int>(255.99 * color.green());
                    int ib = static_cast<int>(255.99 * color.blue());

                    // Protecting the output to the file with a lock
                    std::lock_guard<std::mutex> lock(outputMutex);
                    ppmFile << ir << " " << ig << " " << ib << "\n";
                }
            }));
        }

        // Wait for all futures to complete
        for (auto& future : futures) {
            future.wait();
        }

        std::cout << "\nDone." << std::endl;
        ppmFile.close();
    }


// private:
    int width;
    int height;
    Vec3 position;
    Vec3 lookAt;
    Vec3 upVector;
    double fov;
    double exposure;
    Color backgroundColor;
    int maxDepth;
    
    Vec3 horizontal;
    Vec3 vertical;
    Vec3 lowerLeftCorner;

    void initialize()
    {
        double aspectRatio = static_cast<double>(width) / static_cast<double>(height);
        double theta = degreesToRadians(fov);
        double halfHeight = tan(theta / 2);
        double halfWidth = aspectRatio * halfHeight;

        Vec3 w = (position - lookAt).normalize();
        Vec3 u = upVector.cross(w).normalize();
        Vec3 v = w.cross(u);

        lowerLeftCorner = position - halfWidth * u - halfHeight * v - w;
        horizontal = 2 * halfWidth * u;
        vertical = 2 * halfHeight * v;
    }

    double degreesToRadians(double degrees) const {
        return degrees * (M_PI / 180.0);
    }


    Ray generateRay(double u, double v) const {
        Vec3 direction = lowerLeftCorner + u * horizontal + v * vertical - position;
        return Ray(position, direction.normalize(), 0.0);
    }

    Color binaryRayColor(const Ray& ray, const std::vector<Hittable*>& scene) const {
        HitRecord closestHit;
        double closestHitDistance = std::numeric_limits<double>::infinity();
        Hittable* closestObject = nullptr;

        for (const auto& object : scene) {
            HitRecord currentHit;
            if (object->hit(ray, 0.001, closestHitDistance, currentHit)) {
                if (currentHit.t < closestHitDistance) {
                    closestHitDistance = currentHit.t;
                    closestHit = currentHit;
                    closestObject = object;
                }
                // closestHitDistance = currentHit.t;
                // closestHit = currentHit;
                // closestObject = object;
            }
        }

        if (closestHitDistance < std::numeric_limits<double>::infinity()) {
            return closestObject->binaryColor(); // Use the color of the closest hit
        } else {
            return Color(0.0, 0.0, 0.0); // Black if the ray doesn't hit anything
        }
    }

    // a try to multithread phongraycolor function 
    // Color phongRayColor(const Ray& ray, const std::vector<Hittable*>& objects, const std::vector<Light*>& lights, int depth=0) {
    //     if (depth > maxDepth) {
    //         return Color(0.0, 0.0, 0.0);  // Maximum recursion depth reached
    //     }

    //     HitRecord record;

    //     if (hitAnyObjects(ray, objects, record)) {
    //         // Calculate ambient color
    //         Vec3 ambientColorVec = record.material->ka * record.material->diffusecolor;
    //         Color ambientColor(ambientColorVec.x, ambientColorVec.y, ambientColorVec.z);

    //         // Initialize final color with ambient color
    //         Color finalColor = ambientColor;

    //         // Iterate over lights for diffuse and specular contributions, considering shadows
    //     #pragma omp parallel for
    //         for (int i = 0; i < lights.size(); ++i) {
    //             const auto& light = lights[i];

    //             // Calculate light direction
    //             Vec3 lightDirection = (light->position - record.hit_point).normalize();

    //             // Shadows: Check if there's an object between the hit point and the light source
    //             Ray shadowRay(record.hit_point, lightDirection, 0.0);
    //             HitRecord shadowRecord;
    //             if (!hitAnyObjects(shadowRay, objects, shadowRecord)) {
    //                 Color specularColor(0.0, 0.0, 0.0);
    //                 Color diffuseColor(0.0, 0.0, 0.0);
    //                 if (!mirror) {
    //                     // No shadow, calculate diffuse and specular contributions
    //                     double diffuseFactor = std::max(0.0, record.normal.dot(lightDirection));
    //                     Vec3 diffuseColorVec = record.material->kd * record.material->diffusecolor * light->intensity * diffuseFactor;
    //                     Color temp(diffuseColorVec.x, diffuseColorVec.y, diffuseColorVec.z);
    //                     diffuseColor = temp;

    //                     // Specular reflection (Blinn-Phong model)
    //                     Vec3 viewDirection = -ray.getDirection().normalize();
    //                     Vec3 halfway = (lightDirection + viewDirection).normalize();
    //                     double specularFactor = std::pow(std::max(0.0, record.normal.dot(halfway)), record.material->specularexponent);
    //                     Vec3 specularColorVec = record.material->ks * record.material->specularcolor * light->intensity * specularFactor;
    //                     Color temp2(specularColorVec.x, specularColorVec.y, specularColorVec.z);
    //                     specularColor = temp2;
    //                 }

    //                 // Add diffuse and specular contributions to the final color
    //                 finalColor = finalColor + diffuseColor + specularColor;
    //             }
    //         }

    //         // Recursive reflection
    //         if (record.material->isreflective) {
    //             Vec3 reflectedDirection = Vec3::reflect(ray.getDirection().normalize(), record.normal);
    //             Ray reflectedRay(record.hit_point, reflectedDirection, 0.0);
    //             Color reflectedColor = phongRayColor(reflectedRay, objects, lights, depth + 1);
    //             finalColor = finalColor + record.material->reflectivity * reflectedColor;
    //         }

    //         // Recursive refraction (transparency)
    //         if (record.material->isrefractive) {
    //             Vec3 refractedDirection;
    //             if (Vec3::refract(ray.getDirection().normalize(), record.normal, record.material->refractiveindex, refractedDirection)) {
    //                 Ray refractedRay(record.hit_point, refractedDirection, 0.0);
    //                 Color refractedColor = phongRayColor(refractedRay, objects, lights, depth + 1);
    //                 finalColor = finalColor + (1.0 - record.material->reflectivity) * refractedColor;
    //             }
    //         }

    //         return finalColor;
    //     }

    //     // If no intersection, return the background color
    //     return backgroundColor;
    // }

    // a working function but not so accruate 
    // Color phongRayColor(const Ray& ray, const std::vector<Hittable*>& objects, const std::vector<Light*>& lights) {
    //     HitRecord record;

    //     if (hitAnyObjects(ray, objects, record)) {
    //         // Calculate ambient color
    //         Vec3 ambientColorVec = record.material->ka * record.material->diffusecolor;
    //         Color ambientColor(ambientColorVec.x, ambientColorVec.y, ambientColorVec.z);

    //         // Initialize final color with ambient color
    //         Color finalColor = ambientColor;

    //         // Iterate over lights for diffuse and specular contributions
    //         for (const auto& light : lights) {
    //             // Calculate light direction
    //             Vec3 lightDirection = (light->position - record.hit_point).normalize();

    //             // Diffuse reflection
    //             double diffuseFactor = std::max(0.0, record.normal.dot(lightDirection));
    //             Vec3 diffuseColorVec = record.material->kd * record.material->diffusecolor * light->intensity * diffuseFactor;
    //             Color diffuseColor(diffuseColorVec.x, diffuseColorVec.y, diffuseColorVec.z);

    //             // Specular reflection (Blinn-Phong model)
    //             Vec3 viewDirection = -ray.getDirection().normalize();
    //             Vec3 halfway = (lightDirection + viewDirection).normalize();
    //             double specularFactor = std::pow(std::max(0.0, record.normal.dot(halfway)), record.material->specularexponent);
    //             Vec3 specularColorVec = record.material->ks * record.material->specularcolor * light->intensity * specularFactor;
    //             Color specularColor(specularColorVec.x, specularColorVec.y, specularColorVec.z);

    //             // Add diffuse and specular contributions to the final color
    //             finalColor = finalColor + diffuseColor + specularColor;
    //         }

    //         return finalColor;
    //     }

    //     // If no intersection, return the background color
    //     return backgroundColor;
    // }

    // another try for the function that reflect very good but (bugs)
    bool mirror = false; 
    Color phongRayColor(const Ray& ray, const std::vector<Hittable*>& objects, const std::vector<Light*>& lights, int depth = 0) {
        if (depth > maxDepth) {
            return Color(0.0, 0.0, 0.0);  // Maximum recursion depth reached
        }

        HitRecord record;

        if (hitAnyObjects(ray, objects, record)) {
            // Calculate ambient color
            Vec3 ambientColorVec = record.material->ka * record.material->diffusecolor;
            Color ambientColor(ambientColorVec.x, ambientColorVec.y, ambientColorVec.z);

            // Initialize final color with ambient color
            Color finalColor = ambientColor;

            // Iterate over lights for diffuse and specular contributions, considering shadows
            // for (const auto& light : lights) {
            //     // Calculate light direction
            //     Vec3 lightDirection = (light->position - record.hit_point).normalize();

            //     // Shadows: Check if there's an object between the hit point and the light source
            //     // Ray shadowRay(record.hit_point, lightDirection);
            //     Ray shadowRay(record.hit_point, lightDirection, 0.0);
            //     HitRecord shadowRecord;
            //     if (!hitAnyObjects(shadowRay, objects, shadowRecord)) {


            //         Color specularColor(0.0, 0.0, 0.0);
            //         Color diffuseColor(0.0, 0.0, 0.0);
            //         if (!mirror){
            //         // No shadow, calculate diffuse and specular contributions ( normal is the only variable that changes)
            //         double diffuseFactor = std::max(0.0, record.normal.dot(lightDirection));
            //         Vec3 diffuseColorVec = record.material->kd * record.material->diffusecolor * light->intensity * diffuseFactor;
            //         // diffuseColor(diffuseColorVec.x, diffuseColorVec.y, diffuseColorVec.z);
            //         Color temp(diffuseColorVec.x, diffuseColorVec.y, diffuseColorVec.z);
            //         diffuseColor = temp;

            //         // Specular reflection (Blinn-Phong model) ( normal is the only variable that changes)
                    
            //         Vec3 viewDirection = -ray.getDirection().normalize();
            //         Vec3 halfway = (lightDirection + viewDirection).normalize();
            //         double specularFactor = std::pow(std::max(0.0, record.normal.dot(halfway)), record.material->specularexponent);
            //         Vec3 specularColorVec = record.material->ks * record.material->specularcolor * light->intensity * specularFactor;
            //         Color temp2(specularColorVec.x, specularColorVec.y, specularColorVec.z);
            //         specularColor = temp2;
            //         }

            //         // Add diffuse and specular contributions to the final color
            //         finalColor = finalColor + diffuseColor + specularColor;
            //         // finalColor = finalColor + specularColor;
            //         // finalColor = finalColor + diffuseColor;
            //     }
            // }
            #pragma omp parallel for reduction(+:finalColor) schedule(dynamic)
            for (int i = 0; i < lights.size(); ++i) {
                const auto& light = lights[i];

                // Calculate light direction
                Vec3 lightDirection = (light->position - record.hit_point).normalize();

                // Shadows: Check if there's an object between the hit point and the light source
                Ray shadowRay(record.hit_point, lightDirection, 0.0);
                HitRecord shadowRecord;
                if (!hitAnyObjects(shadowRay, objects, shadowRecord)) {
                    Color specularColor(0.0, 0.0, 0.0);
                    Color diffuseColor(0.0, 0.0, 0.0);
                    if (!mirror) {
                        // No shadow, calculate diffuse and specular contributions (normal is the only variable that changes)
                        double diffuseFactor = std::max(0.0, record.normal.dot(lightDirection));
                        Vec3 diffuseColorVec = record.material->kd * record.material->diffusecolor * light->intensity * diffuseFactor;
                        Color temp(diffuseColorVec.x, diffuseColorVec.y, diffuseColorVec.z);
                        diffuseColor = temp;

                        // Specular reflection (Blinn-Phong model) (normal is the only variable that changes)
                        Vec3 viewDirection = -ray.getDirection().normalize();
                        Vec3 halfway = (lightDirection + viewDirection).normalize();
                        double specularFactor = std::pow(std::max(0.0, record.normal.dot(halfway)), record.material->specularexponent);
                        Vec3 specularColorVec = record.material->ks * record.material->specularcolor * light->intensity * specularFactor;
                        Color temp2(specularColorVec.x, specularColorVec.y, specularColorVec.z);
                        specularColor = temp2;
                    }

                    // Add diffuse and specular contributions to the final color
                    finalColor = finalColor + diffuseColor + specularColor;
                }
            }

            // old Recursive reflection method
            if (record.material->isreflective) {
                Vec3 reflectedDirection = Vec3::reflect(ray.getDirection().normalize(), record.normal);
                // Ray reflectedRay(record.hit_point, reflectedDirection);
                Ray reflectedRay(record.hit_point, reflectedDirection, 0.0);
                Color reflectedColor = phongRayColor(reflectedRay, objects, lights, depth + 1);
                finalColor = finalColor + record.material->reflectivity * reflectedColor;
            }

            // a try for a fersnel reflection method. 
            // New Check if the material is reflective 
            // if (record.material->isreflective && depth < maxDepth) {
            //     // Calculate Fresnel reflection coefficient
            //     double reflectance = Vec3::fresnelReflection(record.normal, ray.getDirection(), record.material->refractiveindex);

            //     // Reflection direction
            //     Vec3 reflectedDirection = Vec3::reflect(ray.getDirection().normalize(), record.normal);

            //     // Recursive reflection
            //     Color reflectedColor = reflectance * phongRayColor(Ray(record.hit_point, reflectedDirection, 0.0), objects, lights, depth + 1);

            //     // Combine with the existing color
            //     finalColor = finalColor + reflectedColor;
            // }

            // old Recursive refraction (transparency)
            if (record.material->isrefractive) {
                Vec3 refractedDirection;
                if (Vec3::refract(ray.getDirection().normalize(), record.normal, record.material->refractiveindex, refractedDirection)) {
                    // Ray refractedRay(record.hit_point, refractedDirection);
                    Ray refractedRay(record.hit_point, refractedDirection, 0.0 /* refractedRay starts outside the object */);
                    Color refractedColor = phongRayColor(refractedRay, objects, lights, depth + 1);
                    finalColor = finalColor + (1.0 - record.material->reflectivity) * refractedColor;
                }
            }

            // a try for a refraction method. 
            // Check if the material is refractive
            // if (record.material->isrefractive && depth < maxDepth) {
            //     // Calculate Fresnel reflection coefficient
            //     double reflectance = Vec3::fresnelReflection(record.normal, ray.getDirection(), record.material->refractiveindex);

            //     // Reflection direction
            //     Vec3 reflectedDirection = Vec3::reflect(ray.getDirection().normalize(), record.normal);

            //     // Recursive reflection
            //     Color reflectedColor = reflectance * phongRayColor(Ray(record.hit_point, reflectedDirection, 0.0), objects, lights, depth + 1);

            //     // Refraction direction
            //     Vec3 refractedDirection;
            //     double transparency = 1.0 - reflectance;
            //     if (Vec3::refract(ray.getDirection().normalize(), record.normal, record.material->refractiveindex, refractedDirection)) {
            //         // Recursive refraction
            //         Color refractedColor = transparency * phongRayColor(Ray(record.hit_point, refractedDirection, 0.0), objects, lights, depth + 1);

            //         // Combine with the existing color
            //         finalColor = finalColor + reflectedColor + refractedColor;
            //     } else {
            //         // Total internal reflection, only reflect
            //         finalColor = finalColor + reflectedColor;
            //     }
            // }

            return finalColor;
        }

        // If no intersection, return the background color
        return backgroundColor;
    }

    // another try for the function for the new texture ... 
    // Color phongRayColor(const Ray& ray, const std::vector<Hittable*>& objects, const std::vector<Light*>& lights, int depth = 0) {
    //     if (depth > maxDepth) {
    //         return Color(0.0, 0.0, 0.0);  // Maximum recursion depth reached
    //     }

    //     HitRecord record;

    //     if (hitAnyObjects(ray, objects, record)) {
    //         // Calculate ambient color
    //         Vec3 ambientColorVec = record.material->ka * record.material->diffusecolor;
    //         Color ambientColor(ambientColorVec.x, ambientColorVec.y, ambientColorVec.z);

    //         // Initialize final color with ambient color
    //         Color finalColor = ambientColor;

    //         // If the material has a texture, get the color from the wooden texture
    //         // Vec3 woodenTextureColor(0,0,0); 
    //         // if (record.material->texture != nullptr) {
    //         //     woodenTextureColor = getWoodenTextureColor(record);
    //         // }

    //         // Iterate over lights for diffuse and specular contributions, considering shadows
    //         for (const auto& light : lights) {
    //             // Calculate light direction
    //             Vec3 lightDirection = (light->position - record.hit_point).normalize();

    //             // Shadows: Check if there's an object between the hit point and the light source
    //             // Ray shadowRay(record.hit_point, lightDirection);
    //             Ray shadowRay(record.hit_point, lightDirection, 0.0);
    //             HitRecord shadowRecord;
    //             if (!hitAnyObjects(shadowRay, objects, shadowRecord)) {
    //                 // No shadow, calculate diffuse and specular contributions ( normal is the only variable that changes)
    //                 double diffuseFactor = std::max(0.0, record.normal.dot(lightDirection));
    //                 Vec3 diffuseColorVec = record.material->kd * record.material->diffusecolor * light->intensity * diffuseFactor;
    //                 Color diffuseColor(diffuseColorVec.x, diffuseColorVec.y, diffuseColorVec.z);

    //                 // Specular reflection (Blinn-Phong model) ( normal is the only variable that changes   )
    //                 Vec3 viewDirection = -ray.getDirection().normalize();
    //                 Vec3 halfway = (lightDirection + viewDirection).normalize();
    //                 double specularFactor = std::pow(std::max(0.0, record.normal.dot(halfway)), record.material->specularexponent);
    //                 Vec3 specularColorVec = record.material->ks * record.material->specularcolor * light->intensity * specularFactor;
    //                 Color specularColor(specularColorVec.x, specularColorVec.y, specularColorVec.z);

    //                 // Add diffuse and specular contributions to the final color
    //                 // finalColor = finalColor + diffuseColor + specularColor;
    //                 finalColor = finalColor + diffuseColor + specularColor;
    //                 // add to the final color the color of the texture if the object has a texture
    //                 if (record.material->texture != nullptr) {
    //                     Vec3 textureColor(record.material->texture->getColor(record.uv.x, record.uv.y));
    //                     finalColor = finalColor + record.material->texture->getColor(record.uv.x, record.uv.y);
    //                 }
    //             }
    //         }
    //         // old Recursive reflection method
    //         if (record.material->isreflective) {
    //             Vec3 reflectedDirection = Vec3::reflect(ray.getDirection().normalize(), record.normal);
    //             // Ray reflectedRay(record.hit_point, reflectedDirection);
    //             Ray reflectedRay(record.hit_point, reflectedDirection, 0.0);
    //             Color reflectedColor = phongRayColor(reflectedRay, objects, lights, depth + 1);
    //             finalColor = finalColor + record.material->reflectivity * reflectedColor;
    //         }
    //         // old Recursive refraction (transparency)
    //         if (record.material->isrefractive) {
    //             Vec3 refractedDirection;
    //             if (Vec3::refract(ray.getDirection().normalize(), record.normal, record.material->refractiveindex, refractedDirection)) {
    //                 // Ray refractedRay(record.hit_point, refractedDirection);
    //                 Ray refractedRay(record.hit_point, refractedDirection, 0.0 /* refractedRay starts outside the object */);
    //                 Color refractedColor = phongRayColor(refractedRay, objects, lights, depth + 1);
    //                 finalColor = finalColor + (1.0 - record.material->reflectivity) * refractedColor;
    //             }
    //         }
    //         return finalColor;
    //     }
    //     // If no intersection, return the background color
    //     return backgroundColor;
    // }

    // bool done = false; // a bool to check if the ray hit any object or not
    bool hitAnyObjects(const Ray& ray, const std::vector<Hittable*>& objects, HitRecord& record) {
        
        mirror = false; 
        bool hitAnything = false;
        // hitRecord closestHit;
        double closestSoFar = std::numeric_limits<double>::infinity();
        double closestHitDistance = std::numeric_limits<double>::infinity();
        Hittable* closestObject = nullptr;

        // a for loop to consider the nearest hit 
        for (const auto& object : objects) {
            HitRecord currentHit;
            if (object->hit(ray, 0.001, closestHitDistance, currentHit)) {
                hitAnything = true;
                if (currentHit.t < closestHitDistance) {
                    closestHitDistance = currentHit.t;
                    record = currentHit;
                    closestObject = object;
                }
            }
        }
        // check if closestobject was a triangle which inherits from hittable
        if (closestObject != nullptr && closestObject->getType() == "triangle"){
            mirror = true; 
        }
        

        return hitAnything;
    }

};

#endif // PINHOLE_CAMERA_H
