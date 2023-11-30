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
#include "basics.h"

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
        samplesPerPixel = 1; // default value
        apertureRadius = 0.0001; // default value
        // apertureRadius = 0.1; // default value

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


    Color toneMap(const Color& color, double exposure) {
        // Reinhard tone mapping operator
        Color mappedColor = color * (exposure / (color + 1.0));
        return mappedColor;
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
            // std::clog << "\rScanlines remaining: " << (y) << ' ' << std::flush;
            for (int x = width - 1; x >= 0; --x) {
                Color finalColor(0.0, 0.0, 0.0);

                // Multi-sampling loop
                for (int s = 0; s < samplesPerPixel; ++s) {
                    double u = (x + random_double()) / (width - 1);
                    double v = (y + random_double()) / (height - 1);

                    Ray ray = generateRay(u, v);
                    Color color = phongRayColor(ray, hittableObjects, lights);

                    finalColor = finalColor + color;
                }

                // Average color over all samples
                finalColor = finalColor / samplesPerPixel;

                // Tone mapping
                finalColor = toneMap(finalColor, 2);

                // Normalize color values to the range [0, 1]
                finalColor.normalize();

                int ir = static_cast<int>(255.99 * finalColor.red());
                int ig = static_cast<int>(255.99 * finalColor.green());
                int ib = static_cast<int>(255.99 * finalColor.blue());

                // Ensure values are in the valid range [0, 255]
                ir = std::min(std::max(ir, 0), 255);
                ig = std::min(std::max(ig, 0), 255);
                ib = std::min(std::max(ib, 0), 255);

                ppmFile << ir << " " << ig << " " << ib << "\n";
            }
        }

        std::cout << "\nDone." << std::endl;
        ppmFile.close();
    }

    // a try to make the phong rendering function multithreaded 
    // void renderRows(const Hittable& hittableObjects,
    //                 const std::vector<Light*>& lights,
    //                 const std::string& binaryOutput,
    //                 int startRow, int endRow) {
    //     std::ofstream outFile(binaryOutput, std::ios::binary | std::ios::app);

    //     for (int y = startRow; y < endRow; ++y) {
    //         for (int x = 0; x < width; ++x) {
    //             // Assuming you have a castRay function
    //             Ray ray = generateRay(x, y); // Adjust this according to your implementation
    //             Vec3 color = castRay(ray, hittableObjects, lights, maxDepth);

    //             // Assuming you have a writeColor function to convert and write color data
    //             writeColor(outFile, color);
    //         }
    //     }

    //     outFile.close();
    // }
    // // Function to convert and write color data to a binary file
    // Vec3 castRay(const Ray& ray, const HittableList& hittableObjects, const std::vector<Light*>& lights, int depth) const {
    //     if (depth <= 0) {
    //         return Vec3(0.0, 0.0, 0.0); // Return black for rays that exceed max depth
    //     }

    //     HitRecord hitRecord;
    //     if (hittableObjects.hit(ray, 0.001, std::numeric_limits<double>::infinity(), hitRecord)) {
    //         Vec3 emitted = hitRecord.material->emitted(hitRecord.u, hitRecord.v, hitRecord.point);
    //         ScatterRecord scatterRecord;
    //         Vec3 attenuation;
    //         if (hitRecord.material->scatter(ray, hitRecord, attenuation, scatterRecord)) {
    //             if (scatterRecord.isSpecular) {
    //                 return attenuation * castRay(scatterRecord.specularRay, hittableObjects, lights, depth - 1);
    //             }

    //             Vec3 diffuseLight = computeDiffuseLight(hitRecord.point, hitRecord.normal, lights);
    //             Vec3 specularLight = computeSpecularLight(hitRecord.point, hitRecord.normal, ray.direction(), lights, hitRecord.material->shininess());

    //             return emitted + attenuation * (diffuseLight + specularLight);
    //         } else {
    //             return emitted;
    //         }
    //     } else {
    //         // Background color
    //         return backgroundColor;
    //     }
    // }
    // // Function to compute diffuse light
    // Vec3 computeDiffuseLight(const Vec3& point, const Vec3& normal, const std::vector<Light*>& lights) const {
    //     Vec3 diffuseLight(0.0, 0.0, 0.0);
    //     for (const auto& light : lights) {
    //         diffuseLight += light->computeDiffuse(point, normal);
    //     }
    //     return diffuseLight;
    // }
    // // Function to compute specular light
    // Vec3 computeSpecularLight(const Vec3& point, const Vec3& normal, const Vec3& viewDirection, const std::vector<Light*>& lights, double shininess) const {
    //     Vec3 specularLight(0.0, 0.0, 0.0);
    //     for (const auto& light : lights) {
    //         specularLight += light->computeSpecular(point, normal, viewDirection, shininess);
    //     }
    //     return specularLight;
    // }
    // void writeColor(std::ofstream& outFile, const Vec3& color) {
    //     // Assuming Color has values in the range [0.0, 1.0]
    //     uint8_t r = static_cast<uint8_t>(255.999 * color.x);
    //     uint8_t g = static_cast<uint8_t>(255.999 * color.y);
    //     uint8_t b = static_cast<uint8_t>(255.999 * color.z);

    //     // Write color data to the binary file
    //     outFile.put(r);
    //     outFile.put(g);
    //     outFile.put(b);
    // }
    // void phongRenderingMultithreaded(const Hittable& hittableObjects,
    //                                  const std::vector<Light*>& lights,
    //                                  const std::string& binaryOutput) {
    //     // Assuming binaryOutput is a string (adjust the type if needed)

    //     int numThreads = std::thread::hardware_concurrency(); // Use the number of available threads

    //     std::vector<std::thread> threads;
    //     int rowsPerThread = height / numThreads;

    //     for (int i = 0; i < numThreads; ++i) {
    //         int startRow = i * rowsPerThread;
    //         int endRow = (i == numThreads - 1) ? height : (startRow + rowsPerThread);

    //         threads.emplace_back(&PinholeCamera::renderRows, this, hittableObjects, lights, binaryOutput, startRow, endRow);
    //     }

    //     for (auto& thread : threads) {
    //         thread.join();
    //     }
    // }


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
    int samplesPerPixel;
    double apertureRadius;
    std::vector<std::vector<Vec3>> image;

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


    // Ray generateRay(double u, double v) const {
    //     Vec3 direction = lowerLeftCorner + u * horizontal + v * vertical - position;
    //     return Ray(position, direction.normalize(), 0.0);
    // }
    Ray generateRay(double u, double v) const {
        // Generate a ray from the camera position to a point on the image plane
        Vec3 direction = lowerLeftCorner + u * horizontal + v * vertical - position;

        // Sample a point on the camera's aperture
        Vec3 apertureSample = randomInUnitDisk() * apertureRadius;
        Vec3 offset = Vec3(u * apertureSample.x, v * apertureSample.y,0.0);

        return Ray(position + offset, direction.normalize(), 0.0);
    }
    Vec3 randomInUnitDisk() const {
        Vec3 p;
        do {
            p = 2.0 * Vec3(random_double(), random_double(), 0) - Vec3(1, 1, 0);
        } while (p.dot(p) >= 1.0);
        return p;
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

    // another try for the function that reflect very good ( with no refraction with no texture)
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

    //         // Check if the material has a texture
    //         if (record.material->texture != nullptr) {
    //             // Apply texture to the diffuse color
    //             record.material->textureColor = applyWoodenTexture(record.hit_point, *(record.material->texture));
    //         }

    //         #pragma omp parallel for reduction(+:finalColor) schedule(dynamic)
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
    //                     // No shadow, calculate diffuse and specular contributions (normal is the only variable that changes)
    //                     double diffuseFactor = std::max(0.0, record.normal.dot(lightDirection));
    //                     Vec3 diffuseColorVec = record.material->kd * record.material->diffusecolor * light->intensity * diffuseFactor;

    //                     Vec3 v = record.material->textureColor;
    //                     if (v.x != 0 && v.y!=0 && v.z!=0) {
    //                         diffuseFactor = std::max(0.0, record.normal.dot(lightDirection));
    //                         diffuseColorVec = record.material->textureColor * record.material->diffusecolor * light->intensity * diffuseFactor;
    //                     }
    //                     Color temp(diffuseColorVec.x, diffuseColorVec.y, diffuseColorVec.z);
    //                     diffuseColor = temp;

    //                     // Specular reflection (Blinn-Phong model) (normal is the only variable that changes)
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

    //         // gamma correctoin 
    //         // finalColor = applyGammaCorrection(finalColor, 2.2);
    //         return finalColor;
    //     }

    //     // If no intersection, return the background color
    //     return backgroundColor;
    // }
    bool hitAnyObjects(const Ray& ray, const std::vector<Hittable*>& objects, HitRecord& record) {
        
        bool hitAnything = false;
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
        
        return hitAnything;
    }

    Color applyGammaCorrection(const Color& color, double gamma) {
        double r = std::pow(color.red(), 1.0 / gamma);
        double g = std::pow(color.green(), 1.0 / gamma);
        double b = std::pow(color.blue(), 1.0 / gamma);
        return Color(r, g, b);
    }

    Vec3 applyWoodenTexture(const Vec3& hitPoint, const WoodTexture& woodTexture) {
        // Extract texture parameters from the base Texture class
        double scale = woodTexture.scale;
        Vec3 colorMin = woodTexture.colorMin;
        Vec3 colorMax = woodTexture.colorMax;
        double ringSpacing = woodTexture.ringSpacing;
        double ringThickness = woodTexture.ringThickness;
        double noise = woodTexture.noise;

        // Calculate 3D coordinates in texture space
        double s = scale * hitPoint.x;
        double t = scale * hitPoint.z;
        double woodPattern = (sinf(s / ringSpacing) * cosf(t / ringSpacing) + 1.0) / 2.0;
        double ring = fmod(woodPattern, ringThickness);

        // Apply color variations and noise
        Vec3 woodColor = colorMin + (colorMax - colorMin) * woodPattern;
        woodColor = woodColor * (1.0 - noise + 2.0 * noise * (rand() / (double)RAND_MAX));

        return woodColor;
    }

    Vec3 applyTexture(const Vec3& hitPoint, const CustomMaterial& material) {
        if (material.texture != nullptr) {
            // Apply the specific texture to the hit point
            return applyWoodenTexture(hitPoint, *(static_cast<WoodTexture*>(material.texture))); // Assuming WoodTexture inherits from Texture
        }
        return Vec3(1.0, 1.0, 1.0);  // Default to white if no texture
    }

    Color phongRayColor(const Ray& ray, const std::vector<Hittable*>& objects, const std::vector<Light*>& lights, int depth = 0) {
        if (depth > maxDepth) {
            return Color(0.0, 0.0, 0.0);  // Maximum recursion depth reached
        }

        HitRecord record;

        if (hitAnyObjects(ray, objects, record)) {
            // Calculate ambient color
            Vec3 ambientColorVec = record.material->ka * applyTexture(record.hit_point, *(record.material));
            Color ambientColor(ambientColorVec.x, ambientColorVec.y, ambientColorVec.z);

            // Initialize final color with ambient color
            Color finalColor = ambientColor;

            // Check if the material has a texture
            if (record.material->texture != nullptr) {
                // Apply texture to the diffuse color
                record.material->textureColor = applyTexture(record.hit_point, *(record.material));
            }

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

                    // No shadow, calculate diffuse and specular contributions
                    double diffuseFactor = std::max(0.0, record.normal.dot(lightDirection));
                    Vec3 diffuseColorVec = record.material->kd * record.material->diffusecolor * light->intensity * diffuseFactor;
                    Color temp(diffuseColorVec.x, diffuseColorVec.y, diffuseColorVec.z);
                    diffuseColor = temp;

                    // Specular reflection (Blinn-Phong model)
                    Vec3 viewDirection = -ray.getDirection().normalize();
                    Vec3 halfway = (lightDirection + viewDirection).normalize();
                    double specularFactor = std::pow(std::max(0.0, record.normal.dot(halfway)), record.material->specularexponent);
                    Vec3 specularColorVec = record.material->ks * record.material->specularcolor * light->intensity * specularFactor;
                    Color temp2(specularColorVec.x, specularColorVec.y, specularColorVec.z);
                    specularColor = temp2;

                    // Add diffuse and specular contributions to the final color
                    finalColor = finalColor + diffuseColor + specularColor;
                }
            }

            // Recursive reflection
            if (record.material->isreflective) {
                Vec3 reflectedDirection = Vec3::reflect(ray.getDirection().normalize(), record.normal);
                Ray reflectedRay(record.hit_point, reflectedDirection, 0.0);
                Color reflectedColor = phongRayColor(reflectedRay, objects, lights, depth + 1);
                finalColor = finalColor + record.material->reflectivity * reflectedColor;
            }

            // if (record.material->isreflective) {
            //     // Use fresnelReflection instead of reflect
            //     double reflectionCoefficient = Vec3::fresnelReflection(record.normal, ray.getDirection().normalize(), record.material->refractiveindex);

            //     Vec3 reflectedDirection = Vec3::reflect(ray.getDirection().normalize(), record.normal);
            //     Ray reflectedRay(record.hit_point, reflectedDirection, 0.0);
            //     Color reflectedColor = phongRayColor(reflectedRay, objects, lights, depth + 1);

            //     finalColor = finalColor + reflectionCoefficient * reflectedColor;
            // }

            // Recursive refraction
            if (record.material->isrefractive) {
                Vec3 refractedDirection;
                if (Vec3::refract(ray.getDirection().normalize(), record.normal, record.material->refractiveindex, refractedDirection)) {
                    Ray refractedRay(record.hit_point, refractedDirection, 0.0);
                    Color refractedColor = phongRayColor(refractedRay, objects, lights, depth + 1);
                    finalColor = finalColor + (1.0 - record.material->reflectivity) * refractedColor;
                }
            }

            // Gamma correction
            finalColor = applyGammaCorrection(finalColor, 0.8);
            return finalColor;
        }

        // If no intersection, return the background color
        return backgroundColor;
    }


};

#endif // PINHOLE_CAMERA_H
