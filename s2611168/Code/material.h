#ifndef CUSTOM_MATERIAL_H
#define CUSTOM_MATERIAL_H

#include "ray.h"
#include "hittable.h"
#include "vec3.h"
#include "texture.h"

#include <cmath>
#include "json_read.h" // Assuming you have the json_read.h file for parsing JSON

// forward declaration
// struct HitRecord;

class CustomMaterial{
public:
    CustomMaterial(const json& materialData, const int nbounces)
        : ks(materialData.at("ks").get<double>()),
          kd(materialData.at("kd").get<double>()),
          specularexponent(materialData.at("specularexponent").get<double>()),
          diffusecolor(Vec3(materialData.at("diffusecolor")[0], materialData.at("diffusecolor")[1], materialData.at("diffusecolor")[2])),
          specularcolor(Vec3(materialData.at("specularcolor")[0], materialData.at("specularcolor")[1], materialData.at("specularcolor")[2])),
          isreflective(materialData.at("isreflective").get<bool>()),
          reflectivity(materialData.at("reflectivity").get<double>()),
          isrefractive(materialData.at("isrefractive").get<bool>()),
          refractiveindex(materialData.at("refractiveindex").get<double>()),
          maxDepth(nbounces)
          {
            // if texture is present, then create a texture object
            // if (materialData.at("texture").is_null()) {
            //     texture = nullptr;
            // }
            // else {
            //     texture = new Texture(materialData.at("texture"));
            // }

            emittedColor = Vec3(0.0, 0.0, 0.0);
            ka = 0.2; // Default value
            // maxDepth = 1; // Default value
            
          }

    // another constructor for empty material. 
    CustomMaterial() {
        ks = 0.0;
        kd = 0.0;
        specularexponent = 0.0;
        diffusecolor = Vec3(0.0, 0.0, 0.0);
        specularcolor = Vec3(0.0, 0.0, 0.0);
        isreflective = false;
        reflectivity = 0.0;
        isrefractive = false;
        refractiveindex = 0.0;
        emittedColor = Vec3(0.0, 0.0, 0.0);
        maxDepth = 0;
    }

    // Function to compute the color of the material
    // This function is called from the shade function in the Ray class
    // Vec3 shade(const Ray& ray, const HitRecord& hitRecord) const {
    //     // Implementation of the shade function...
    //     return Vec3(0.0, 0.0, 0.0);  // Default color
    // }
    
    //
    // Implement any additional functions related to rendering, scattering, etc.
    // bool scatter(const Ray& ray, const HitRecord& record, Vec3& attenuation, Ray& scattered) const {
    //     // Lambertian Diffuse Reflection
    //     Vec3 target = record.hit_point + record.normal + random_in_unit_sphere();
    //     scattered = Ray(record.hit_point, target - record.hit_point, ray.getTime());  // Pass the time value
    //     attenuation = kd * diffusecolor / M_PI;

    //     // Specular Reflection
    //     if (isreflective && random_double() < ks) {
    //         Vec3 reflected = reflect(ray.getDirection().normalize(), record.normal);
    //         scattered = Ray(record.hit_point, reflected, ray.getTime());  // Pass the time value
    //         attenuation = attenuation + ks * specularcolor;
    //     }

    //     // Refraction (Transparency)
    //     if (isrefractive && random_double() < reflectivity) {
    //         Vec3 refracted;
    //         double reflect_prob = schlick(ray.getDirection(), record.normal, refractiveindex);
    //         if (random_double() >= reflect_prob && refract(ray.getDirection(), record.normal, refractiveindex, refracted)) {
    //             scattered = Ray(record.hit_point, refracted, ray.getTime());  // Pass the time value
    //             attenuation = (1.0 - reflectivity) * ks * specularcolor;
    //         }
    //     }

    //     // Adjust attenuation for Lambertian diffuse reflection
    //     attenuation = attenuation / M_PI;
    //     // Handling depth for recursive calls
    //     if (record.depth < maxDepth) {
    //         // Recursive calls for additional bounces
    //         // Example:
    //         // Vec3 additionalAttenuation;
    //         // Ray additionalScattered;
    //         // HitRecord additionalRecord = record;
    //         // additionalRecord.depth += 1;
    //         // if (additionalMaterial.scatter(scattered, additionalRecord, additionalAttenuation, additionalScattered)) {
    //         //     // Combine the result from additional bounces
    //         //     attenuation *= additionalAttenuation;
    //         //     scattered = additionalScattered;
    //         // }
    //     }

    //     return true;
    // }

    // bool scatter(const Ray& ray, const HitRecord& record, Vec3& attenuation, Ray& scattered) const {
    //     // Lambertian Diffuse Reflection
    //     Vec3 target = record.hit_point + record.normal + random_in_unit_sphere();
    //     scattered = Ray(record.hit_point, target - record.hit_point, ray.getTime());  
    //     attenuation = kd * diffusecolor / M_PI;

    //     // Specular Reflection
    //     if (isreflective && random_double() < ks) {
    //         Vec3 reflected = reflect(ray.getDirection().normalize(), record.normal);
    //         scattered = Ray(record.hit_point, reflected, ray.getTime());
    //         attenuation = attenuation + ks * specularcolor;
    //     }

    //     // Refraction (Transparency)
    //     if (isrefractive && random_double() < reflectivity) {
    //         Vec3 refracted;
    //         double reflect_prob = schlick(ray.getDirection(), record.normal, refractiveindex);
    //         if (random_double() >= reflect_prob && refract(ray.getDirection(), record.normal, refractiveindex, refracted)) {
    //             scattered = Ray(record.hit_point, refracted, ray.getTime());
    //             attenuation = (1.0 - reflectivity) * ks * specularcolor;
    //         }
    //     }

    //     return true;
    // }

    // make a static reflect function 
    static Vec3 reflect(const Vec3& incident, const Vec3& normal, const Vec3& hit_point) {
        return incident - 2 * incident.dot(normal) * normal;
    }

    Vec3 reflect(const Vec3& incident, const Vec3& normal) const {
        return incident - 2 * incident.dot(normal) * normal;
    }

    bool refract(const Vec3& incident, const Vec3& normal, double refractiveIndex, Vec3& refracted) const {
        Vec3 unit_incident = incident.normalize();
        double cos_theta = fmin(-unit_incident.dot(normal), 1.0);
        double sin_theta = std::sqrt(1.0 - cos_theta * cos_theta);

        Vec3 ref_direction = refractiveIndex * (unit_incident + cos_theta * normal) - sin_theta * normal;

        refracted = ref_direction.normalize();
        return true;
    }

    // make a static refract function
    static bool refract(const Vec3& incident, const Vec3& normal, const Vec3& hit_point, double refractiveIndex, Vec3& refracted) {
        Vec3 unit_incident = incident.normalize();
        double cos_theta = fmin(-unit_incident.dot(normal), 1.0);
        double sin_theta = std::sqrt(1.0 - cos_theta * cos_theta);

        Vec3 ref_direction = refractiveIndex * (unit_incident + cos_theta * normal) - sin_theta * normal;

        refracted = ref_direction.normalize();
        return true;
    }

    double schlick(const Vec3& incident, const Vec3& normal, double refractiveIndex) const {
        double cos_theta = fmin(-incident.dot(normal), 1.0);
        double r0 = (1.0 - refractiveIndex) / (1.0 + refractiveIndex);
        r0 = r0 * r0;

        return r0 + (1.0 - r0) * std::pow((1.0 - cos_theta), 5);
    }

    Vec3 random_in_unit_sphere() const {
        Vec3 p;
        do {
            p = 2.0 * Vec3(random_double(), random_double(), random_double()) - Vec3(1, 1, 1);
        } while (p.squared_length() >= 1.0);
        return p;
    }

    double random_double() const {
        // Generate a random number between 0 and 1
        return rand() / (RAND_MAX + 1.0);
    }

    // Function to compute Phong color
    // Vec3 phongColor(const Ray& ray, const HitRecord& hitRecord) const {
    //     // Ambient color (replace with actual ambient color)
    //     Vec3 ambientColor(0.1, 0.1, 0.1);

    //     // Diffuse color
    //     // Vec3 lightDirection = ;
    //     double diffuseIntensity = std::max(0.0, hitRecord.normal.dot(lightDirection));
    //     Vec3 diffuseColor = diffuseIntensity * diffusecolor;

    //     // Specular color
    //     // Vec3 viewDirection = ;
    //     Vec3 reflectedLight = reflect(-lightDirection, hitRecord.normal);
    //     double specularIntensity = std::pow(std::max(0.0, viewDirection.dot(reflectedLight)), specularexponent);
    //     Vec3 specularColor = specularIntensity * specularcolor;

    //     // Final Phong color
    //     Vec3 phongcolor = ambientColor + kd * diffuseColor + ks * specularColor;

    //     return phongcolor;
    // }

    Vec3 emitted() const {
        return emittedColor;
    }

    void print(){
        std::cout << "********CustomMaterial********" << std::endl;
        std::cout << "Ks: " << ks << std::endl;
        std::cout << "Kd: " << kd << std::endl;
        std::cout << "Ka: " << ka << std::endl;
        std::cout << "Specular Exponent: " << specularexponent << std::endl;
        std::cout << "Diffuse Color: " << diffusecolor << std::endl;
        std::cout << "Specular Color: " << specularcolor << std::endl;
        std::cout << "Is Reflective: " << isreflective << std::endl;
        std::cout << "Reflectivity: " << reflectivity << std::endl;
        std::cout << "Is Refractive: " << isrefractive << std::endl;
        std::cout << "Refractive Index: " << refractiveindex << std::endl;
        std::cout << "Emitted Color: " << emittedColor << std::endl;
        std::cout << "Max Depth: " << maxDepth << std::endl;
        // std::cout << "Texture: " << texture << std::endl;
        std::cout << "*****************************" << std::endl;

    }

// private:
    double ks;
    double kd;
    double ka; // ambient coefficient
    double specularexponent;
    Vec3 diffusecolor;
    Vec3 specularcolor;
    bool isreflective;
    double reflectivity;
    bool isrefractive;
    double refractiveindex;
    Vec3 emittedColor; 
    int maxDepth; // Maximum allowed depth for recursive ray tracing

    // Texture *texture;
};

// call the static reflect functoin 
// Vec3 CustomMaterial::reflect(const Vec3& incident, const Vec3& normal, const Vec3& hit_point) {
//     return incident - 2 * incident.dot(normal) * normal;
// }

#endif // CUSTOM_MATERIAL_H
