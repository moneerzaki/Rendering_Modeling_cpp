#ifndef TRIANGLE_H
#define TRIANGLE_H



#include "hittable.h"
// #include "hitrecord.h"

// #include "material.h"  // Assuming you have a Material class
#include "json_read.h" // Assuming you have the json_read.h file for parsing JSON

class Triangle : public Hittable {
public:
    Triangle(const json& triangleData, CustomMaterial* material)
        : v0(Vec3(triangleData.at("v0")[0], triangleData.at("v0")[1], triangleData.at("v0")[2])),
          v1(Vec3(triangleData.at("v1")[0], triangleData.at("v1")[1], triangleData.at("v1")[2])),
          v2(Vec3(triangleData.at("v2")[0], triangleData.at("v2")[1], triangleData.at("v2")[2])),
          material(material) {}


    // working for one triangle only.
    virtual bool hit(const Ray& ray, double t_min, double t_max, HitRecord& record) const override {
        Vec3 e1 = v1 - v0;
        Vec3 e2 = v2 - v0;
        Vec3 h = ray.getDirection().cross(e2);
        double a = e1.dot(h);

        if (a > -EPSILON && a < EPSILON) {
            return false;  // Ray is parallel to the triangle
        }

        double f = 1.0 / a;
        Vec3 s = ray.getOrigin() - v0;
        double u = f * s.dot(h);

        if (u < 0.0 || u > 1.0) {
            return false;
        }

        Vec3 q = s.cross(e1);
        double v = f * ray.getDirection().dot(q);

        if (v < 0.0 || u + v > 1.0) {
            return false;
        }

        double t = f * e2.dot(q);

        if (t > t_min && t < t_max) {
            record.t = t;
            record.hit_point = ray.at(t);
            record.normal = e1.cross(e2).normalize();
            record.material = material;
            return true;
        }

        return false;
    }

    // get type 
    virtual std::string getType() const override {
        return "triangle";
    }

    bool scatter(const Ray& ray, const HitRecord& record, Vec3& attenuation, Ray& scattered) const override {
        // Lambertian Diffuse Reflection
        Vec3 target = record.hit_point + record.normal + random_in_unit_sphere();
        scattered = Ray(record.hit_point, target - record.hit_point, ray.getTime());
        attenuation = record.material->kd * record.material->diffusecolor / M_PI;

        // Specular Reflection
        if (record.material->isreflective && random_double() < record.material->ks) {
            Vec3 reflected = record.material->reflect(ray.getDirection().normalize(), record.normal);
            scattered = Ray(record.hit_point, reflected, ray.getTime());
            attenuation = attenuation + record.material->ks * record.material->specularcolor;
        }

        // Refraction (Transparency)
        if (record.material->isrefractive && random_double() < record.material->reflectivity) {
            Vec3 refracted;
            double reflect_prob = record.material->schlick(ray.getDirection(), record.normal, record.material->refractiveindex);
            if (random_double() >= reflect_prob && record.material->refract(ray.getDirection(), record.normal, record.material->refractiveindex, refracted)) {
                scattered = Ray(record.hit_point, refracted, ray.getTime());
                attenuation = (1.0 - record.material->reflectivity) * record.material->ks * record.material->specularcolor;
            }
        }

        // Adjust attenuation for Lambertian diffuse reflection
        attenuation = attenuation / M_PI;

        return true;
    }
    
    Color binaryColor() const {
        // Return the specific color for the sphere in binary rendering
        return Color(0.0, 0.0, 1.0);  // Blue color
    }

    // Vec3 phongColor(const Ray& ray, const HitRecord& hitRecord, const std::vector<Light*>& lights, const std::vector<Hittable*>& objects) const override{
    //     Vec3 ambientColor(0.1, 0.1, 0.1);  // Ambient color
    //     Vec3 diffuseColor(0.8, 0.8, 0.8);  // Diffuse color
    //     Vec3 specularColor(1.0, 1.0, 1.0); // Specular color

    //     Vec3 totalColor = ambientColor * hitRecord.material->diffusecolor;

    //     for (const auto& light : lights) {
    //         // Direction from hit point to light source
    //         Vec3 lightDirection = (light->getPosition() - hitRecord.hit_point).normalize();

    //         // Lambertian Diffuse Reflection
    //         double diffuseIntensity = std::max(0.0, hitRecord.normal.dot(lightDirection));
    //         Vec3 diffuseContrib = diffuseIntensity * diffuseColor * hitRecord.material->diffusecolor;

    //         // Specular Reflection
    //         Vec3 viewDirection = (ray.getOrigin() - hitRecord.hit_point).normalize();
    //         Vec3 reflectedLight = hitRecord.material->reflect(-lightDirection, hitRecord.normal);
    //         double specularIntensity = std::pow(std::max(0.0, viewDirection.dot(reflectedLight)), hitRecord.material->specularexponent);
    //         Vec3 specularContrib = specularIntensity * specularColor * hitRecord.material->specularcolor;

    //         // Combine ambient, diffuse, and specular contributions
    //         totalColor = totalColor + light->getIntensity() * (diffuseContrib + specularContrib);
    //     }

    //     return totalColor;
    // }


    void print() const override{
        std::cout << "******Triangle*******: " << std::endl;
        std::cout << "v0: " << v0.x << " " << v0.y << " " << v0.z << std::endl;
        std::cout << "v1: " << v1.x << " " << v1.y << " " << v1.z << std::endl;
        std::cout << "v2: " << v2.x << " " << v2.y << " " << v2.z << std::endl;
        // std::cout << "Material: " << std::endl;
        if (material == nullptr) {
            std::cout << "Material is null" << std::endl;
        } else {
            // std::cout << "Material is not null" << std::endl;
            material->print();
        }
        
    }

private:
    Vec3 v0, v1, v2;
    CustomMaterial* material;
};

#endif // TRIANGLE_H
