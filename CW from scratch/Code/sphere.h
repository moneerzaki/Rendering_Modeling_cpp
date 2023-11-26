#ifndef SPHERE_H
#define SPHERE_H

#include "hittable.h"
// #include "hitrecord.h"
#include "ray.h"       // Assuming you have a Ray class
#include "vec3.h"      // Assuming you have a Vec3 class

// #include "material.h"  // Assuming you have a Material class
#include "json_read.h" // Assuming you have the json_read.h file for parsing JSON

class Sphere : public Hittable {
public:
    Sphere(const json& sphereData, CustomMaterial* material)
        : center(Vec3(sphereData.at("center")[0], sphereData.at("center")[1], sphereData.at("center")[2])),
          radius(sphereData.at("radius").get<double>()),
          material(material) {}

    virtual std::string getType() const override {
        return "sphere";
    }

    // this is a very good working function for the sphere ....
    bool hit(const Ray& ray, double t_min, double t_max, HitRecord& record) const override {
        Vec3 oc = ray.getOrigin() - center;
        double a = ray.getDirection().length_squared();
        // double half_b = dot(oc, ray.getDirection());
        double half_b = oc.dot(ray.getDirection());
        double c = oc.length_squared() - radius * radius;

        double discriminant = half_b * half_b - a * c;

        if (discriminant > 0) {
            double root = (-half_b - sqrt(discriminant)) / a;
            if (root < t_max && root > t_min) {
                record.t = root;
                record.hit_point = ray.at(root);
                record.normal = (record.hit_point - center) / radius;
                record.material = material;
                return true;
            }

            root = (-half_b + sqrt(discriminant)) / a;
            if (root < t_max && root > t_min) {
                record.t = root;
                record.hit_point = ray.at(root);
                record.normal = (record.hit_point - center) / radius;
                record.material = material;
                return true;
            }
        }

        return false;
    }

    // hit version with the uv edit 
    // bool hit(const Ray& ray, double t_min, double t_max, HitRecord& record) const override {
    //     Vec3 oc = ray.getOrigin() - center;
    //     double a = ray.getDirection().length_squared();
    //     double half_b = oc.dot(ray.getDirection());
    //     double c = oc.length_squared() - radius * radius;

    //     double discriminant = half_b * half_b - a * c;

    //     if (discriminant > 0) {
    //         double root = (-half_b - sqrt(discriminant)) / a;
    //         if (root < t_max && root > t_min) {
    //             record.t = root;
    //             record.hit_point = ray.at(root);
    //             record.normal = (record.hit_point - center) / radius;

    //             // Calculate UV coordinates for the wooden texture
    //             Vec3 sphereToHit = record.hit_point - center;
    //             double phi = atan2(sphereToHit.z, sphereToHit.x);
    //             double theta = asin(sphereToHit.y / radius);
    //             double u = 1.0 - (phi + M_PI) / (2.0 * M_PI);
    //             double v = (theta + 0.5 * M_PI) / M_PI;

    //             record.uv = Vec3(u, v, 0.0);  // Assuming Vec3 is used for UV coordinates
    //             record.material = material;
    //             return true;
    //         }

    //         root = (-half_b + sqrt(discriminant)) / a;
    //         if (root < t_max && root > t_min) {
    //             record.t = root;
    //             record.hit_point = ray.at(root);
    //             record.normal = (record.hit_point - center) / radius;

    //             // Calculate UV coordinates for the wooden texture
    //             Vec3 sphereToHit = record.hit_point - center;
    //             double phi = atan2(sphereToHit.z, sphereToHit.x);
    //             double theta = asin(sphereToHit.y / radius);
    //             double u = 1.0 - (phi + M_PI) / (2.0 * M_PI);
    //             double v = (theta + 0.5 * M_PI) / M_PI;

    //             record.uv = Vec3(u, v, 0.0);  // Assuming Vec3 is used for UV coordinates
    //             record.material = material;
    //             return true;
    //         }
    //     }
    //     return false;
    // }


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
        return Color(1.0, 0.0, 0.0);  // Blue color
    }

    // Vec3 phongColor(const Ray& ray, const HitRecord& hitRecord, const std::vector<Light*>& lights, const std::vector<Hittable*>& objects) const override {
    //     Vec3 ambientColor(0.1, 0.1, 0.1); // Ambient color

    //     Vec3 totalColor = ambientColor * hitRecord.material->diffusecolor; // Ambient component

    //     for (const auto& light : lights) {
    //         Vec3 lightDirection = (light->position - hitRecord.hit_point).normalize();

    //         // Diffuse component
    //         double diffuseIntensity = std::max(0.0, hitRecord.normal.dot(lightDirection));
    //         Vec3 diffuseColor = diffuseIntensity * hitRecord.material->diffusecolor;

    //         // Specular component
    //         Vec3 viewDirection = (ray.getOrigin() - hitRecord.hit_point).normalize();
    //         Vec3 reflectedLight = hitRecord.material->reflect(-lightDirection, hitRecord.normal);
    //         double specularIntensity = std::pow(std::max(0.0, viewDirection.dot(reflectedLight)), hitRecord.material->specularexponent);
    //         Vec3 specularColor = specularIntensity * hitRecord.material->specularcolor;

    //         // Shadows: Check if there's an obstruction between the hit point and the light source
    //         Ray shadowRay(hitRecord.hit_point, lightDirection, 0.0);
    //         HitRecord shadowHit;
    //         bool inShadow = false;

    //         // 
    //         for (const auto& object : objects) {
    //             if (object->hit(shadowRay, 0.001, std::numeric_limits<double>::infinity(), shadowHit)) {
    //                 inShadow = true;
    //                 break;
    //             }
    //         }

    //         if (!inShadow) {
    //             totalColor = totalColor + diffuseColor + specularColor;
    //         }
    //     }

    //     return totalColor;
    // }



    void print() const override{
        std::cout << "********* Sphere: *********" << std::endl;
        std::cout << "Center: " << center.x << " " << center.y << " " << center.z << std::endl;
        std::cout << "Radius: " << radius << std::endl;
        if (material == nullptr) {
            std::cout << "Material is null" << std::endl;
        } else {
            // std::cout << "Material is not null" << std::endl;
            material->print();
        }
    }

    void setPosition(const Vec3& position) {
        center = position;
    }

private:
    Vec3 center;
    double radius;
    CustomMaterial* material;
};

#endif // SPHERE_H
