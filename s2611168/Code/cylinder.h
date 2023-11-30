// cylinder.h
#ifndef CYLINDER_H
#define CYLINDER_H

#include "hittable.h"
// #include "hitrecord.h"
// #include "material.h" // Assuming you have a CustomMaterial class
#include <cmath>


class Cylinder : public Hittable {
public:
    Cylinder(const json& cylinderData, CustomMaterial* mat) : material(mat) {
        center = Vec3(cylinderData.at("center")[0], cylinderData.at("center")[1], cylinderData.at("center")[2]);
        axis = Vec3(cylinderData.at("axis")[0], cylinderData.at("axis")[1], cylinderData.at("axis")[2]).normalize();
        radius = cylinderData.at("radius").get<double>();
        height = cylinderData.at("height").get<double>();

        // double the height and modify the center to the midpoint of the height. 
        height = 2.0 * height;
        center = center - axis * height / 2.0;
    }

    // get type
    virtual std::string getType() const override {
        return "cylinder";
    }

    ~Cylinder() {}

    // ahit function that works for the a hollow cylinder
    // a very good working function for the cylinder .... 
    bool hit(const Ray& ray, double t_min, double t_max, HitRecord& record) const override {
        
        // Shift the ray to the local coordinate system of the cylinder
        Vec3 oc = ray.getOrigin() - center;
        Vec3 oc_cross_axis = cross(oc, axis);
        Vec3 dir_cross_axis = cross(ray.getDirection(), axis);

        double a = dir_cross_axis.length_squared();
        double b = 2.0 * oc_cross_axis.dot(dir_cross_axis);
        double c = oc_cross_axis.length_squared() - radius * radius;

        double discriminant = b * b - 4 * a * c;

        if (discriminant > 0) {
            double root1 = (-b - sqrt(discriminant)) / (2.0 * a);
            double root2 = (-b + sqrt(discriminant)) / (2.0 * a);

            for (double root : {root1, root2}) {
                if (root < t_max && root > t_min) {
                    Vec3 hit_point = ray.at(root);
                    // double projection = dot(hit_point - center, axis);
                    double projection = (hit_point - center).dot( axis);

                    if (projection > 0 && projection < height) {
                        // Valid intersection inside the hollow cylinder
                        record.t = root;
                        record.hit_point = hit_point;
                        record.normal = (oc + root * dir_cross_axis / a - projection * axis / height).normalize();
                        record.material = material;
                        return true;
                    }
                }
            }
        }

        return false;
    }


    Vec3 computeNormal(const Vec3& point) const {
        // Calculate the normal vector at the given point on the cylinder
        Vec3 normal(point.x - center.x, 0, point.z - center.z);
        return normal.normalize();
    }

    Color binaryColor() const override {
        // Return green color for the cylinder
        return Color(0.0, 1.0, 0.0);
    }


    bool scatter(const Ray& ray, const HitRecord& record, Vec3& attenuation, Ray& scattered) const override {
        // Lambertian Diffuse Reflection
        Vec3 target = record.hit_point + record.normal + random_in_unit_sphere();
        scattered = Ray(record.hit_point, target - record.hit_point, ray.getTime());
        attenuation = record.material->kd * record.material->diffusecolor / M_PI;


        // // Specular Reflection
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



    void print() const override {
        std::cout << "Cylinder Attributes:" << std::endl;
        std::cout << "Center: " << center << std::endl;
        std::cout << "Axis: " << axis << std::endl;
        std::cout << "Radius: " << radius << std::endl;
        std::cout << "Height: " << height << std::endl;
        std::cout << "Material: ";
        material->print();
        std::cout << std::endl;
    }

    void move(const Vec3& displacement) {
        center = center + displacement;
    }

    void rotate (const Vec3& axis, double angle) {
        center = center.rotate(angle, axis);
    }

    
    
private:
    Vec3 center;
    Vec3 axis;
    double radius;
    double height;
    CustomMaterial* material;
};

#endif // CYLINDER_H
