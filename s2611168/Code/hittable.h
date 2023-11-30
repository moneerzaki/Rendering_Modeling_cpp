#ifndef HITTABLE_H
#define HITTABLE_H

#include "ray.h"  // Assuming you have a Ray class
#include "vec3.h" // Assuming you have a Vec3 class
#include "color.h" // Assuming you have a Color class
#include "basics.h"
#include "lights.h"
#include "material2.h"
#include <vector>



// class to store hit information
class HitRecord {
public: 

    double t;           // Parameter along the ray
    Vec3 hit_point;     // Point of intersection
    Vec3 normal;        // Surface normal at the point of intersection
    Vec3 uv;  // UV coordinates of the hit point
    CustomMaterial* material;  // Pointer to the material of the object
    int depth; // Depth information

    // print all the information of the hit record
    void print() const {
        std::cout << "-----------HitRecord: ---------------" << std::endl;
        std::cout << "t: " << t << std::endl;
        std::cout << "hit_point: ";
        hit_point.print();
        std::cout << "normal: ";
        normal.print();
        std::cout << "depth: " << depth << std::endl;
        std::cout << "material: " << std::endl;
        // material->print();
    }
    // Constructor & overloaded constructor
    HitRecord() {
        t = 0.0;
        hit_point = Vec3(0.0, 0.0, 0.0);
        normal = Vec3(0.0, 0.0, 0.0);
        material = nullptr;
        depth = 0;
    }
    HitRecord(double t, Vec3 hit_point, Vec3 normal, CustomMaterial* material, int depth) {
        this->t = t;
        this->hit_point = hit_point;
        this->normal = normal;
        this->material = material;
        this->depth = depth;
    }
    // constructor has the uv edition
    HitRecord(double t, Vec3 hit_point, Vec3 normal, Vec3 uv, CustomMaterial* material, int depth) {
        this->t = t;
        this->hit_point = hit_point;
        this->normal = normal;
        this->uv = uv;
        this->material = material;
        this->depth = depth;
    }
    HitRecord(const HitRecord& record) {
        this->t = record.t;
        this->hit_point = record.hit_point;
        this->normal = record.normal;
        this->material = record.material;
        this->depth = record.depth;
    }
    HitRecord& operator=(const HitRecord& record) {
        this->t = record.t;
        this->hit_point = record.hit_point;
        this->normal = record.normal;
        this->material = record.material;
        this->depth = record.depth;
        return *this;
    }
    ~HitRecord() {
        // std::cout << "HitRecord destructor called" << std::endl;
    }


};

// class CustomMaterial;  // Forward declaration, assuming Material class will be defined later


// hittable class
class Hittable {
public:
    virtual ~Hittable() {}

    // get type 
    virtual std::string getType() const = 0;

    virtual bool hit(const Ray& ray, double t_min, double t_max, HitRecord& record) const = 0;

    virtual Color binaryColor() const =0;

    // virtual Vec3 phongColor(const Ray& ray, const HitRecord& hitRecord, const std::vector<Light*>& lights) const = 0;

    // Pure virtual function for scattering. Subclasses should override this.
    virtual bool scatter(const Ray& ray, const HitRecord& record, Vec3& attenuation, Ray& scattered) const = 0;

    // Pure virtual function for emitting light. Subclasses should override this.
    // virtual Vec3 emitted() const = 0;

    // Pure virtual function for printing the object. Subclasses should override this.
    virtual void print() const = 0;

};


#endif // HITTABLE_H
