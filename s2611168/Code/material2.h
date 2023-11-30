#ifndef CUSTOM_MATERIAL_H
#define CUSTOM_MATERIAL_H

#include "ray.h"
#include "hittable.h"
#include "vec3.h"
#include "texture.h"

#include <cmath>
#include "json_read.h" // Assuming you have the json_read.h file for parsing JSON
#include <nlohmann/json.hpp>

using json = nlohmann::json;

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
            texture = nullptr;
            // std::cout << materialData<< std::endl;

            try {
                // std::cout << "here" << std::endl;
                // std::string texture_value = materialData.at("texture");
                // std::cout << "Value of 'texture': " << texture_value << std::endl;
                if (materialData.contains("texture") && materialData["texture"].contains("type")) {
                    std::string textureType = materialData["texture"]["type"];
                    // Rest of your code based on textureType
                } else {
                    // Handle the case when either "texture" or "type" key is not found
                    std::cout << "Error: 'texture' or 'type' key not found in materialData" << std::endl;
                }
                if (materialData["texture"]["type"] == "woodwood") {
                    // std::cout << "not here" << std::endl;
                    texture = new WoodTexture(materialData["texture"]);
                    // texture->print(); 
                    // std::cout << "Wooden texture created" << std::endl;
                }
                else 
                {   
                    texture = nullptr;
                    // std::cout << "No texture created" << std::endl;
                }
            } catch (const nlohmann::json::out_of_range& e) {
                // Handle the case when the key is not found
                std::cout << "Error: " << e.what() << std::endl;
            } catch (const nlohmann::json::exception& e) {
                // Handle other JSON parsing errors
                std::cout << "JSON parsing error: " << e.what() << std::endl;
            } catch (const std::exception& e) {
                // Handle other general exceptions
                std::cout << "Error: " << e.what() << std::endl;
            }


            emittedColor = Vec3(0.0, 0.0, 0.0);
            ka = 0.05; // Default value
            // maxDepth = 1; // Default value
            
            textureColor = Vec3(0.0, 0.0, 0.0);
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

        texture = nullptr;
    }

    // function to edit all material data based on the new wooden texture 
    void edit_wooden_material(){}


    Vec3 reflect(const Vec3& incident, const Vec3& normal) const {
        return incident - 2 * incident.dot(normal) * normal;
    }

    bool refract(const Vec3& incident, const Vec3& normal, double refractiveIndex, Vec3& refracted) const {
        Vec3 unit_incident = incident.normalize();
        double cos_theta = fmin(-unit_incident.dot(normal), 1.0);
        double sin_theta = std::sqrt(1.0 - cos_theta * cos_theta);

        double critical_angle = asin(1.0 / refractiveIndex);

        if (sin_theta > critical_angle) {
            // Total internal reflection
            return false;
        }

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

    WoodTexture *texture;

    Vec3 textureColor; // color of the texture
};

#endif // CUSTOM_MATERIAL_H
