// lights.h
#ifndef LIGHTS_H
#define LIGHTS_H

#include <iostream>
#include <nlohmann/json.hpp>
#include "vec3.h"

using json = nlohmann::json;

class Light {
public:
    // Constructor that takes a JSON object and initializes the light
    Light(const json& lightData) :
        type(lightData.at("type").get<std::string>()),
        position(Vec3(lightData.at("position")[0], lightData.at("position")[1], lightData.at("position")[2])),
        intensity(Vec3(lightData.at("intensity")[0], lightData.at("intensity")[1], lightData.at("intensity")[2]))
        // Initialize other light parameters...
    {
        // Other initialization...
    }

    // Accessor methods to retrieve light properties
    std::string getType() const {
        return type;
    }

    Vec3 getPosition() const {
        return position;
    }

    Vec3 getIntensity() const {
        return intensity;
    }

    void print() const {
        std::cout << "Light type: " << type << std::endl;
        std::cout << "Light position: ";
        position.print();
        std::cout << "Light intensity: ";
        intensity.print();
        // Print other light properties...
    }

// private:
    // Private member variables
    std::string type;
    Vec3 position;
    Vec3 intensity;
    // Add more private fields as needed
};

#endif // LIGHTS_H
