// json_read.h
#ifndef JSON_READ_H
#define JSON_READ_H

#include <iostream>
#include <fstream>
#include <nlohmann/json.hpp>

using json = nlohmann::json;

class JsonData {
public:
    // Constructor that reads the JSON file
    JsonData(const std::string& filePath) {
        std::ifstream file(filePath);

        if (!file.is_open()) {
            throw std::runtime_error("Failed to open the JSON file");
        }
        
        json jsonData;
        file >> jsonData;

        try {
            // Check if 'nbounces' key exists
            // if (jsonData.find("nbounces") == jsonData.end()) {
            //     std::cout << "nbounces not found, using default value." << std::endl;
            //     nbounces = 8;  // Set a default value or choose another approach
            // } else {
            //     nbounces = jsonData["nbounces"];
            // }
            nbounces = 8; // Set a default value or choose another approach

            rendermode = jsonData.at("rendermode").get<std::string>();
            camera = jsonData.at("camera");
            scene = jsonData.at("scene");

        } catch (const json::exception& e) {
            throw std::runtime_error("Error parsing JSON: " + std::string(e.what()));
        }
    }

    // Public member functions

    // Getters
    int getNbounces() const {
        return nbounces;
    }
    std::string getRenderMode() const {
        return rendermode;
    }
    json getCamera() const {
        return camera;
    }

    json getScene() const {
        return scene;
    }
    // Function to print the data
    void printData() const {
        std::cout << "Nbounces: " << nbounces << std::endl;
        std::cout << "Render Mode: " << rendermode << std::endl;

        // Print all camera data
        std::cout << "Camera Type: " << camera.at("type").get<std::string>() << std::endl;
        std::cout << "Camera Width: " << camera.at("width") << std::endl;
        std::cout << "Camera Height: " << camera.at("height") << std::endl;
        std::cout << "Camera Position: " << camera.at("position") << std::endl;
        std::cout << "Camera LookAt: " << camera.at("lookAt") << std::endl;
        std::cout << "Camera UpVector: " << camera.at("upVector") << std::endl;
        std::cout << "Camera Fov: " << camera.at("fov") << std::endl;
        std::cout << "Camera Exposure: " << camera.at("exposure") << std::endl;

        // Print all scene data
        std::cout << "Scene Background Color: " << scene.at("backgroundcolor") << std::endl;
        std::cout << "Scene Lights: " << std::endl;
        if (scene.find("lightsources") == scene.end()) {
            std::cout << "lightsources not found" << std::endl;
        } else {
            for (const auto& light : scene.at("lightsources")) {
                std::cout << "Light Type: " << light.at("type").get<std::string>() << std::endl;
                std::cout << "Light Position: " << light.at("position") << std::endl;
                std::cout << "Light Intensity: " << light.at("intensity") << std::endl;
            }
        }
        
        std::cout << "Scene shapes: " << std::endl;
        for (const auto& shape : scene.at("shapes")) {
            std::cout << "Shape Type: " << shape.at("type").get<std::string>() << std::endl;

            if (shape.at("type").get<std::string>() == "sphere") {
                std::cout << "Sphere Center: " << shape.at("center") << std::endl;
                std::cout << "Sphere Radius: " << shape.at("radius") << std::endl;
            } else if (shape.at("type").get<std::string>() == "triangle") {
                std::cout << "Triangle Vertex 1: " << shape.at("v0") << std::endl;
                std::cout << "Triangle Vertex 2: " << shape.at("v1") << std::endl;
                std::cout << "Triangle Vertex 3: " << shape.at("v2") << std::endl;
            } else if (shape.at("type").get<std::string>() == "cylinder") {
                std::cout << "Cylinder Center: " << shape.at("center") << std::endl;
                std::cout << "Cylinder Axis  : " << shape.at("axis"  ) << std::endl;
                std::cout << "Cylinder Radius: " << shape.at("radius") << std::endl;
                std::cout << "Cylinder Height: " << shape.at("height") << std::endl;
            }

            if (shape.find("material") == shape.end()) {
                std::cout << "material not found" << std::endl;
            } else {
                // std::cout << "Shape Material: " << shape.at("material") << std::endl;
            }
            // std::cout << "Shape Transformations: " << std::endl;
            // for (const auto& transformation : shape.at("transformations")) {
            //     std::cout << "Transformation Type: " << transformation.at("type").get<std::string>() << std::endl;
            //     std::cout << "Transformation Parameters: " << transformation.at("parameters") << std::endl;
            // }
        }
        
        // Print other scene fields as needed
    }

private:
    // Private member variables
    int nbounces;
    std::string rendermode;
    json camera;
    json scene;

    // Add more private fields as needed
};

#endif // JSON_READ_H
