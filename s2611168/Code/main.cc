#include "json_read.h"
#include "camera5.h"
#include "lights.h"
#include "hittable.h"
#include "sphere.h"
#include "triangle.h"
#include "cylinder.h"
#include "material2.h"
#include "motion2.h"


#include <nlohmann/json.hpp>

using json = nlohmann::json;

// // Function to move and render the scene


int main(void) {

    // read and print the JSON file.
    try {
        std::string filePath1 = "../TestSuite/binary_primitves.json";  // Adjust the path as needed
        std::string filePath2 = "../TestSuite/mirror_image.json";  // Adjust the path as needed
        std::string filePath3 = "../TestSuite/scene_with_texture.json";  // Adjust the path as needed
        std::string filePath4 = "../TestSuite/scene.json";  // Adjust the path as needed
        std::string filePath5 = "../TestSuite/simple_phong_with_texture.json";  // Adjust the path as needed
        std::string filePath6 = "../TestSuite/simple_phong.json";  // Adjust the path as needed

        std::string binaryOutput = "./out";
        // std::string binaryOutput1 = "outputoutputoutput1";

        JsonData jsonData(filePath3);

        // Use the member functions to access and print the contents
        // jsonData.printData();

        

        // creating light sources 
        // checking if light sources are present in the scene
        std::vector<Light*> lights;

        if (jsonData.getScene()["lightsources"].is_null()) {
            std::cout << "lightsources not found" << std::endl;
        }
        else {
            // std::cout << "lightsources found" << std::endl;
            // for (const auto& lightData : jsonData.getScene()["lightsources"]) {
            for (int i = 0; i < jsonData.getScene()["lightsources"].size(); i++) {
                json lightData = jsonData.getScene()["lightsources"][i];
                // std::cout << lightData << std::endl;
                // push_back() creates a copy of the object
                lights.push_back(new Light(lightData));
            }
            // cout << "lightsources size: " << lights.size() << endl;
        }
        // print_lightsources(lights);



        // Creating the hittable objects vector
        std::vector<Hittable*> hittableObjects;
        // Creating the shapes
        json shapes = jsonData.getScene().at("shapes");
        for (int i = 0; i < shapes.size(); i++) {

            // Creating the material of the shape first
            CustomMaterial* customMaterial = nullptr;
            if (shapes[i]["material"].is_null()) {
                std::cout << "material not found" << std::endl;
            }
            else {
                // std::cout << "material found" << std::endl;
                json materialData = shapes[i]["material"];
                customMaterial = new CustomMaterial(materialData, jsonData.getNbounces());
            }
            
            // Creating the shape based on the type
            json shapeData = shapes[i];
            Hittable* shape = nullptr;

            if (shapeData.at("type") == "sphere") {
                shape = new Sphere(shapeData, customMaterial);
            } 
            else if (shapeData.at("type") == "triangle") {
                // Assuming you have a Triangle class
                shape = new Triangle(shapeData, customMaterial);
            } 
            else if (shapeData.at("type") == "cylinder") {
                // Assuming you have a Cylinder class
                shape = new Cylinder(shapeData, customMaterial);
            } 
            else {
                std::cout << "Error: Unknown shape type." << std::endl;
                // return 1;
            }
            // Add the shape to the vector of hittable objects
            if (shape != nullptr) {
                hittableObjects.push_back(shape);
            }
        }


        // print_hittableObjects(hittableObjects);
        // for (int i = 0; i < hittableObjects.size(); i++) {
        //     // cout << "Hittable object " << i << endl;
        //     hittableObjects[i]->print();
        // }


        // rendering the scene 
        PinholeCamera camera(jsonData.getCamera(), jsonData.getScene().at("backgroundcolor"), jsonData.getNbounces());
        // camera.print();
        
            // camera.binaryRendering(hittableObjects, binaryOutput);
            // camera.binaryRenderingMT(hittableObjects, binaryOutput, 2);

            // auto start = std::chrono::high_resolution_clock::now();
            camera.phongRendering(hittableObjects, lights, binaryOutput);
            // auto stop = std::chrono::high_resolution_clock::now();
            // auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
            // std::cout << std::endl << "Time taken without MT: " << duration.count() << " microseconds" << std::endl;

            // start = std::chrono::high_resolution_clock::now();
            // camera.phongRenderingMultithreaded(hittableObjects, lights, binaryOutput);
            // stop = std::chrono::high_resolution_clock::now();
            // duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
            // std::cout << std::endl << "Time taken with MT: " << duration.count() << " microseconds" << std::endl;

        // Function to change positions of hittable objects and render the scene
        // MoveAndRender(camera,hittableObjects, lights,binaryOutput);


    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    return 0;
}
