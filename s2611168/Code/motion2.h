#ifndef MOTION_H
#define MOTION_H

#include "hittable.h"
#include "camera5.h"
#include "lights.h"
#include "material2.h"
#include <vector>
#include <string>

void MoveAndRender(PinholeCamera& camera,
                    std::vector<Hittable*>& hittableObjects, 
                    const std::vector<Light*>& lights,
                    const std::string& binaryOutput) {
    const int numFrames = 300; // Adjust as needed
    const double animationDuration = 10.0; // Adjust as needed

    // Small changes in light intensity and position
    Vec3 deltaIntensity(0.1, 0.1, 0.1);
    Vec3 deltaPosition(0.01, 0, 0);
    

    // Calculate the step for each frame
    double step = animationDuration / numFrames;

    for (int frame = 0; frame < numFrames; ++frame) {
        double t = frame * step;

        // light sources 
        // Update light position and intensity
        for (auto& light : lights) {
            light->position = light->position +  deltaPosition;
            // light->intensity = light->intensity + deltaIntensity ;
        }

        // Camera Simple diagonal motion 
        double cameraX = 0.5 * std::sin(t); // Adjust the factor to control the range of motion
        double cameraZ = 0.5 * std::cos(t); // Adjust the factor to control the range of motion
        // Update camera position
        camera.position = Vec3(cameraX, 0.75, -cameraZ);

        // Example: Rotate the cylinder around its center to keep its position fixed
        double rotationAngle = 0.05 * t; // Adjust the factor to control the speed of rotation
        if (hittableObjects[1]->getType() == "cylinder") 
        {dynamic_cast<Cylinder*>(hittableObjects[1])->rotate(Vec3(0, 1, 0), rotationAngle);}

        // Example: Move the cylinder in a circular motion within the camera range
        // double cylinderX = 0.2 * std::cos(0.05 * t); // Adjust the factors to control the speed and range of motion
        // double cylinderZ = 0.2 * std::sin(0.05 * t); // Adjust the factors to control the speed and range of motion
        // if (hittableObjects[1]->getType() == "cylinder") {
        //     dynamic_cast<Cylinder*>(hittableObjects[1])->move(Vec3(cylinderX, 0.19, cylinderZ));
        // }

        // Example: Move the first sphere in a slow vertical motion with reversal
        double sphere1Y = 0.1 * std::sin(0.02 * t); // Adjust the factors to control the speed and range of motion
        // Revert the direction when t crosses a certain point (e.g., halfway through the animation)
        if (t > animationDuration / 2.0) 
        {sphere1Y *= -1.0;}
        if (hittableObjects[4]->getType() == "sphere") 
        {dynamic_cast<Sphere*>(hittableObjects[4])->move(Vec3(0, sphere1Y, 0));}
        // Example: Make the sphere jump up and down to rest at the top of the cylinder height
        double sphere2Y = 0.1 * std::sin(0.02 * t); // Adjust the factors to control the speed and range of motion
        // Revert the direction when t crosses a certain point (e.g., halfway through the animation)
        // if (frame < numFrames-numFrames/4){
        if (t > animationDuration * 3 / 4) 
        {sphere2Y *= -1.0;}
        if (hittableObjects[5]->getType() == "sphere") 
        {dynamic_cast<Sphere*>(hittableObjects[5])->move(Vec3(0, sphere2Y, 0));}
        // }
        // make a mroe interesting motion to the sphere
        


        // Example: Move the falling sphere
        // if (hittableObjects[6]->getType() == "sphere") {
        //     double fallingSphereY = std::max(0.0, 0.5 - 0.1 * t);
        //     dynamic_cast<Sphere*>(hittableObjects[6])->setPosition(Vec3(0.3, fallingSphereY, 1));
        // }

        // Example: Apply random motion to the two triangles
        if (hittableObjects[2]->getType() == "triangle") {
            dynamic_cast<Triangle*>(hittableObjects[2])->moveRandomly(0.02);
        }
        if (hittableObjects[3]->getType() == "triangle") {
            dynamic_cast<Triangle*>(hittableObjects[3])->moveRandomly(0.02);
        }

        // Render the scene
        auto start = std::chrono::steady_clock::now();
        std::string final = binaryOutput + std::to_string(frame);
        // camera.binaryRendering(hittableObjects, final);
        // phong rendering 
        camera.phongRendering(hittableObjects, lights, final);
        auto end = std::chrono::steady_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

        std::cout << "Frame " << frame + 1 << "/" << numFrames << " rendered in " << duration << " ms." << std::endl;
    }
}

#endif // MOTION_H
