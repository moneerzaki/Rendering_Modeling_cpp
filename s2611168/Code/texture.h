#ifndef TEXTURE_H
#define TEXTURE_H

#include "vec3.h"
#include <string>
#include <iostream> // for debugging purposes
#include <nlohmann/json.hpp>

using json = nlohmann::json;

class Texture {
public:
    virtual Vec3 getColor(double u, double v) const = 0;
};

class SolidColorTexture : public Texture {
public:
    SolidColorTexture(const Vec3& color) : color(color) {}

    virtual Vec3 getColor(double u, double v) const override {
        // For a solid color texture, ignore UV coordinates and return the constant color
        return color;
    }

private:
    Vec3 color;
};

class ImageTexture : public Texture {
public:
    ImageTexture(const std::string& filename) {
        // Load image and initialize data structures
        // Implementation details depend on the image loading library you are using
        // Example: stb_image, OpenCV, etc.
        // ...

        std::cout << "Loaded image from file: " << filename << std::endl;
    }

    virtual Vec3 getColor(double u, double v) const override {
        // Map UV coordinates to pixel coordinates in the image
        int i = static_cast<int>(u * width);
        int j = static_cast<int>((1 - v) * height);

        // Retrieve color from the image at the specified pixel
        // Implementation details depend on the image loading library
        // Example: return image[i][j];
        // ...

        return Vec3(1.0, 1.0, 1.0); // Placeholder color for now
    }

// private:
    int width;  // Image width
    int height; // Image height
    // Additional data structures to store image data
    // ...
};


// Assume WoodTexture inherits from the Texture class
class WoodTexture : public Texture {
public:
    // Constructor using json data 
    WoodTexture(json j) {
        scale = j["scale"];
        // std::cout << "color variations: " << j["colorVariations"] << std::endl;
        colorMin = Vec3(j["colorVariations"]["min"][0], j["colorVariations"]["min"][1], j["colorVariations"]["min"][2]);
        colorMax = Vec3(j["colorVariations"]["max"][0], j["colorVariations"]["max"][1], j["colorVariations"]["max"][2]);
        // // colorMin = Vec3(j["colorMin"][0], j["colorMin"][1], j["colorMin"][2]);
        // // colorMax = Vec3(j["colorMax"][0], j["colorMax"][1], j["colorMax"][2]);
        ringSpacing = j["ringSpacing"];
        ringThickness = j["ringThickness"];
        noise = j["noise"];

        // print();
    }
    WoodTexture(double scale, const Vec3& colorMin, const Vec3& colorMax,
                double ringSpacing, double ringThickness, double noise)
        : scale(scale), colorMin(colorMin), colorMax(colorMax),
          ringSpacing(ringSpacing), ringThickness(ringThickness), noise(noise) {}


    // getcolor function depending on u and v
    virtual Vec3 getColor(double u, double v) const override {
        // Apply the wood texture algorithm based on UV coordinates
        double woodPattern = std::sin(u * scale * M_PI) * std::sin(v * scale * M_PI);
        woodPattern = woodPattern * 0.5 + 0.5;  // Map [-1, 1] to [0, 1]

        // Add ring pattern
        double ring = std::fmod(woodPattern + ringSpacing, 1.0);
        ring = std::clamp((ring - ringThickness) / (ringThickness * 0.5), 0.0, 1.0);

        // Apply color variations
        Vec3 color = mix(colorMin, colorMax, ring);

        // Add noise
        double noiseValue = (std::rand() % 1000 / 1000.0) * 2.0 - 1.0;
        color = mix(color, Vec3(1.0, 1.0, 1.0), noiseValue * noise);

        return color;
    }

    // print function for all attributes. 
    void print() {
        std::cout << "scale: " << scale << std::endl;
        std::cout << "colorMin: " << colorMin << std::endl;
        std::cout << "colorMax: " << colorMax << std::endl;
        std::cout << "ringSpacing: " << ringSpacing << std::endl;
        std::cout << "ringThickness: " << ringThickness << std::endl;
        std::cout << "noise: " << noise << std::endl;
    }

    // getcolor function depending on u and v 
    // Vec3 getColor(double u, double v) {
    //     // Scale the UV coordinates
    //     u *= scale;
    //     v *= scale;

    //     // Calculate wood texture based on parameters
    //     double rings = std::sin(u / ringSpacing) + std::sin(v / ringSpacing);
    //     rings = std::abs(rings) + noise * random_double(-1, 1);

    //     // Apply ring thickness
    //     double woodColor = smoothstep(0.5 - ringThickness, 0.5, rings);

    //     // Interpolate between colorMin and colorMax based on woodColor
    //     Vec3 finalColor = lerp(colorMin, colorMax, woodColor);

    //     return finalColor;
    // }

// private:
    double scale;
    Vec3 colorMin;
    Vec3 colorMax;
    double ringSpacing;
    double ringThickness;
    double noise;

    // Utility function for linear interpolation between two values
    Vec3 mix(const Vec3& a, const Vec3& b, double t) const {
        return a * (1.0 - t) + b * t;
    }

    // Utility function for linear interpolation
    static double lerp(double a, double b, double t) {
        return a + t * (b - a);
    }

    // Utility function for linear interpolation between Vec3
    static Vec3 lerp(const Vec3& a, const Vec3& b, double t) {
        return Vec3(lerp(a.x, b.x, t), lerp(a.y, b.y, t), lerp(a.z, b.z, t));
    }

    // Utility function for smoothstep interpolation
    static double smoothstep(double edge0, double edge1, double x) {
        x = std::clamp((x - edge0) / (edge1 - edge0), 0.0, 1.0);
        return x * x * (3 - 2 * x);
    }
};



class BrickTexture : public Texture {
public:
    BrickTexture(double mortarWidth) : mortarWidth(mortarWidth) {}

    virtual Vec3 getColor(double u, double v) const override {
        // Brick texture implementation based on UV coordinates
        int brickRow = static_cast<int>(std::floor(v / 0.5)); // Each row is 0.5 units tall
        int brickCol = static_cast<int>(std::floor(u / 1.0)); // Each column is 1.0 units wide

        // Use alternating colors for adjacent bricks
        if ((brickRow + brickCol) % 2 == 0) {
            return Vec3(0.8, 0.3, 0.3);  // Red brick color
        } else {
            return Vec3(0.9, 0.9, 0.9);  // White mortar color
        }
    }

private:
    double mortarWidth;
};
#endif // TEXTURE_H
