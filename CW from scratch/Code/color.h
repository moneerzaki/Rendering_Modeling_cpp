#ifndef COLOR_H
#define COLOR_H

#include "Vec3.h"

class Color {
public:
    Color() : r(0.0), g(0.0), b(0.0) {}
    Color(json colorData) : r(colorData[0]), g(colorData[1]), b(colorData[2]) {}
    Color(double red, double green, double blue) : r(red), g(green), b(blue) {}
    // copy constructor
    Color(const Color& other) : r(other.r), g(other.g), b(other.b) {}

    double red() const { return r; }
    double green() const { return g; }
    double blue() const { return b; }
    // Some color operations
    Color operator*(double scalar) const {
        return Color(r * scalar, g * scalar, b * scalar);
    }
    Color operator*(const Color& other) const {
        return Color(r * other.r, g * other.g, b * other.b);
    }
    // operator overloading for division
    Color operator/(double scalar) const {
        double invScalar = 1.0 / scalar;
        return Color(r * invScalar, g * invScalar, b * invScalar);
    }

    // operator overloading for adding a color and a vector
    Color operator+(const Vec3& other) const {
        return Color(r + other.x, g + other.y, b + other.z);
    }

    // operator overloading for adding two colors
    Color operator+(const Color& other) const {
        return Color(r + other.r, g + other.g, b + other.b);
    }

    // operator overloading for subtracting two colors
    Color operator-(const Color& other) const {
        return Color(r - other.r, g - other.g, b - other.b);
    }

    // operator overloading for this funciton pixel_color += binaryRayColor(r, world);
    Color& operator+=(const Color& other) {
        r += other.r;
        g += other.g;
        b += other.b;
        return *this;
    }
    
    // You can add more color operations as needed
    // Color operator<<(const Color& other) const {
    //     return Color(r << other.r, g << other.g, b << other.b);
    // }

    void print() const {
        std::cout << "Color: (" << r << ", " << g << ", " << b << ")" << std::endl;
    }

private:
    double r, g, b;
};


inline Color operator*(double scalar, const Color& color) {
    return color * scalar;
}
#endif // COLOR_H
