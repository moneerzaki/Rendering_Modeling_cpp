// vec3.h

#ifndef VEC3_H
#define VEC3_H

#include <cmath>
#include <iostream>
#include "basics.h"

class Vec3 {
public:
    double x, y, z;

    Vec3() : x(0), y(0), z(0) {}
    Vec3(double x, double y, double z) : x(x), y(y), z(z) {}

    // Addition
    Vec3 operator+(const Vec3& other) const {
        return Vec3(x + other.x, y + other.y, z + other.z);
    }

    // Fresnel reflection coefficient 
    static double fresnelReflection(const Vec3& normal, const Vec3& incident, double refractiveIndex) {
        double cosTheta = fmin(-incident.dot(normal), 1.0);
        double sinTheta = sqrt(1.0 - cosTheta * cosTheta);
        double etaI = 1.0;
        double etaT = refractiveIndex;

        if (cosTheta > 0) {
            std::swap(etaI, etaT);
        }

        double sint = etaI / etaT * sinTheta;
        if (sint >= 1.0) {
            // Total internal reflection
            return 1.0;
        } else {
            double cost = sqrt(1.0 - sint * sint);
            cosTheta = fabs(cosTheta);
            double Rs = ((etaT * cosTheta) - (etaI * cost)) / ((etaT * cosTheta) + (etaI * cost));
            double Rp = ((etaI * cosTheta) - (etaT * cost)) / ((etaI * cosTheta) + (etaT * cost));
            return 0.5 * (Rs * Rs + Rp * Rp);
        }
    }

    // static Vec3 reflect(const Vec3& incident, const Vec3& normal) const {
    static Vec3 reflect(const Vec3& incident, const Vec3& normal) {
        return incident - (2 * incident.dot(normal) * normal);
        // return incident - mult(2 * incident.dot(normal), normal);
    } 

    // Vec3 reflect(const Vec3& normal) const {
    //     return *this - 2 * normal.dot(*this) * normal;
    // }

    // old refract function
    static bool refract(const Vec3& incident, const Vec3& normal, double refractiveIndex, Vec3& refracted) {
        Vec3 unit_incident = incident.normalize();
        double cos_theta = fmin(-unit_incident.dot(normal), 1.0);
        double sin_theta = std::sqrt(1.0 - cos_theta * cos_theta);

        Vec3 ref_direction = refractiveIndex * (unit_incident + cos_theta * normal) - sin_theta * normal;

        refracted = ref_direction.normalize();
        return true;
    }

    // new refract function
    // static bool refract(const Vec3& incident, const Vec3& normal, double refractiveIndex, Vec3& refracted) {
    //     Vec3 unitIncident = incident.normalize();
    //     double cosTheta = fmin(-unitIncident.dot(normal), 1.0);
    //     double sinTheta = sqrt(1.0 - cosTheta * cosTheta);
    //     double etaI = 1.0;
    //     double etaT = refractiveIndex;

    //     if (cosTheta > 0) {
    //         std::swap(etaI, etaT);
    //     }

    //     double sint = etaI / etaT * sinTheta;
    //     if (sint >= 1.0) {
    //         // Total internal reflection
    //         return false;
    //     } else {
    //         double cost = sqrt(1.0 - sint * sint);
    //         refracted = (etaI / etaT) * unitIncident + ((etaI / etaT) * cosTheta - cost) * normal;
    //         return true;
    //     }
    // }

    // Subtraction
    Vec3 operator-(const Vec3& other) const {
        return Vec3(x - other.x, y - other.y, z - other.z);
    }

    // Element-wise multiplication
    Vec3 operator*(const Vec3& other) const {
        return Vec3(x * other.x, y * other.y, z * other.z);
    }

    // Scalar multiplication
    Vec3 operator*(double scalar) const {
        return Vec3(x * scalar, y * scalar, z * scalar);
    }

    // Scalar multiplication (left operand)
    friend Vec3 operator*(double scalar, const Vec3& vector);

    // static Vec3 operator*(double scalar, const Vec3& vector) {
    //     return vector * scalar;
    // }

    // Scalar division
    Vec3 operator/(double scalar) const {
        double invScalar = 1.0 / scalar;
        return Vec3(x * invScalar, y * invScalar, z * invScalar);
    }

    Vec3 operator-() const { return Vec3(-x, -y, -z); }
    
    double length_squared() const {
        return x * x + y * y + z * z;
    }

    // Dot product
    double dot(const Vec3& other) const {
        return x * other.x + y * other.y + z * other.z;
    }

    // Cross product
    Vec3 cross(const Vec3& other) const {
        return Vec3(y * other.z - z * other.y, z * other.x - x * other.z, x * other.y - y * other.x);
    }

    // Length of the vector
    double length() const {
        return std::sqrt(x * x + y * y + z * z);
    }

    // Squared length of the vector
    double squared_length() const {
        return x * x + y * y + z * z;
    }

    // Normalize the vector
    Vec3 normalize() const {
        double len = length();
        if (len != 0) {
            return Vec3(x / len, y / len, z / len);
        } else {
            return *this; // Avoid division by zero
        }
    }

    // Print the vector
    void print() const {
        std::cout << "(" << x << ", " << y << ", " << z << ")" << std::endl;
    }
};

// Scalar multiplication (left operand)
inline Vec3 operator*(double scalar, const Vec3& vector) {
    return vector * scalar;
}

// Output stream operator
inline std::ostream& operator<<(std::ostream& out, const Vec3& vector) {
    out << "(" << vector.x << ", " << vector.y << ", " << vector.z << ")";
    return out;
}

inline Vec3 random_in_unit_disk() {
    while (true) {
        auto p = Vec3(random_double(-1,1), random_double(-1,1), 0);
        if (p.length_squared() < 1)
            return p;
    }
}

inline Vec3 random_in_unit_sphere() {
    Vec3 p;
    do {
        p = 2.0 * Vec3(random_double(), random_double(), random_double()) - Vec3(1, 1, 1);
    } while (p.squared_length() >= 1.0);
    return p;
}

inline Vec3 cross(const Vec3& u, const Vec3& v) {
    return u.cross(v);
}
// inline vec3 cross(const vec3& u, const vec3& v) {
//     return vec3(u.e[1] * v.e[2] - u.e[2] * v.e[1],
//                 u.e[2] * v.e[0] - u.e[0] * v.e[2],
//                 u.e[0] * v.e[1] - u.e[1] * v.e[0]);
// }

#endif // VEC3_H
