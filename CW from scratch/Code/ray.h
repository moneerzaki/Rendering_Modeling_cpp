#ifndef RAY_H
#define RAY_H

#include "vec3.h"  // Assuming you have a Vec3 class


class Ray {
public:
    Ray() {}
    Ray(const Vec3& origin, const Vec3& direction) : origin(origin), direction(direction) { t = 0.0;}
    Ray(const Vec3& origin, const Vec3& direction, double time = 0.0) : origin(origin), direction(direction), t(time) {}


    Vec3 getOrigin() const { return origin; }
    Vec3 getDirection() const { return direction.normalize(); }
    double getTime() const { return t; }

    // Get a point along the ray at parameter t
    Vec3 at(double t) const {
        return origin + t * direction;
    }
    Vec3 at(double t, double time) const {
        return origin + t * direction + timeOffset(time);
    }

private:
    Vec3 origin;
    Vec3 direction;
    double t = 0.0;

    Vec3 timeOffset(double time) const {
        return time * Vec3(1.0, 1.0, 1.0);
    }

};

#endif // RAY_H
