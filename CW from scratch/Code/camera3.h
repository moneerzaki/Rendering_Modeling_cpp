#ifndef CAMERA_H
#define CAMERA_H
//==============================================================================================
// Originally written in 2016 by Peter Shirley <ptrshrl@gmail.com>
//
// To the extent possible under law, the author(s) have dedicated all copyright and related and
// neighboring rights to this software to the public domain worldwide. This software is
// distributed without any warranty.
//
// You should have received a copy (see file COPYING.txt) of the CC0 Public Domain Dedication
// along with this software. If not, see <http://creativecommons.org/publicdomain/zero/1.0/>.
//==============================================================================================

#include <iostream>
#include <fstream>
#include <nlohmann/json.hpp>
#include "json_read.h"
#include "vec3.h"
#include "ray.h"
#include "hittable.h"
#include "color.h"
#include "basics.h"
#include "interval.h"

using json = nlohmann::json;

class camera {
  public:
    double aspect_ratio      = 1.0;  // Ratio of image width over height
    int    image_width       = 100;  // Rendered image width in pixel count
    int    image_height      = 100;  // Rendered image height in pixel count
    int    samples_per_pixel = 1;   // Count of random samples for each pixel
    int    max_depth         = 10;   // Maximum number of ray bounces into scene

    double vfov     = 90;              // Vertical view angle (field of view)
    Vec3 lookfrom = Vec3(0,0,-1);  // Point camera is looking from
    Vec3 lookat   = Vec3(0,0,0);   // Point camera is looking at
    Vec3   vup      = Vec3(0,1,0);     // Camera-relative "up" direction

    double defocus_angle = 0;  // Variation angle of rays through each pixel
    double focus_dist = 10;    // Distance from camera lookfrom point to plane of perfect focus

    // void render(const Hittable& world) {
    void binaryRendering(std::vector<Hittable*> world, std::string ppmFile) {
        initialize();

        std::ofstream output_file(ppmFile);
        output_file << "P3\n" << image_width << ' ' << image_height << "\n255\n";
        // std::cout << "P3\n" << image_width << ' ' << image_height << "\n255\n";

        for (int j = 0; j < image_height; j++) {
            std::clog << "\rScanlines remaining: " << (image_height - j) << ' ' << std::flush;
            for (int i = 0; i < image_width; i++) {
                Color pixel_color(0,0,0);
                for (int sample = 0; sample < samples_per_pixel; sample++) {
                    Ray r = get_ray(i, j);
                    // pixel_color += binaryRayColor(r, max_depth, world);
                    pixel_color += binaryRayColor(r, world);
                }
                // write_color(output_file, pixel_color, samples_per_pixel);
                int ir = static_cast<int>(255.99 * pixel_color.red());
                int ig = static_cast<int>(255.99 * pixel_color.green());
                int ib = static_cast<int>(255.99 * pixel_color.blue());

                output_file << ir << " " << ig << " " << ib << "\n";
            }
        }

        std::clog << "\rDone.                 \n";
        output_file.close();
    }

    // void binaryRendering(const std::vector<Hittable*>& objects, const std::string& outputFileName) {
    //     initialize(); // initialize the camera
    //     std::ofstream ppmFile(outputFileName);

    //     if (!ppmFile.is_open()) {
    //         std::cerr << "Error: Unable to open output file." << std::endl;
    //         return;
    //     }

    //     ppmFile << "P3\n" << width << " " << height << "\n255\n";  // PPM header
    //     for (int y = height - 1; y >= 0; --y) {
    //         // std::cout << "\rScanlines remaining: " << y << ' ' << std::flush;
    //         for (int x = 0; x < width; ++x) {

    //             double u = double(x) / (width - 1);
    //             double v = double(y) / (height - 1);

    //             Ray ray = generateRay(u, v, 0.0);
    //             if ((x == 0 && y == 0) || (x == width - 1 && y == height - 1)) {
    //                 std::cout << "Ray origin: " << ray.getOrigin() << std::endl;
    //                 std::cout << "Ray direction: " << ray.getDirection() << std::endl;
    //             }
    //             // Ray ray = get_ray(x, y, 1, 1);
    //             Color color = binaryRayColor(ray, objects); 

    //             // Color color = binaryRayColor(ray, scene);
    //             int ir = static_cast<int>(255.99 * color.red());
    //             int ig = static_cast<int>(255.99 * color.green());
    //             int ib = static_cast<int>(255.99 * color.blue());

    //             ppmFile << ir << " " << ig << " " << ib << "\n";
    //         }
    //     }

    //     std::cout << "\nDone." << std::endl;
    //     // ppmFile.flush();
    //     ppmFile.close();
    // }

    // Ray generateRay(double x, double y, double time) const {
    //     // Convert (x, y) from pixel coordinates to normalized device coordinates
    //     double ndcX = (2.0 * x - width) / width;
    //     double ndcY = (height - 2.0 * y) / height;

    //     // Calculate the direction vector in camera space
    //     Vec3 directionInCameraSpace(ndcX, ndcY, -1.0);

    //     // Transform the direction to world space
    //     Vec3 direction = transformDirection(directionInCameraSpace);

    //     // Apply defocus effects
    //     Vec3 defocus_offset = defocus_disk_u * (random_double(0, 1) - 0.5) + defocus_disk_v * (random_double(0, 1) - 0.5);
    //     Vec3 new_origin = position + defocus_offset;

    //     // Calculate the adjusted direction
    //     Vec3 adjusted_direction = (direction * focus_dist - defocus_offset).normalize();

    //     // std::cout << "Generated Ray: origin = " << new_origin << ", direction = " << (direction * focus_dist - defocus_offset).normalize() << std::endl;

    //     return Ray(new_origin, adjusted_direction, time);
    // }

    // Vec3 transformDirection(const Vec3& directionInCameraSpace) const {
    //     // Transform the direction from camera space to world space
    //     Vec3 w = (position - lookAt).normalize();
    //     Vec3 u = upVector.cross(w).normalize();
    //     Vec3 v = w.cross(u);

    //     // Calculate the world space direction
    //     Vec3 worldSpaceDirection = (u * directionInCameraSpace.x) +
    //                             (v * directionInCameraSpace.y) +
    //                             (w * directionInCameraSpace.z);

    //     return worldSpaceDirection;
    // }

  private:
    // int    image_height;    // Rendered image height
    Vec3 center;            // Camera center
    Vec3 pixel00_loc;       // Location of pixel 0, 0
    Vec3   pixel_delta_u;   // Offset to pixel to the right
    Vec3   pixel_delta_v;   // Offset to pixel below
    Vec3   u, v, w;         // Camera frame basis vectors
    Vec3   defocus_disk_u;  // Defocus disk horizontal radius
    Vec3   defocus_disk_v;  // Defocus disk vertical radius
    Color backgroundColor;  // 

    void initialize() {

        // background color 
        backgroundColor = Color(0.5,0.5,0.5);
        // image_height = int(image_width / aspect_ratio);
        // image_height = 800; //
        image_height = (image_height < 1) ? 1 : image_height;

        center = lookfrom;

        // Determine viewport dimensions.
        auto theta = degrees_to_radians(vfov);
        auto h = tan(theta/2);
        auto viewport_height = 2 * h * focus_dist;
        auto viewport_width = viewport_height * (double(image_width)/image_height);

        // Calculate the u,v,w unit basis vectors for the camera coordinate frame.
        w = (lookfrom - lookat).normalize();
        u = (cross(vup, w)).normalize();
        v = cross(w, u);

        // Calculate the vectors across the horizontal and down the vertical viewport edges.
        Vec3 viewport_u = viewport_width * u;    // Vector across viewport horizontal edge
        Vec3 viewport_v = viewport_height * -v;  // Vector down viewport vertical edge

        // Calculate the horizontal and vertical delta vectors to the next pixel.
        pixel_delta_u = viewport_u / image_width;
        pixel_delta_v = viewport_v / image_height;

        // Calculate the location of the upper left pixel.
        auto viewport_upper_left = center - (focus_dist * w) - viewport_u/2 - viewport_v/2;
        pixel00_loc = viewport_upper_left + 0.5 * (pixel_delta_u + pixel_delta_v);

        // Calculate the camera defocus disk basis vectors.
        auto defocus_radius = focus_dist * tan(degrees_to_radians(defocus_angle / 2));
        defocus_disk_u = u * defocus_radius;
        defocus_disk_v = v * defocus_radius;
    }

    Color binaryRayColor(const Ray& ray, const std::vector<Hittable*>& objects) {
        HitRecord record;
        for (const auto& object : objects) {
            if (object->hit(ray, 0.001, std::numeric_limits<double>::infinity(), record)) {
                // If the ray hits an object, return its binary color
                return object->binaryColor();
            }
        }
        // If no intersection, return the default background color
        return backgroundColor;
    }

    Ray get_ray(int i, int j) const {
        // Get a randomly-sampled camera ray for the pixel at location i,j, originating from
        // the camera defocus disk.

        auto pixel_center = pixel00_loc + (i * pixel_delta_u) + (j * pixel_delta_v);
        auto pixel_sample = pixel_center + pixel_sample_square();

        auto ray_origin = (defocus_angle <= 0) ? center : defocus_disk_sample();
        auto ray_direction = pixel_sample - ray_origin;

        return Ray(ray_origin, ray_direction, 0.0);
    }

    Vec3 pixel_sample_square() const {
        // Returns a random point in the square surrounding a pixel at the origin.
        auto px = -0.5 + random_double();
        auto py = -0.5 + random_double();
        return (px * pixel_delta_u) + (py * pixel_delta_v);
    }

    Vec3 pixel_sample_disk(double radius) const {
        // Generate a sample from the disk of given radius around a pixel at the origin.
        auto p = radius * random_in_unit_disk();
        // return (p[0] * pixel_delta_u) + (p[1] * pixel_delta_v);
        return (p.x * pixel_delta_u) + (p.y * pixel_delta_v);
    }

    Vec3 defocus_disk_sample() const {
        // Returns a random point in the camera defocus disk.
        auto p = random_in_unit_disk();
        return center + (p.x * defocus_disk_u) + (p.y * defocus_disk_v);
    }

    // Color ray_color(const Ray& r, int depth, const Hittable& world) const {
    //     // If we've exceeded the ray bounce limit, no more light is gathered.
    //     if (depth <= 0)
    //         return Color(0,0,0);

    //     HitRecord rec;
    //     for (const auto& object : world) {
    //         if (object->hit(r, interval(0.001, infinity), rec)) {
    //             Ray scattered;
    //             Color attenuation;
    //             if (rec.mat->scatter(r, rec, attenuation, scattered))
    //                 return attenuation * ray_color(scattered, depth-1, world);
    //             return Color(0,0,0);
    //         }
    //     }
    //     // if (world.hit(r, interval(0.001, infinity), rec)) {
    //     //     Ray scattered;
    //     //     Color attenuation;
    //     //     if (rec.mat->scatter(r, rec, attenuation, scattered))
    //     //         return attenuation * ray_color(scattered, depth-1, world);
    //     //     return Color(0,0,0);
    //     // }

    //     Vec3 unit_direction = (r.getDirection()).normalize();
    //     auto a = 0.5*(unit_direction.y + 1.0);
    //     return (1.0-a)*Color(1.0, 1.0, 1.0) + a*Color(0.5, 0.7, 1.0);
    // }
    
};


#endif
