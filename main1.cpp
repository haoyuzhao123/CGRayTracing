// The main cpp file for the graphics project
#include <cstdlib> 
#include <cstdio> 
#include <cmath> 
#include <fstream> 
#include <vector> 
#include <iostream> 
#include <cassert> 

#include "headers/vec3.h"
#include "headers/objects.h"

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

//define constants
const int width = 1024;
const int height = 768;

char image_data[width * height * 3];

void write_image(char * filename) {
    stbi_write_png(filename, width, height, 3, image_data, width * 3);
}

#define MAX_RAY_DEPTH 5 
 
float mix(const float &a, const float &b, const float &mix) 
{ 
    return b * mix + a * (1 - mix); 
} 

Vec3 trace( 
    const Vec3 &rayorig, 
    const Vec3 &raydir, 
    const std::vector<Sphere> &spheres, 
    const int &depth) 
{ 
    //if (raydir.length() != 1) std::cerr << "Error " << raydir << std::endl;
    double tnear = INFINITY; 
    const Object* sphere = NULL; 
    // find intersection of this ray with the sphere in the scene
    for (unsigned i = 0; i < spheres.size(); ++i) { 
        double len = INFINITY; 
        if (spheres[i].intersect(rayorig, raydir, len)) { 
            if (len < tnear) { 
                tnear = len; 
                sphere = &spheres[i]; 
            } 
        } 
    } 
    // if there's no intersection return black or background color
    if (!sphere) return Vec3(1); 
    Vec3 surfaceColor = 0; // color of the ray/surfaceof the object intersected by the ray 
    Vec3 phit = rayorig + raydir * tnear; // point of intersection 
    Vec3 nhit = sphere -> normalvec(phit); // normal at the intersection point 
    nhit.normalize(); // normalize normal direction 
    // If the normal and the view direction are not opposite to each other
    // reverse the normal direction. That also means we are inside the sphere so set
    // the inside bool to true. Finally reverse the sign of IdotN which we want
    // positive.
    float bias = 1e-4; // add some bias to the point from which we will be tracing 
    bool inside = false; 
    if (raydir.dot(nhit) > 0) nhit = -nhit, inside = true; 
    if ((sphere->getTransparency() > 0 || sphere->getReflection() > 0) && depth < MAX_RAY_DEPTH) { 
        float facingratio = -raydir.dot(nhit); 
        // change the mix value to tweak the effect
        float fresneleffect = mix(pow(1 - facingratio, 3), 1, 0.1); 
        // compute reflection direction (not need to normalize because all vectors
        // are already normalized)
        Vec3 refldir = raydir - nhit * 2 * raydir.dot(nhit); 
        refldir.normalize(); 
        Vec3 reflection = trace(phit + nhit * bias, refldir, spheres, depth + 1); 
        Vec3 refraction = 0; 
        // if the sphere is also transparent compute refraction ray (transmission)
        if (sphere->getTransparency()) { 
            float ior = 1.1, eta = (inside) ? ior : 1 / ior; // are we inside or outside the surface? 
            float cosi = -nhit.dot(raydir); 
            float k = 1 - eta * eta * (1 - cosi * cosi); 
            Vec3 refrdir = raydir * eta + nhit * (eta *  cosi - sqrt(k)); 
            refrdir.normalize(); 
            refraction = trace(phit - nhit * bias, refrdir, spheres, depth + 1); 
        } 
        // the result is a mix of reflection and refraction (if the sphere is transparent)
        surfaceColor = ( 
            reflection * fresneleffect + 
            refraction * (1 - fresneleffect) * sphere -> getTransparency()) * sphere->getSurfaceColor(); 
    } 
    else { 
        // it's a diffuse object, no need to raytrace any further
        for (unsigned i = 0; i < spheres.size(); ++i) { 
            if (spheres[i].getEmissionColor().x > 0) { 
                // this is a light
                Vec3 transmission = 1; 
                Vec3 lightDirection = spheres[i].getLightDirection(phit); 
                lightDirection.normalize(); 
                for (unsigned j = 0; j < spheres.size(); ++j) { 
                    if (i != j) { 
                        double len; 
                        if (spheres[j].intersect(phit + nhit * bias, lightDirection, len)) { 
                            transmission = 0; 
                            break; 
                        } 
                    } 
                } 
                surfaceColor = surfaceColor + sphere->getSurfaceColor() * transmission * 
                std::max(0.0, nhit.dot(lightDirection)) * spheres[i].getEmissionColor(); 
            } 
        } 
    } 
 
    return surfaceColor + sphere->getEmissionColor(); 
} 

void render(const std::vector<Sphere> &spheres) { 
    Vec3 *image = new Vec3[width * height], *pixel = image; 
    /*
    float invWidth = 1 / float(width), invHeight = 1 / float(height); 
    float fov = 30, aspectratio = width / float(height); 
    float angle = tan(M_PI * 0.5 * fov / 180.); 
    */
    // Trace rays
    for (unsigned y = 0; y < height; ++y) { 
        for (unsigned x = 0; x < width; ++x, ++pixel) { 
            //float xx = (2 * ((x + 0.5) * invWidth) - 1) * angle * aspectratio; 
            //float yy = (1 - 2 * ((y + 0.5) * invHeight)) * angle; 
            //the width (y-axis) of the picture is [-10,10], the picture lies at z=0
            double xx = (2 * double(x) / width - 1) * 10;
            double yy = (2 * double(y) / height - 1) * 10 * height / width;
            // position of the camera is (0,0,zz)
            double zz = 10;
            Vec3 raydir(xx, yy, -10); 
            raydir.normalize(); 
            *pixel = trace(Vec3(0,0,zz), raydir, spheres, 0); 
        } 
    } 
    // Save result to a PPM image (keep these flags if you compile under Windows)
    /*
    std::ofstream ofs("./untitled.ppm", std::ios::out | std::ios::binary); 
    ofs << "P6\n" << width << " " << height << "\n255\n"; 
    for (unsigned i = 0; i < width * height; ++i) { 
        ofs << (unsigned char)(std::min(float(1), image[i].x) * 255) << 
               (unsigned char)(std::min(float(1), image[i].y) * 255) << 
               (unsigned char)(std::min(float(1), image[i].z) * 255); 
    } 
    ofs.close(); 
    delete [] image; 
    */
    unsigned idx = 3 * width * (height+1);
    for (unsigned i = 0; i < width * height; ++i) { 
        if (idx % (3 * width) == 0) {
            idx -= 2 * 3 * width;
        }
        image_data[idx++] =  (unsigned char)(std::min(1.0, image[i].x) * 255);
        image_data[idx++] =  (unsigned char)(std::min(1.0, image[i].y) * 255);
        image_data[idx++] =  (unsigned char)(std::min(1.0, image[i].z) * 255); 
        /*
        if (i < width * height - i) {
            image_data[3*i] =  0;
            image_data[3*i+1] =  0;
            image_data[3*i+2] =  0; 
        } else {
            image_data[3*i] =  255;
            image_data[3*i+1] =  255;
            image_data[3*i+2] =  255; 
        }
        */
    }
    delete [] image; 
    stbi_write_png("test.png", width, height, 3, image_data, width * 3);
} 

int main(int argc, char **argv) 
{ 
    srand(13); 
    std::vector<Sphere> spheres; 
    // position, radius, surface color, reflectivity, transparency, emission color
    //spheres.push_back(Sphere(Vec3( 10100, 0, -0), 10000, Vec3(0.20, 0.20, 0.20), 0, 0.0));
    //spheres.push_back(Sphere(Vec3( -10100, 0, -0), 10000, Vec3(0.20, 0.20, 0.20), 0, 0.0));
    //spheres.push_back(Sphere(Vec3( 0.0, 0, 10100), 10000, Vec3(0.20, 0.20, 0.20), 0, 0.0));
    //spheres.push_back(Sphere(Vec3( 0.0, 0, -10100), 10000, Vec3(0.20, 0.20, 0.20), 0, 0.0)); 
    spheres.push_back(Sphere(Vec3( 0.0, 10100, -0), 10000, Vec3(0.1, 0.1, 0.1), 0, 0.0)); 
    spheres.push_back(Sphere(Vec3( 0.0, -10008, -0), 10000, Vec3(0.20, 0.20, 0.20), 0, 0.0)); 
    spheres.push_back(Sphere(Vec3( 0.0,      -4, -10),     4, Vec3(1.00, 0.32, 0.36), 1, 0.5)); 
    spheres.push_back(Sphere(Vec3( 5.0,     -5, -5),     2, Vec3(0.90, 0.76, 0.46), 1, 0.0)); 
    spheres.push_back(Sphere(Vec3( 5.0,      -4, -15),     3, Vec3(0.65, 0.77, 0.97), 1, 0.0)); 
    spheres.push_back(Sphere(Vec3(-5.5,      -4, -5),     3, Vec3(0.90, 0.90, 0.90), 1, 0.0)); 
    // light
    spheres.push_back(Sphere(Vec3( 0.0,     5, -5),     3, Vec3(0.00, 0.00, 0.00), 0, 0.0, Vec3(3))); 
    render(spheres); 
 
    return 0; 
} 
