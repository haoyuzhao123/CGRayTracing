#ifndef _SAMPLING_H_
#define _SAMPLING_H_

#include <cstdlib>
#include <cmath>
#include <cstdio>

#include "vec3.h"

// las vegas algorithm to generate the uniform distribution on a sphere
Vec3 uniform_sampling_sphere() {
    while(true) {
        double x = (double)rand() / RAND_MAX * 2.0 - 1;
        double y = (double)rand() / RAND_MAX * 2.0 - 1;
        double z = (double)rand() / RAND_MAX * 2.0 - 1;
        if (x * x + y * y + z * z <= 1) {
            return Vec3(x,y,z).normalize();
        }
    }
}

Vec3 uniform_sampling_halfsphere(const Vec3 &dir) {
    while(true) {
        Vec3 sample = uniform_sampling_sphere();
        if (sample.dot(dir) > 0) {
            return sample;
        }
    }
}

double uniform_sampling_zeroone() {
    return (double)rand() / RAND_MAX;
}

#endif