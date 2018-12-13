#ifndef _RAY_H_
#define _RAY_H_

#include "vec3.h"

class Ray {
    public:
        Vec3 org;
        Vec3 dir;
        Ray() {
            org = Vec3();
            dir = Vec3();
        }
        Ray(Vec3 org, Vec3 dir) {
            this -> org = org;
            this -> dir = dir;
        }
};

#endif