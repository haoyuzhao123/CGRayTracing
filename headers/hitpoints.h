#ifndef _HITPOINTS_H_
#define _HITPOINTS_H_

#include "vec3.h"

class Hitpoint {
    public:
        Vec3 f;
        Vec3 pos;
        Vec3 normal;
        Vec3 flux;
        double r2;
        int n;
        int h;
        int w;

        Hitpoint() : f(), pos(), normal(), flux(), r2(), n(), h(), w() {
            // empty
        }
};

#endif