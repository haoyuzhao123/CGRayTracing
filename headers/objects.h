#ifndef _OBJECTS_H_
#define _OBJECTS_H_

#include "vec3.h"

class Object {
    public:
        virtual bool intersect(const Vec3 &rayorig, const Vec3 &raydir, Vec3 * intersect);
        virtual Vec3 normalvec(const Vec3 &point);
};

class Sphere {
    public:
        Vec3 center;
        double radius;
        Vec3 surfaceColor;
        Vec3 emmisionColor;
        double transparency;
        double reflection;

        bool intersect(const Vec3 &rayorig, const Vec3 &raydir, Vec3 * intersect) {
            Vec3 l = this -> center - rayorig;
            if (l.norm() < this -> radius) {
                return true;
            }
            double len1 = l.dot(raydir);
            if(len1 < 0) {
                return false;
            }
        }

        Vec3 normalvec(const Vec3 &point) {
            Vec3 ret = point - this -> center;
            return ret.normalize();
        }
};

#endif