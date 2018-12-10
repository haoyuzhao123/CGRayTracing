#ifndef _OBJECTS_H_
#define _OBJECTS_H_

#include "vec3.h"

class Object {
    public:
        virtual bool intersects(const Vec3 &rayorig, const Vec3 &raydir, Vec3 * intersect);
        virtual Vec3 normalvec(const Vec3 &point);
        virtual bool intersect(const Vec3 &rayorig, const Vec3 &raydir, float &t0, float &t1) const;
};

class Sphere : public Object{
    public:
        Vec3 center;
        double radius;
        double radius2;
        Vec3 surfaceColor;
        Vec3 emissionColor;
        double transparency;
        double reflection;

        Sphere( 
        const Vec3 &c, 
        const double &r, 
        const Vec3 &sc, 
        const double &refl = 0, 
        const double &transp = 0, 
        const Vec3 &ec = 0) : 
        center(c), radius(r), radius2(r * r), surfaceColor(sc), emissionColor(ec), 
        transparency(transp), reflection(refl) 
        { /* empty */ }

        bool intersects(const Vec3 &rayorig, const Vec3 &raydir, Vec3 * intersect) {
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

        bool intersect(const Vec3 &rayorig, const Vec3 &raydir, float &t0, float &t1) const 
        { 
        Vec3 l = center - rayorig; 
        float tca = l.dot(raydir); 
        if (tca < 0) return false; 
        float d2 = l.dot(l) - tca * tca; 
        if (d2 > radius2) return false; 
        float thc = sqrt(radius2 - d2); 
        t0 = tca - thc; 
        t1 = tca + thc; 
 
        return true; 
        } 
};

#endif