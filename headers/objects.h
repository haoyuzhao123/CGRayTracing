#ifndef _OBJECTS_H_
#define _OBJECTS_H_

#include "vec3.h"

class Object {
    public:
        //virtual Vec3 normalvec(const Vec3 &point) const {};
        virtual bool intersect(const Vec3 &rayorig, const Vec3 &raydir, double &len, Vec3 & normalvector) const {};
        virtual double getTransparency() const {};
        virtual double getReflection() const {};
        virtual Vec3 getSurfaceColor() const {};
};

class Sphere : public Object{
    public:
        Sphere( 
            const Vec3 &c, 
            const double &r, 
            const Vec3 &sc, 
            const double &refl = 0, 
            const double &transp = 0, 
            const Vec3 &ec = 0) : 
            center(c), radius(r), radius2(r * r), surfaceColor(sc), 
            transparency(transp), reflection(refl) { 
            /* empty */ 
        }

        Vec3 normalvec(const Vec3 &point) const {
            Vec3 ret = point - this -> center;
            return ret.normalize();
        }

        bool intersect(const Vec3 &rayorig, const Vec3 &raydir, double &len, Vec3 &normalvector) const { 
            Vec3 l = center - rayorig; 
            double tca = l.dot(raydir); 
            double l2 = l.dot(l);
            if (tca < 0 && l2 > radius2) {
                return false; 
            }
            double d2 = l.dot(l) - tca * tca; 
            if (d2 > radius2) {
                return false; 
            }
            double thc = sqrt(radius2 - d2); 
            double t0 = tca - thc; 
            double t1 = tca + thc; 
            if (t0 < 0) {
                len = t1;
            }
            else {
                len = t0;
            }
            Vec3 intersection = rayorig + raydir * len;
            normalvector = normalvec(intersection);
            return true; 
        } 

        double getTransparency() const {
            return this -> transparency;
        }

        double getReflection() const {
            return this -> reflection;
        }

        Vec3 getSurfaceColor() const {
            return (this -> surfaceColor).copy();
        }

    private:
        Vec3 center;
        double radius;
        double radius2;
        Vec3 surfaceColor;
        double transparency;
        double reflection;
};

class Triangle : public Object {
    public:
        Triangle(const Vec3 &a, const Vec3 &b, const Vec3 &c, const Vec3 &sc, 
            const double &refl = 0, const double &transp = 0, const Vec3 &ec = 0) : pa(a), pb(b), pc(c), surfaceColor(sc), transparency(transp), reflection(refl) {
            }

        bool intersect(const Vec3 &rayorig, const Vec3 &raydir, double &len, Vec3 &normalvector) const {
            Vec3 res();
            Vec3 e1 = pa - pb;
            Vec3 e2 = pa - pc;
            Vec3 s = pa - rayorig;
            double det1 = det(raydir, e1, e2);
            double det2 = det(s, e1, e2);
            double det3 = det(raydir, s, e2);
            double det4 = det(raydir, e1, s);
            if (det2 / det1 > 0.0 && det3 / det1 >= 0.0 && det3 / det1 >= 0.0 && (det3 + det4) / det1 <= 1.0) {
                len = det2 / det1;
                normalvector = ((pa - pb).cross(pa - pc)).normalize();
                return true;
            }
            return false;
        }
        double getTransparency() const {
            return this -> transparency;
        }

        double getReflection() const {
            return this -> reflection;
        }

        Vec3 getSurfaceColor() const {
            return (this -> surfaceColor).copy();
        }
    private:
        Vec3 pa;
        Vec3 pb;
        Vec3 pc;
        Vec3 surfaceColor;
        double transparency;
        double reflection;
};

#endif