#ifndef _OBJECTS_H_
#define _OBJECTS_H_

#include "vec3.h"

class Object {
    public:
        virtual Vec3 normalvec(const Vec3 &point) const;
        virtual bool intersect(const Vec3 &rayorig, const Vec3 &raydir, double &len) const;
        virtual double getTransparency() const;
        virtual double getReflection() const;
        virtual Vec3 getSurfaceColor() const;
        virtual Vec3 getEmissionColor() const;
        virtual Vec3 getLightDirection(const Vec3 &point) const;
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
            center(c), radius(r), radius2(r * r), surfaceColor(sc), emissionColor(ec), 
            transparency(transp), reflection(refl) { 
            /* empty */ 
        }

        Vec3 normalvec(const Vec3 &point) const {
            Vec3 ret = point - this -> center;
            return ret.normalize();
        }

        bool intersect(const Vec3 &rayorig, const Vec3 &raydir, double &len) const { 
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

        Vec3 getEmissionColor() const {
            return (this -> emissionColor).copy();
        }

        Vec3 getLightDirection(const Vec3 &point) const {
            return (this -> center - point).copy();
        }

    private:
        Vec3 center;
        double radius;
        double radius2;
        Vec3 surfaceColor;
        Vec3 emissionColor;
        double transparency;
        double reflection;
};

#endif