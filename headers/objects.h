#ifndef _OBJECTS_H_
#define _OBJECTS_H_

#include "vec3.h"

const double doubleINF = (double)1e10;

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

class Triangle {
    public:
        Triangle(const Vec3 &a, const Vec3 &b, const Vec3 &c) : pa(a), pb(b), pc(c){
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
            if (det2 / det1 > 0.0 && det3 / det1 >= 0.0 && det4 / det1 >= 0.0 && (det3 + det4) / det1 <= 1.0) {
                len = det2 / det1;
                normalvector = ((pa - pb).cross(pa - pc)).normalize();
                return true;
            }
            return false;
        }

        Vec3 normalvec() {
            // return a normalized normal vector of this triangle
            // the direction is not guaranteed
            return ((pa - pb).cross(pa - pc)).normalize();
        }
    private:
        Vec3 pa;
        Vec3 pb;
        Vec3 pc;
};

class TriangleMesh : public Object {
    public:
        TriangleMesh(char * filename, double a, const Vec3 &b, const Vec3 &sc, 
            const double &refl = 0, const double &transp = 0, const Vec3 &ec = 0) :
            surfaceColor(sc), transparency(transp), reflection(refl) {
                freopen(filename, "r", stdin);
                double ax,ay,az,bx,by,bz,cx,cy,cz;
                int i;
                while ((i = scanf("begin\nvertex %lf %lf %lf\nvertex %lf %lf %lf\nvertex %lf %lf %lf\nend\n\n", &ax, &ay, &az, &bx, &by, &bz, &cx, &cy, &cz)) != EOF ) {
                    //printf("begin\nvertex %lf %lf %lf\nvertex %lf %lf %lf\nvertex %lf %lf %lf\nend\n\n", ax, ay, az, bx, by, bz, cx, cy, cz);
                    Triangle tri(Vec3(ax,ay,az) * a + b, Vec3(bx,by,bz) * a + b, Vec3(cx,cy,cz) * a + b);
                    //printf("pa: %lf, %lf, %lf\n", tri.pa.x, tri.pa.y, tri.pa.z);
                    //printf("pb: %lf, %lf, %lf\n", tri.pb.x, tri.pb.y, tri.pb.z);
                    //printf("pc: %lf, %lf, %lf\n", tri.pc.x, tri.pc.y, tri.pc.z);
                    triangles.push_back(tri);
                }
                fclose(stdin);
        }

        bool intersect(const Vec3 &rayorig, const Vec3 &raydir, double &len, Vec3 &normalvector) const {
            double len_temp;
            Vec3 normalvector_temp;
            int counter = 0;
            bool res = false;
            len = doubleINF;
            for(int i = 0; i < triangles.size(); i++) {
                if(triangles[i].intersect(rayorig, raydir, len_temp, normalvector_temp)) {
                    if (len_temp < len) {
                        len = len_temp;
                        normalvector = normalvector_temp;
                        res = true;
                        counter++;
                    }
                }
            }
            if (res) {
                if (counter % 2 == 0) {
                    // orig outside object
                    normalvector = normalvector * ((normalvector.dot(raydir)<0)?1:-1);
                } else {
                    // orig inside object
                    normalvector = normalvector * ((normalvector.dot(raydir)<0)?-1:1);
                }
            }
            return res;
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
        vector<Triangle> triangles;
        Vec3 surfaceColor;
        double transparency;
        double reflection;
};

#endif