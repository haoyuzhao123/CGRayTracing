#ifndef _VEC3_H_
#define _VEC3_H_

#include <cmath>

class Vec3 {
    public:
        double x,y,z;

        Vec3() {
            x = 0;
            y = 0;
            z = 0;
        }
        Vec3(double x, double y, double z) {
            this -> x = x;
            this -> y = y;
            this -> z = z;
        }

        double norm() {
            return sqrt(x * x + y * y + z * z);
        }

        Vec3 normalize() {
            double len = norm();
            if (len > 0) {
                return Vec3(x/len, y/len, z/len);
            }
            else {
                return Vec3();
            }
        }
        //multiply with double
        Vec3 operator * (const double & f) const {
            return Vec3(x * f, y * f, z * f);
        }
        //elementwise product
        Vec3 operator * (const Vec3 &v) const {
            return Vec3(x * v.x, y * v.y, z * v.z);
        }
        //dot product
        double dot(const Vec3 &v) {
            return x * v.x + y * v.y + z * v.z;
        }

        Vec3 operator + (const Vec3 &v) const {
            return Vec3(x + v.x, y + v.y, z + v.z);
        }
        Vec3 operator - (const Vec3 &v) const {
            return Vec3(x - v.x, y - v.y, z - v.z);
        }
        Vec3 operator - () const {
            return Vec3(-x, -y, -z);
        }
};

#endif