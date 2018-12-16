#ifndef _VEC3_H_
#define _VEC3_H_

#include <cmath>

class Vec3 {
    public:
        double x,y,z;
        /*
        Vec3() {
            x = 0;
            y = 0;
            z = 0;
        }
        Vec3(double x) {
            this -> x = x;
            this -> y = x;
            this -> z = x;
        }*/
        Vec3(double x = 0, double y = 0, double z = 0) {
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
                this -> x *= 1 / len;
                this -> y *= 1 / len;
                this -> z *= 1 / len;
            }
            return *this;
        }
        // copy
        Vec3 copy() const {
            return Vec3(this -> x, this -> y, this -> z);
        }
        // operators
        // multiply with double
        Vec3 operator * (const double & f) const {
            return Vec3(x * f, y * f, z * f);
        }
        // elementwise product
        Vec3 operator * (const Vec3 &v) const {
            return Vec3(x * v.x, y * v.y, z * v.z);
        }
        Vec3 mul(const Vec3 &v) const {
            return Vec3(x * v.x, y * v.y, z * v.z);
        }
        // dot product
        double dot(const Vec3 &v) const {
            return x * v.x + y * v.y + z * v.z;
        }

        Vec3 operator + (const Vec3 &v) const {
            return Vec3(x + v.x, y + v.y, z + v.z);
        }
        Vec3 operator + (double b) const {
            return Vec3(x + b, y + b, z + b);
        }
        Vec3 operator - (const Vec3 &v) const {
            return Vec3(x - v.x, y - v.y, z - v.z);
        }
        Vec3 operator - (double b) const {
            return Vec3(x - b, y - b, z - b);
        }
        Vec3 operator - () const {
            return Vec3(-x, -y, -z);
        }
        // cross product
        Vec3 cross (Vec3 &b) {
            return Vec3(y*b.z-z*b.y,z*b.x-x*b.z,x*b.y-y*b.x);
        }
};

#endif