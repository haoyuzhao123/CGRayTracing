#ifndef _VEC3_H_
#define _VEC3_H_

#include <cmath>
#include <cstdio>

using namespace std;

const double doubleeps = 1e-4;

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
        Vec3 cross (const Vec3 &b) const {
            return Vec3(y*b.z-z*b.y,z*b.x-x*b.z,x*b.y-y*b.x);
        }

        double det(const Vec3 &a, const Vec3 &b, const Vec3 &c) {
            return (a.x * b.y * c.z + b.x * c.y * a.z + c.x * a.y * b.z - a.x * c.y * b.z - b.x * a.y * c.z - c.x * b.y * a.z);
        }

        void print() {
            printf("x: %.6lf, y: %.6lf, z: %.6lf\n", x, y, z);
        }
};

// computing the determinent
double det(const Vec3 &a, const Vec3 &b, const Vec3 &c) {
    return (a.x * b.y * c.z + b.x * c.y * a.z + c.x * a.y * b.z - a.x * c.y * b.z - b.x * a.y * c.z - c.x * b.y * a.z);
}

Vec3 matrixVectorProduct(const Vec3 &a, const Vec3 &b, const Vec3 &c, const Vec3 &d) {
    return a * d.x + b * d.y + c * d.z;
}

bool inv(const Vec3 &a, const Vec3 &b, const Vec3 &c, Vec3 &resa, Vec3 &resb, Vec3 &resc) {
    double d = det(a, b, c);
    if (d < doubleeps && d > -doubleeps) {
        // determinant is 0, no inverse
        return false;
    }
    resa.x = (b.y * c.z - b.z * c.y) / d;
    resa.y = (c.y * a.z - c.z * a.y) / d;
    resa.z = (a.y * b.z - a.z * b.y) / d;
    resb.x = (c.x * b.z - c.z * b.x) / d;
    resb.y = (a.x * c.z - a.z * c.x) / d;
    resb.z = (b.x * a.z - b.z * a.x) / d;
    resc.x = (b.x * c.y - c.x * b.y) / d;
    resc.y = (c.x * a.y - c.y * a.x) / d;
    resc.z = (a.x * b.y - a.y * b.x) / d;
    return true;
}

#endif