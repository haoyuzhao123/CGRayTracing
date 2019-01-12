#ifndef _BEZIER_H_
#define _BEZIER_H_

#include<vector>
#include<cmath>
#include<cstdio>

using namespace std;

#include "vec3.h"
#include "objects.h"
#include "sampling.h"

// a bezier curve rotate from y-axis
// the curve should not have degree larger than 6

const double Cni[7][7] = {{1,0,0,0,0,0,0},
                    {1,1,0,0,0,0,0},
                    {1,2,1,0,0,0,0},
                    {1,3,3,1,0,0,0},
                    {1,4,6,4,1,0,0},
                    {1,5,10,10,5,1,0},
                    {1,6,15,20,15,6,1}};

const int NEWTON_MAX_ITER = 100;
const double NEWTON_STOP_EPS = 1e-6;
const int num_of_samples_newton = 10;
//const double doubleeps = 1e-4;

double B(int n, int i, double t) {
    if (i > n || i < 0) {
        return 0;
    }
    return Cni[n][i] * pow(1-t, n-i) * pow(t,i);
}

double dB(int n, int i, double t) {
    return B(n-1, i-1, t) * (double)i - B(n-1, i, t) * (double)(n-i);
    //return Cni[n][k] * (k * pow(u, k - 1) * pow(1 - u, n - k) - (n - k) * pow(u, k) * pow(1 - u, n - k - 1));
}

class Bezier : public Object {
    public:
        Bezier(vector<Vec3> points, const Vec3 &pos, const Vec3 &sc, 
            const double &refl = 0, const double &transp = 0, int typeofdata = 0, const Vec3 &ec = 0) : cpoints(points), position(pos), surfaceColor(sc), transparency(transp), reflection(refl) {
            if(cpoints.size() > 6) {
                printf("Num of control points exceeds.\n");
            }
            //build bounding box
            double max_z = -doubleINF;
            double max_y = -doubleINF;
            double min_y = doubleINF;
            for (int i = 0; i < cpoints.size(); i++) {
                if(cpoints[i].z > max_z) {
                    max_z = cpoints[i].z;
                }
                if(cpoints[i].y > max_y) {
                    max_y = cpoints[i].y;
                }
                if(cpoints[i].y < min_y) {
                    min_y = cpoints[i].y;
                }
            }
            xmax = max_z + position.x;
            xmin = -max_z + position.x;
            ymax = max_y + position.y;
            ymin = min_y + position.y;
            zmax = max_z + position.z;
            zmin = -max_z + position.z;
            printf("max_z: %lf, max_y: %lf, min_y: %lf\n", max_z, max_y, min_y);
        }
        bool intersect_with_box(const Vec3 &rayorig, const Vec3 &raydir, double &len) const {
            len = doubleINF;
            double t;
            bool flag = false;
            Vec3 point;
            t = (xmax - rayorig.x) / raydir.x;
            point = rayorig + raydir * t;
            if (t > 0 && point.y >= ymin - epsdouble && point.y <= ymax + epsdouble && point.z >= zmin - epsdouble && point.z <= zmax + epsdouble) {
                if (t < len) {
                    len = t;
                    flag = true;
                }
            }
            t = (xmin - rayorig.x) / raydir.x;
            point = rayorig + raydir * t;
            if (t > 0 && point.y >= ymin - epsdouble && point.y <= ymax + epsdouble && point.z >= zmin - epsdouble && point.z <= zmax + epsdouble) {
                if (t < len) {
                    len = t;
                    flag = true;
                }
            }
            t = (ymax - rayorig.y) / raydir.y;
            point = rayorig + raydir * t;
            if (t > 0 && point.x >= xmin - epsdouble && point.x <= xmax + epsdouble && point.z >= zmin - epsdouble && point.z <= zmax + epsdouble) {
                if (t < len) {
                    len = t;
                    flag = true;
                }
            }
            t = (ymin - rayorig.y) / raydir.y;
            point = rayorig + raydir * t;
            if (t > 0 && point.x >= xmin - epsdouble && point.x <= xmax + epsdouble && point.z >= zmin - epsdouble && point.z <= zmax + epsdouble) {
                if (t < len) {
                    len = t;
                    flag = true;
                }
            }
            t = (zmax - rayorig.z) / raydir.z;
            point = rayorig + raydir * t;
            if (t > 0 && point.x >= xmin - epsdouble && point.x <= xmax + epsdouble && point.y >= ymin - epsdouble && point.y <= ymax + epsdouble) {
                if (t < len) {
                    len = t;
                    flag = true;
                }
            }
            t = (zmin - rayorig.z) / raydir.z;
            point = rayorig + raydir * t;
            if (t > 0 && point.x >= xmin - epsdouble && point.x <= xmax + epsdouble && point.y >= ymin - epsdouble && point.y <= ymax + epsdouble) {
                if (t < len) {
                    len = t;
                    flag = true;
                }
            }
            return flag;
        }
        Vec3 valueP(double u) const {
            Vec3 res = Vec3();
            int n = cpoints.size();
            for (int i = 0; i < n; i++) {
                res = res + cpoints[i] * B(n-1, i, u);
            }
            return res;
        }
        Vec3 gradP(double u) const {
            Vec3 res = Vec3();
            int n = cpoints.size();
            for (int i = 0; i < n; i++) {
                res = res + cpoints[i] * dB(n-1, i, u);
            }
            return res;
        }
        // paras = (t,u,theta)
        Vec3 funcValue(const Vec3 & paras, const Vec3 &rayorig, const Vec3 &raydir) const {
            Vec3 temp = valueP(paras.y);
            temp.x = temp.z * sin(paras.z);
            temp.z *= cos(paras.z);
            return rayorig + raydir * paras.x - position - temp;
        }
        void gradValue(const Vec3 & paras, const Vec3 &rayorig, const Vec3 &raydir, Vec3 &resa, Vec3 &resb, Vec3 &resc) const {
            resa = raydir;
            Vec3 temp1 = gradP(paras.y);
            Vec3 temp2 = valueP(paras.y);
            resb.x = -sin(paras.z) * temp1.z;
            resb.y = -temp1.y;
            resb.z = -cos(paras.z) * temp1.z;
            resc.x = -cos(paras.z) * temp2.z;
            resc.y = 0;
            resc.z = sin(paras.z) * temp2.z;
            //resb = gradP(paras.y) * Vec3(-sin(paras.z), -1, -cos(paras.z));
            //resc = valueP(paras.y) * Vec3(-cos(paras.z), 0, sin(paras.z));
        }
        Vec3 newtonMethod(const Vec3 &initial, const Vec3 &rayorig, const Vec3 &raydir) const {
            Vec3 res = initial;
            int counter = 0;
            Vec3 a,b,c,d,e,f;
            //printf("test\n");
            Vec3 funcval = funcValue(res, rayorig, raydir);
            //funcval.print();
            while(funcval.norm() > NEWTON_STOP_EPS && counter < NEWTON_MAX_ITER) {
                counter++;
                gradValue(res, rayorig, raydir, a, b, c);
                /*
                printf("newiter\n");
                a.print();
                b.print();
                c.print();
                gradP(res.y).print();
                */
                bool flag = inv(a,b,c,d,e,f);
                if(!flag) {
                    //printf("no inverse for the Jacobi matrix.\n");
                    res = res + Vec3(uniform_sampling_zeroone(),uniform_sampling_zeroone(), uniform_sampling_zeroone()) * 0.2 - 0.1;
                    //a.print();
                    //b.print();
                    //c.print();
                    //res.print();
                    //break;
                }
                /*
                printf("print inverse matrix\n");
                d.print();
                e.print();
                f.print();
                res.print();
                */
                res = res - matrixVectorProduct(d,e,f,funcval);
                funcval = funcValue(res, rayorig, raydir);
                //printf("iter: %d, err: %lf\n", counter, funcval.norm());
                //res.print();
                if (funcval.norm() < NEWTON_STOP_EPS) {
                    /*
                    printf("test\n");
                    funcval.print();
                    valueP(res.y).print();
                    printf("sin: %lf, cos: %lf\n", sin(res.z), cos(res.z));
                    funcValue(res,rayorig,raydir).print();
                    printf("%lf", sin(10));
                    */
                }
            }
            //res.print();
            return res;
        }
        Vec3 normalvec(const Vec3 &paras) const {
            Vec3 resp = gradP(paras.y).normalize();
            Vec3 res;
            res.x = resp.x;
            res.z = resp.y;
            res.y = -resp.z;
            res.x = res.z * sin(paras.z);
            res.z = res.z * cos(paras.z);
            return res;
        }
        bool intersect(const Vec3 &rayorig, const Vec3 &raydir, double &len, Vec3 &normalvector) const {
            double t;
            if (!intersect_with_box(rayorig, raydir, t)) {
                return false;
            }
            bool flag = false;
            double a,b,c;
            len = doubleINF;
            for (int i = 0; i < num_of_samples_newton; i++) {
            //a = uniform_sampling_zeroone() * 2 - 1;
            //a = 30;
            b = uniform_sampling_zeroone();
            Vec3 valp = valueP(b);
            double h = valp.y;
            t = 20 + 10 * uniform_sampling_zeroone();
            Vec3 point = rayorig + raydir * t;
            point = point - position;
            double theta;
            if(point.z < 0) {
                theta = 3.14159265 + atan(point.x / point.z);
            } else {
                theta = atan(point.x / point.z);
            }
            Vec3 initial = Vec3(t, b, theta);
            Vec3 res = newtonMethod(initial, rayorig, raydir);
            /*
            if (funcValue(res, rayorig, raydir).norm() < doubleeps) {
            printf("print res and funcval\n");
            res.print();
            funcValue(res,rayorig,raydir).print();
            }
            */
            if((funcValue(res, rayorig, raydir).norm() < doubleeps) && (res.x > 0) && (res.y <= 1) && (res.y >= 0)) {
                //printf("Test\n");
                //printf("%lf %lf\n", len, res.x);
                if (res.x < len) {
                    len = res.x;
                    normalvector = normalvec(res);
                    flag = true;
                    //printf("test\n");
                }
            //len = 30;
            //normalvector = Vec3(0,0,-1);
            //printf("intersect!\n");
            //normalvector.print();
                }
            }
            normalvector = normalvector * ((normalvector.dot(raydir) < 0) ? 1 : -1);
            /*
            if(flag) {
                printf("lenth: %lf\n", len);
            }
            if(flag)
                printf("%d\n", flag);
                */
            return flag;
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
        vector<Vec3> cpoints;
        Vec3 position;
        Vec3 surfaceColor;
        double transparency;
        double reflection;
        double xmax;
        double xmin;
        double ymax;
        double ymin;
        double zmax;
        double zmin;
};

#endif