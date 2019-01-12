#ifndef _OBJECTS_H_
#define _OBJECTS_H_

#include <vector>
#include <algorithm>

using namespace std;

#include "vec3.h"
//#include "kdtree.h"


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
        // helper functions for the bounding boxes in KD-Tree
        double max_x() {
            return max(pa.x, pb.x, pc.x);
        }
        double max_y() {
            return max(pa.y, pb.y, pc.y);
        }
        double max_z() {
            return max(pa.z, pb.z, pc.z);
        }
        double min_x() {
            return min(pa.x, pb.x, pc.x);
        }
        double min_y() {
            return min(pa.y, pb.y, pc.y);
        }
        double min_z() {
            return min(pa.z, pb.z, pc.z);
        }

        Vec3 pa;
        Vec3 pb;
        Vec3 pc;
};

const int Minkdsize = 10;
const double epsdouble = 1e-4;
//const double doubleINF = 1e10;

class KDNode {
    public:
        double xmax;
        double xmin;
        double ymax;
        double ymin;
        double zmax;
        double zmin;

        int left, right; // use a vector to store the KDNodes in KDTree, use index to represent the children

        vector< pair<int, Triangle> > triangleList;

        KDNode(const vector< pair<int, Triangle> > &tl) {
            triangleList = tl;
            left = -1;
            right = -1;
        }

        bool intersect(const Vec3 &rayorig, const Vec3 &raydir) const {
            double t;
            Vec3 point;
            t = (xmax - rayorig.x) / raydir.x;
            point = rayorig + raydir * t;
            if (t > 0 && point.y >= ymin - epsdouble && point.y <= ymax + epsdouble && point.z >= zmin - epsdouble && point.z <= zmax + epsdouble) {
                return true;
            }
            t = (xmin - rayorig.x) / raydir.x;
            point = rayorig + raydir * t;
            if (t > 0 && point.y >= ymin - epsdouble && point.y <= ymax + epsdouble && point.z >= zmin - epsdouble && point.z <= zmax + epsdouble) {
                return true;
            }
            t = (ymax - rayorig.y) / raydir.y;
            point = rayorig + raydir * t;
            if (t > 0 && point.x >= xmin - epsdouble && point.x <= xmax + epsdouble && point.z >= zmin - epsdouble && point.z <= zmax + epsdouble) {
                return true;
            }
            t = (ymin - rayorig.y) / raydir.y;
            point = rayorig + raydir * t;
            if (t > 0 && point.x >= xmin - epsdouble && point.x <= xmax + epsdouble && point.z >= zmin - epsdouble && point.z <= zmax + epsdouble) {
                return true;
            }
            t = (zmax - rayorig.z) / raydir.z;
            point = rayorig + raydir * t;
            if (t > 0 && point.x >= xmin - epsdouble && point.x <= xmax + epsdouble && point.y >= ymin - epsdouble && point.y <= ymax + epsdouble) {
                return true;
            }
            t = (zmin - rayorig.z) / raydir.z;
            point = rayorig + raydir * t;
            if (t > 0 && point.x >= xmin - epsdouble && point.x <= xmax + epsdouble && point.y >= ymin - epsdouble && point.y <= ymax + epsdouble) {
                return true;
            }
            return false;
        }
};

bool compareX(pair<int, Triangle> p1, pair<int, Triangle> p2) {
    return p1.second.max_x() < p2.second.max_x();
}
bool compareY(pair<int, Triangle> p1, pair<int, Triangle> p2) {
    return p1.second.max_y() < p2.second.max_y();
}
bool compareZ(pair<int, Triangle> p1, pair<int, Triangle> p2) {
    return p1.second.max_z() < p2.second.max_z();
}

class KDTree {
    public:
        vector<KDNode> kdnodes;

        void buildKdTree(vector<pair<int, Triangle> > sublist, int parID, bool isLeft, int curDIM, bool isRoot = false) {
            int curID = kdnodes.size();
            if (!isRoot) {
                if (isLeft) {
                    kdnodes[parID].left = curID;
                } else {
                    kdnodes[parID].right = curID;
                }
            }
            kdnodes.push_back(KDNode(sublist));
            kdnodes[curID].xmax = -doubleINF;
            kdnodes[curID].ymax = -doubleINF;
            kdnodes[curID].zmax = -doubleINF;
            kdnodes[curID].xmin = doubleINF;
            kdnodes[curID].ymin = doubleINF;
            kdnodes[curID].zmin = doubleINF;
            for (int i = 0; i < sublist.size(); i++) {
                Triangle a = sublist[i].second;
                double max_x, max_y, max_z, min_x, min_y, min_z;
                max_x = a.max_x();
                max_y = a.max_y();
                max_z = a.max_z();
                min_x = a.min_x();
                min_y = a.min_y();
                min_z = a.min_z();
                if (kdnodes[curID].xmax < max_x) kdnodes[curID].xmax = max_x;
                if (kdnodes[curID].ymax < max_y) kdnodes[curID].ymax = max_y;
                if (kdnodes[curID].zmax < max_z) kdnodes[curID].zmax = max_z;
                if (kdnodes[curID].xmin > min_x) kdnodes[curID].xmin = min_x;
                if (kdnodes[curID].ymin > min_y) kdnodes[curID].ymin = min_y;
                if (kdnodes[curID].zmin > min_z) kdnodes[curID].zmin = min_z;
            }
            kdnodes[curID].left = -1;
            kdnodes[curID].right = -1;
            if (sublist.size() < Minkdsize) {
                return;
            }
            if (curDIM == 0) {
                sort(sublist.begin(), sublist.end(), compareX);
            } else if (curDIM == 1) {
                sort(sublist.begin(), sublist.end(), compareY);
            } else {
                sort(sublist.begin(), sublist.end(), compareZ);
            }
            vector< pair<int, Triangle> > leftsublist(sublist.begin(), sublist.begin() + sublist.size() / 2);
            vector< pair<int, Triangle> > rightsublist(sublist.begin() + sublist.size() / 2, sublist.end());
            int newDIM = (curDIM + 1) % 3;
            buildKdTree(leftsublist, curID, true, newDIM);
            buildKdTree(rightsublist, curID, false, newDIM);
            return;
        }

        int intersect_subtree(const Vec3 &rayorig, const Vec3 &raydir, double & len, Vec3 & normalvector, int curID) const {
            if (! (kdnodes[curID].intersect(rayorig, raydir))) {
                return 0;
            }
            if (kdnodes[curID].triangleList.size() < Minkdsize) {
                double len_temp;
                Vec3 normalvector_temp;
                int counter = 0;
                bool res = false;
                len = doubleINF;
                for(int i = 0; i < kdnodes[curID].triangleList.size(); i++) {
                    if(kdnodes[curID].triangleList[i].second.intersect(rayorig, raydir, len_temp, normalvector_temp)) {
                        if (len_temp < len) {
                            len = len_temp;
                            normalvector = normalvector_temp;
                            res = true;
                            counter++;
                        }
                    }
                }
                return counter;
            } else {
                double lenleft, lenright;
                Vec3 normalvecleft, normalvecright;
                int numleft = intersect_subtree(rayorig, raydir, lenleft, normalvecleft, kdnodes[curID].left);
                int numright = intersect_subtree(rayorig, raydir, lenright, normalvecright, kdnodes[curID].right);
                if (numleft > 0) {
                    if(numright > 0) {
                        if(lenleft < lenright) {
                            len = lenleft;
                            normalvector = normalvecleft;
                        } else {
                            len = lenright;
                            normalvector = normalvecright;
                        }
                    } else {
                        len = lenleft;
                        normalvector = normalvecleft;
                    }
                } else {
                    if(numright > 0) {
                        len = lenright;
                        normalvector = normalvecright;
                    }
                }
                return numleft + numright;
            }
        }

        bool intersect(const Vec3 &rayorig, const Vec3 &raydir, double & len, Vec3 & normalvector) const {
            int counter = intersect_subtree(rayorig, raydir, len, normalvector, 0);
            if (counter > 0) {
                if (counter % 2 == 0) {
                    // orig outside object
                    normalvector = normalvector * ((normalvector.dot(raydir)<0)?1:-1);
                } else {
                    // orig inside object
                    normalvector = normalvector * ((normalvector.dot(raydir)<0)?-1:1);
                }
                return true;
            } else {
                return false;
            }
        }
};

class TriangleMesh : public Object {
    public:
        TriangleMesh(char * filename, double a, const Vec3 &b, const Vec3 &sc, 
            const double &refl = 0, const double &transp = 0, int typeofdata = 0, const Vec3 &ec = 0) :
            surfaceColor(sc), transparency(transp), reflection(refl) {
                objtype = typeofdata;
                freopen(filename, "r", stdin);
                if (typeofdata == 0) {
                double ax,ay,az,bx,by,bz,cx,cy,cz;
                int i;
                while ((i = scanf("begin\nvertex %lf %lf %lf\nvertex %lf %lf %lf\nvertex %lf %lf %lf\nend\n\n", &ax, &ay, &az, &bx, &by, &bz, &cx, &cy, &cz)) != EOF ) {
                    //printf("begin\nvertex %lf %lf %lf\nvertex %lf %lf %lf\nvertex %lf %lf %lf\nend\n\n", ax, ay, az, bx, by, bz, cx, cy, cz);
                    Triangle tri(Vec3(ax,ay,-az) * a + b, Vec3(bx,by,-bz) * a + b, Vec3(cx,cy,-cz) * a + b);
                    //printf("pa: %lf, %lf, %lf\n", tri.pa.x, tri.pa.y, tri.pa.z);
                    //printf("pb: %lf, %lf, %lf\n", tri.pb.x, tri.pb.y, tri.pb.z);
                    //printf("pc: %lf, %lf, %lf\n", tri.pc.x, tri.pc.y, tri.pc.z);
                    triangles.push_back(pair<int, Triangle>(triangles.size(), tri));
                }
                } else if (typeofdata == 1) {
                int num;
                int id1,id2,id3;
                double x,y,z;
                scanf("%d\n", &num);
                printf("number of vertices: %d\n", num);
                vector<Vec3> vertices;
                for (int i = 0; i < num; i++) {
                    int j = scanf("v  %lf %lf %lf\n", &x, &y, &z);
                    //printf("%d\n", j);
                    //printf("%lf %lf %lf\n", x, y, z);
                    vertices.push_back(Vec3(x,y,-z));
                }
                scanf("%d\n", &num);
                printf("number of faces: %d\n", num);
                for (int i = 0; i < num; i++) {
                    scanf("f %d %d %d \n", &id1, &id2, &id3);
                    triangles.push_back(pair<int, Triangle>(triangles.size(), Triangle(vertices[id1-1] * a + b, vertices[id2-1] * a + b, vertices[id3-1] * a + b)));
                }
                } else if (typeofdata == 2) {
                int num;
                int id1,id2,id3,i1,i2,i3,i4,i5,i6;
                double x,y,z,u,v;
                scanf("%d\n", &num);
                printf("number of vertices: %d\n", num);
                vector<Vec3> vertices;
                for (int i = 0; i < num; i++) {
                    int j = scanf("v %lf %lf %lf\n", &x, &y, &z);
                    //printf("%d\n", j);
                    //printf("%lf %lf %lf\n", x, y, z);
                    vertices.push_back(Vec3(x,y,-z));
                }
                
                for (int i = 0; i < num; i++) {
                    scanf("vn %lf %lf %lf\n", &x, &y, &z);
                }
                for (int i = 0; i < num; i++) {
                    scanf("vt %lf %lf\n", &x, &y);
                }
                
                scanf("%d\n", &num);
                printf("number of faces: %d\n", num);
                for (int i = 0; i < num; i++) {
                    scanf("f %d/%d/%d %d/%d/%d %d/%d/%d \n", &id1, &i1, &i2, &id2, &i3, &i4, &id3, &i5, &i6);
                    triangles.push_back(pair<int, Triangle>(triangles.size(), Triangle(vertices[id1-1] * a + b, vertices[id2-1] * a + b, vertices[id3-1] * a + b)));
                }
                }
                fclose(stdin);
                kdtree.buildKdTree(triangles, 0, false, 0, true);
        }

        bool intersect(const Vec3 &rayorig, const Vec3 &raydir, double &len, Vec3 &normalvector) const {
            /*
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
            */
            bool res =  kdtree.intersect(rayorig, raydir, len, normalvector);
            if (objtype == 2) {
                normalvector = normalvector * ((normalvector.dot(Vec3(0,1,0)) > 0) ? 1 : -1);
            }
            /*
            if ((rayorig - Vec3(0,0,-10)).norm() < eps && objtype == 2) {
                // light from back of the image and the object is water
                //fprintf(stderr, "\rtest");
                double len_temp;
                len_temp = (0.0 - rayorig.z) / raydir.z;
                Vec3 inter_temp = rayorig + raydir * len_temp;
                double len_temp2;
                Vec3 normalvec_temp2;
                if (!intersect(inter_temp, Vec3(0,-1,0), len_temp2, normalvec_temp2)) {
                    // the intersection point is under the water.
                    len = len_temp;
                    normalvector = Vec3(0,0,-1);
                    res = true;
                }
            }
            */
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
        vector<pair<int, Triangle> > triangles;
        Vec3 surfaceColor;
        double transparency;
        double reflection;
        KDTree kdtree;
        int objtype;
};

#endif