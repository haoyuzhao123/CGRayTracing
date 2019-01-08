#ifndef _KDTREE_H_
#define _KDTREE_H_

// kd-tree optimization for the tracing process
// will biuld kd-tree for a Triangle Mesh

#include <vector>
#include <algorithm>

using namespace std;

#include "vec3.h"
#include "objects.h"
#include "util.h"

const int Minkdsize = 10;
const double epsdouble = 1e-4;
const double doubleINF = 1e10;

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

        bool intersect(const Vec3 &rayorig, const Vec3 &raydir) {
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

        int intersect_subtree(const Vec3 &rayorig, const Vec3 &raydir, double & len, Vec3 & normalvector, int curID) {
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

#endif