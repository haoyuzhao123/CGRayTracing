#ifndef _TEXTURE_H_
#define _TEXTURE_H_

#include<vector>
#include<cmath>

using namespace std;

#include "vec3.h"
#include "objects.h"

const double texteps = 1e-2;

class Texture {
    public:
        Texture() {
            htexture = false;
        }
        Texture(const vector<vector<Vec3> > &d, const Vec3 &n, const Vec3 &p, double lx, double ly, bool flag = false) {
            data = d;
            normal = n;
            position = p;
            lenx=  lx;
            leny = ly;
            htexture = true;
            isbump = flag;
            height = vector<vector<double> >(d.size(), vector<double>(d[0].size(), 0.0));
            double coeff = 0.5;
            if(isbump) {
                for (int i = 0; i < d.size(); i++) {
                    for (int j = 0; j < d[0].size(); j++) {
                        height[i][j] = (0.299 * d[i][j].x + 0.587 * d[i][j].y + 0.114 * d[i][j].z);
                        height[i][j] = 1 - exp(-3.3 * height[i][j]);
                        height[i][j] *= coeff;
                    }
                }
            }
        }
        bool color(const Vec3 &point, Vec3 &color) const {
            if(htexture == false) {
                return false;
            }
            Vec3 d = point - position;
            d = d - normal * (d.dot(normal));
            if(d.x < texteps && d.x > -texteps) {
                if(0 < d.y && d.y < lenx && 0 < d.z && d.z < leny) {
                    int id1 = (int)floor(d.y / lenx * data.size());
                    int id2 = (int)floor(d.z / leny * data[0].size());
                    color = data[id1][id2];
                    return true;
                }
                return false;
            } else if(d.y < texteps && d.y > -texteps) {
                if(0 < d.x && d.x < lenx && 0 < d.z && d.z < leny) {
                    int id1 = (int)floor(d.x / lenx * data[0].size());
                    int id2 = (int)floor(d.z / leny * data.size());
                    color = data[id2][id1];
                    return true;
                }
                return false;
            } else if(d.z < texteps && d.z > -texteps) {
                if(0 < d.x && d.x < lenx && 0 < d.y && d.y < leny) {
                    int id1 = (int)floor(d.x / lenx * data[0].size());
                    int id2 = (int)floor(d.y / leny * data.size());
                    color = data[data.size() - 1 - id2][id1];
                    return true;
                }
                return false;
            } else {
                return false;
            }
        }

        bool isbump;
        vector<vector<double> > height;
        vector<vector<Vec3> > data;
        Vec3 normal;
        Vec3 position;
        double lenx;
        double leny;
        bool htexture;
};

#endif