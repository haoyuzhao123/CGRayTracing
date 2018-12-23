#ifndef _HASH_H_
#define _HASH_H_

#include <vector>
#include <cmath>
using namespace std;

#include "hitpoints.h"

//for simplicity, we will assume that all of the hitpoints will be bounded in a 70 * 70 * 70 box, and the box is [-35,35] * [-15,55] * [-25,45]
const double SIZE_OF_SCENE = 70.0;

const double XMAX = 35.0;
const double XMIN = -35.0;
const double YMAX = 55.0;
const double YMIN = -15.0;
const double ZMAX = 45.0;
const double ZMIN = -25.0;

class Hashtable {
    public:
        Hashtable(int hashsize, double celllength) {
            this -> hashsize = hashsize;
            this -> celllength = celllength;
            this -> num_of_cell_per_dim = (int)(ceil(SIZE_OF_SCENE / celllength));
            this -> celllength = SIZE_OF_SCENE / num_of_cell_per_dim;
            hashtable = vector<vector<Hitpoint> >(hashsize, vector<Hitpoint>());

            this -> debug_counter = 0;
        }
        // the spatial hash function
        // this hash function appears in smallppm.cpp
        // this hash function is also recommended on stackoverflow
        // https://stackoverflow.com/questions/5928725/hashing-2d-3d-and-nd-vectors
        unsigned int hash(const int ix, const int iy, const int iz) {
	        return (unsigned int)((ix*73856093)^(iy*19349663)^(iz*83492791))%hashsize;
        }
        void compute_coord(double x, double y, double z, int &ix, int &iy, int &iz) {
            ix = (int)floor((x - XMIN) / celllength);
            iy = (int)floor((y - YMIN) / celllength);
            iz = (int)floor((z - ZMIN) / celllength);
        }
        void insert(Hitpoint hpoint) {
            //printf("begin insert\n");
            double x = hpoint.pos.x;
            double y = hpoint.pos.y;
            double z = hpoint.pos.z;
            int ix, iy, iz;
            compute_coord(x,y,z,ix,iy,iz);
            //printf("pass test1\n");
            //if (hash(ix, iy, iz) >= hashsize || hash(ix, iy, iz) < 0)  { printf("bad index%d\n", ++debug_counter); printf("hash: %d\n", hash(ix, iy, iz));}
            hashtable[hash(ix,iy,iz)].push_back(hpoint);
            //printf("end insert\n");
        }
        vector<Hitpoint> * getVectorPtr(int ix, int iy, int iz) {
            return &hashtable[hash(ix, iy, iz)];
        }
        int getSize() {
            return this -> hashsize;
        }
        vector<Hitpoint> * getPtr() {
            return &(hashtable[0]);
        }
        // the data
        vector<vector<Hitpoint> > hashtable;
        int hashsize;
        int num_of_cell_per_dim;
        double celllength;
        int debug_counter;
};

#endif