#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <vector>

using namespace std;

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#include "headers/vec3.h"
#include "headers/util.h"
#include "headers/objects.h"
#include "headers/sampling.h"
#include "headers/hitpoints.h"

const double eps = 1e-4;
const double INF = 1e10;
const double PI = 3.14159265358979;

const int width = 128;
const int height = 96;
// image_data : the data for generating PNG
char image_data[width * height * 3];
// image : the image data in the computation
Vec3 image[height][width];

const int MAX_DEPTH = 5;
const double alpha = 0.7;

const Vec3 background = Vec3(0.0, 0.0, 0.0);

int debug_counter = 0;

void trace(const Vec3 &org, const Vec3 &dir, const vector<Object *> &objs, Vec3 flux, Vec3 adj, bool flag, int depth, vector<Hitpoint> &hitpoints, int x, int y, bool debug = false) {
	// org is the origin of the ray
	// dir is the direction of the ray, should be normalized
	// depth is the current depth of the ray tracing algorithm
	if(debug) printf("depth: %d\n", depth);
	if (depth >= MAX_DEPTH) {
		return;
	}
	if(debug) printf("test1\n");
	// find the nearest intersection point
	double len;
	int id = -1;
	double nearest = INF;
	for (int i = 0; i < objs.size(); i++) {
		if (objs[i]->intersect(org, dir, len)) {
			if (len < nearest) {
				id = i;
				nearest = len;
			}
		}
	}
	if(debug) printf("test2\n");
	if (id == -1) { // no intersection with the objects
		return;
	}
	const Object * obj = objs[id];
	Vec3 intersection = org + dir * nearest;
	Vec3 normalvec = obj->normalvec(intersection);
	// inside the object, change the direction of the normal vector
	if (normalvec.dot(dir) > 0) {
		normalvec = -normalvec;
	}
	Vec3 f = obj->getSurfaceColor();

	double p = max(f.x, f.y, f.z);

	if(debug) printf("test3\n");

	if (obj->getReflection() < eps && obj->getTransparency() < eps) {
		if(debug) printf("test4\n");
		// diffusion
		if (flag) {
			// the ray from the eye
			Hitpoint hp = Hitpoint();
			hp.f = f * adj;
			hp.pos = intersection;
			hp.normal = normalvec;
			hp.w = x;
			hp.h = y;
			hp.flux = Vec3();
			hp.r2 = 4 * (100.0 / height) * (100.0 / height);
			hp.n = 0;
			hitpoints.push_back(hp);
		}
		else {
			if(debug) printf("test5\n");
			// the photon ray from the light source
			for (int i = 0; i < hitpoints.size(); i++) {
				Vec3 d = hitpoints[i].pos - intersection;
				if ((hitpoints[i].normal.dot(normalvec) > eps) && (d.dot(d) <= hitpoints[i].r2)) {
					debug_counter++;
					double g = (hitpoints[i].n * alpha + alpha) / (hitpoints[i].n * alpha + 1.0);
					hitpoints[i].r2 *= g;
					hitpoints[i].n++;
					hitpoints[i].flux = (hitpoints[i].flux + hitpoints[i].f.mul(flux) * (1.0 / PI)) * g;
				}
			}
			if(debug) printf("test6\n");
			//sample the next direction
			if(debug) printf("%f\n", normalvec.norm());
			if(debug) printf("%f %f %f\n", normalvec.x, normalvec.y, normalvec.z);
			Vec3 newdir = uniform_sampling_halfsphere(normalvec);
			//printf("recursice call\n");
			trace(intersection, newdir, objs, f * flux * (1.0 / p), adj, flag, depth+1, hitpoints, x, y, debug);
			//printf("end recursice call\n");
		}
	}
	if(debug) printf("out\n");
}

void render(const vector<Object *> &objs) {
	// the axis: x, left to right(width)
	//           y, bottom to top(height)
	//           z, satisfies the right hand rule
	// the cam should be at (0,0,-10)
	// the image has z = 0
	// the x axis of the image range from (-10,10)
	// in this project, we assume that there is only 1 point light source
	// the light source locate that (0,50,20)
	Vec3 lightorg = Vec3(0,20,20);
	Vec3 camorg = Vec3(0,0,-10);
	vector<Hitpoint> hitpoints;
	for (int h = 0; h < height; h++) {
		fprintf(stderr, "\rHitPointPass %5.2f%%", 100.0 * h / height);
		for (int w = 0; w < width; w++) {
			double x = (2.0 * ((double)w/width)-1) * 10.0;
			double y = (2.0 * ((double)h/height)-1) * 10.0 * height / width;
			Vec3 dir = (Vec3(x,y,0) - camorg).normalize();
			trace(camorg, dir, objs, Vec3(), Vec3(1,1,1), true, 0, hitpoints, w, h);
		}
	}
	fprintf(stderr,"\n"); 
	int num_photon = 100000;

	for (int i = 0; i < num_photon; i++) {
		double p = 100.0 * (i+1) / num_photon;
		fprintf(stderr, "\rPhotonPass %5.2f%%",p);
		Vec3 dir = uniform_sampling_sphere();
		trace(lightorg, dir, objs, Vec3(2500,2500,2500)*(PI*4.0), Vec3(1,1,1), false, 0, hitpoints, 0, 0, false);
	}

	for (int i = 0; i < hitpoints.size(); i++) {
		Hitpoint hp = hitpoints[i];
		image[hp.h][hp.w] = image[hp.h][hp.w] + hp.flux * (1.0/(PI*hp.r2*num_photon));
	}

	printf("\nhitpoints: %d\n", hitpoints.size());
}

int main(int argc, char *argv[]) {
	srand(19);
	//initialize the image data
	for (int i = 0; i < height; i++) {
		for (int j = 0; j < width; j++) {
			image[i][j] = Vec3(0,0,0);
		}
	}

	vector<Object *> objs;
	vector<Sphere> sphs;
	sphs.push_back(Sphere(Vec3(0.0, -10020, 0.0), 10000, Vec3(0.75, 0.25, 0.25), 0.0, 0.0));
	sphs.push_back(Sphere(Vec3(10020, 0.0, 0.0), 10000, Vec3(0.25, 0.25, 0.75), 0.0, 0.0));
	sphs.push_back(Sphere(Vec3(-10020, 0.0, 0.0), 10000, Vec3(0.75, 0.75, 0.75), 0.0, 0.0));
	sphs.push_back(Sphere(Vec3(0.0, 0.0, 10070), 10000, Vec3(0.57, 0.75, 0.75), 0.0, 0.0));
	sphs.push_back(Sphere(Vec3(0.0, 10030, 0.0), 10000, Vec3(0.5, 0.5, 0.5), 0.0, 0.0));
	

	Object * obj;
	for (int i = 0; i < sphs.size(); i++) {
		obj = &sphs[i];
		objs.push_back(obj);
	}

	// render
	render(objs);
	// generate PNG
	int counter = 0;
	for (int i = 0; i < height; i++) {
		for (int j = 0; j < width; j++) {
			image_data[3 * counter] = gammaCorr(image[height - i][j].x);
        	image_data[3 * counter + 1] = gammaCorr(image[height - i][j].y);
        	image_data[3 * counter + 2] = gammaCorr(image[height - i][j].z);
			counter++;
		}
	}
    stbi_write_png("test.png", width, height, 3, image_data, width * 3);
    return 0;
}
