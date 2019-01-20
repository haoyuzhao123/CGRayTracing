#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <vector>
#include <ctime>
#include <omp.h>

using namespace std;

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#include "headers/vec3.h"
#include "headers/util.h"
#include "headers/objects.h"
#include "headers/sampling.h"
#include "headers/hitpoints.h"
#include "headers/hash.h"
#include "headers/bezier.h"
#include "headers/texture.h"

const double eps = 1e-4;
const double INF = 1e10;
const double PI = 3.14159265358979;

const int width = 1024;
const int height = 768;
// image_data : the data for generating PNG
char image_data[width * height * 3];
// image : the image data in the computation
Vec3 image[height][width];

const int MAX_DEPTH = 5;
const double alpha = 0.7;

const Vec3 background = Vec3(0.0, 0.0, 0.0);

int debug_counter = 0;

void trace(const Vec3 &org, const Vec3 &dir, const vector<Object *> &objs, Vec3 flux, Vec3 adj, bool flag, int depth, Hashtable & htable, int x, int y) {
	// org is the origin of the ray
	// dir is the direction of the ray, should be normalized
	// depth is the current depth of the ray tracing algorithm
	if (depth >= MAX_DEPTH) {
		return;
	}
	// find the nearest intersection point
	double len;
	int id = -1;
	Vec3 normalvec;
	Vec3 temp;
	double nearest = INF;
	for (int i = 0; i < objs.size(); i++) {
		if (objs[i]->intersect(org, dir, len, temp)) {
			if (len < nearest) {
				id = i;
				nearest = len;
				normalvec = temp;
			}
		}
	}
	if (id == -1) { // no intersection with the objects
		return;
	}
	const Object * obj = objs[id];
	Vec3 intersection = org + dir * nearest;
	//Vec3 normalvec = obj->normalvec(intersection);
	// inside the object, change the direction of the normal vector
	bool into = true;
	Vec3 normalvec_old = normalvec;
	if (normalvec.dot(dir) > 0) {
		normalvec = -normalvec;
		into = false;
	}
	Vec3 f = obj->getSurfaceColor(intersection);

	double p = max(f.x, f.y, f.z);


	if (obj->getReflection() < eps && obj->getTransparency() < eps) {
		// diffusion
		double r = 200.0 / height;
		if (flag) {
			// the ray from the eye
			Hitpoint hp = Hitpoint();
			hp.f = f * adj;
			hp.pos = intersection;
			hp.normal = normalvec;
			hp.w = x;
			hp.h = y;
			hp.flux = Vec3();
			hp.r2 = r * r;
			hp.n = 0;
			//hitpoints.push_back(hp);
			//printf("begin insert\n");
			htable.insert(hp);
			//printf("end insert\n");
		}
		else {
			// the photon ray from the light source
			int ix, iy, iz;
			htable.compute_coord(intersection.x, intersection.y, intersection.z, ix, iy, iz);
			ix -= 1;
			iy -= 1;
			iz -= 1;
			int idx = 0, idy = 0, idz = 0;
			// search for the adjacent grids
			for (idx = 0; idx < 3; idx++)
			for (idy = 0; idy < 3; idy++)
			for (idz = 0; idz < 3; idz++) {
				int hashid = htable.hash(ix + idx, iy + idy, iz + idz);
				for (int i = 0; i < htable.hashtable[hashid].size(); i++) {
					Vec3 d = htable.hashtable[hashid][i].pos - intersection;
					if ((htable.hashtable[hashid][i].normal.dot(normalvec) > eps) && (d.dot(d) <= htable.hashtable[hashid][i].r2)) {
						// the equations here comes from the equations in smallppm.cpp
						// in smallppm.cpp, n = N / alpha
						double g = (htable.hashtable[hashid][i].n * alpha + alpha) /		(htable.hashtable[hashid][i].n * alpha + 1.0);
						htable.hashtable[hashid][i].r2 *= g;
						htable.hashtable[hashid][i].n++;
						htable.hashtable[hashid][i].flux = (htable.hashtable[hashid][i].flux + htable.hashtable[hashid][i].f.mul(flux) * (1.0 / PI)) * g;
					}
				}
			}
			Vec3 newdir = uniform_sampling_halfsphere(normalvec);
			trace(intersection, newdir, objs, f * flux * (1.0 / p), adj, flag, depth+1, htable, x, y);
		}
	} else if (obj -> getTransparency() < eps) {
		// mirror
		Vec3 newdir = dir - normalvec * 2.0 * normalvec.dot(dir);
		double refl = obj -> getReflection();
		intersection = intersection + normalvec * eps; // prevent problem created by double number precision
		trace(intersection, newdir, objs, f * flux * refl, f * adj * refl, flag, depth+1, htable, x, y);
	} else {
		
		// refraction
		//Ray lr(x,r.d-n*2.0*n.dot(r.d)); 
		//bool into = (normalvec_old.dot(nl)>0.0);
		double nc = 1.0, nt=1.33, nnt = into?nc/nt:nt/nc, ddn = dir.dot(normalvec), cos2t;
		Vec3 refl_dir = dir - normalvec_old * 2.0 * normalvec_old.dot(dir);

		// total internal reflection
		if ((cos2t=1-nnt*nnt*(1-ddn*ddn))<0) return trace(intersection + normalvec * eps, refl_dir, objs, flux, adj, flag, depth+1, htable, x, y);

		Vec3 refr_dir = (dir * nnt - normalvec_old * ((into?1:-1)*(ddn*nnt+sqrt(cos2t)))).normalize();
		if (!into && refr_dir.dot(normalvec_old) < 0) {
			fprintf(stderr, "error\n");
		}
		double a=nt-nc, b=nt+nc, R0=a*a/(b*b), c = 1-(into?-ddn:refr_dir.dot(normalvec_old));
		double Re=R0+(1-R0)*c*c*c*c*c,P=Re;
		
		Vec3 fa=f.mul(adj);
		if (flag) {
			// eye ray (trace both rays)
			trace(intersection + normalvec * eps, refl_dir, objs, flux, fa * Re , flag, depth+1, htable, x, y);
			trace(intersection - normalvec * eps, refr_dir, objs, flux, fa * (1 - Re), flag, depth+1, htable, x, y);
		} else {
			// photon ray (pick one via Russian roulette)
			if (uniform_sampling_zeroone() < 0.5) {
				trace(intersection + normalvec * eps, refl_dir, objs, flux, fa * Re * 0.3, flag, depth+1, htable, x, y);
			} else {
				trace(intersection - normalvec * eps, refr_dir, objs, flux, fa * (1 - Re * 0.3), flag, depth+1, htable, x, y);
			}
		}
	}
}

void render(const vector<Object *> &objs) {
	// the axis: x, left to right(width)
	//           y, bottom to top(height)
	//           z, satisfies the right hand rule
	// the cam should be at (0,0,-10)
	// the image has z = 0
	// the x axis of the image range from (-10,10)
	// in this project, we assume that there is only 1 point light source
	int num_of_samples = 1;
	double focus_plane = 20.0;
	double radius = 1.5;
	Vec3 lightorg = Vec3(0,19.999,20);
	Vec3 camorg = Vec3(0,0,-10);
	//vector<Hitpoint> hitpoints;
	double r = 200.0 / height;
	Hashtable htable = Hashtable(1000001,r);
	for (int h = 0; h < height; h++) {
		fprintf(stderr, "\rHitPointPass %5.2f%%", 100.0 * (h+1) / height);
		for (int w = 0; w < width; w++) {
			double x = (2.0 * ((double)w/width)-1) * 10.0;
			double y = (2.0 * ((double)h/height)-1) * 10.0 * height / width;
			double x1 = (2.0 * ((double)(w+0.5)/width)-1) * 10.0;
			double y1 = (2.0 * ((double)h/height)-1) * 10.0 * height / width;
			double x2 = (2.0 * ((double)(w-0.5)/width)-1) * 10.0;
			double y2 = (2.0 * ((double)h/height)-1) * 10.0 * height / width;
			double x3 = (2.0 * ((double)w/width)-1) * 10.0;
			double y3 = (2.0 * ((double)(h+0.5)/height)-1) * 10.0 * height / width;
			double x4 = (2.0 * ((double)w/width)-1) * 10.0;
			double y4 = (2.0 * ((double)(h-0.5)/height)-1) * 10.0 * height / width;
			Vec3 dir = (Vec3(x,y,0) - camorg).normalize();
			Vec3 dir1 = (Vec3(x1,y1,0) - camorg).normalize();
			Vec3 dir2 = (Vec3(x2,y2,0) - camorg).normalize();
			Vec3 dir3 = (Vec3(x3,y3,0) - camorg).normalize();
			Vec3 dir4 = (Vec3(x4,y4,0) - camorg).normalize();
			Vec3 point_on_focus = dir * ((focus_plane - camorg.z) / dir.z) + camorg;
			for (int j = 0; j < num_of_samples; j++) {
				Vec3 neworg = camorg + uniform_sampling_circle(radius);
				Vec3 newdir = (point_on_focus - neworg).normalize();
				//trace(neworg, newdir, objs, Vec3(), Vec3(1,1,1), true, 0, htable, w, h);
				
				trace(camorg, dir, objs, Vec3(), Vec3(1,1,1), true, 0, htable, w, h);
				/*
				trace(camorg, dir1, objs, Vec3(), Vec3(1,1,1), true, 0, htable, w, h);
				trace(camorg, dir2, objs, Vec3(), Vec3(1,1,1), true, 0, htable, w, h);
				trace(camorg, dir3, objs, Vec3(), Vec3(1,1,1), true, 0, htable, w, h);
				trace(camorg, dir4, objs, Vec3(), Vec3(1,1,1), true, 0, htable, w, h);
				 */
			}
			//trace(camorg, dir, objs, Vec3(), Vec3(1,1,1), true, 0, htable, w, h);
		}
	}
	fprintf(stderr,"\n");
	// multi-thread for the photon pass
	// total photon = num_photon * num_threads
	int num_photon = 2560000;
	int num_threads = 8;
	omp_set_num_threads(num_threads);
	#pragma omp parallel 
	{
		// need to set the rand() seed in each thread
		srand(int(time(NULL)) ^ omp_get_thread_num());
		#pragma omp parallel for
		for (int i = 0; i < num_photon; i++) {
			
			if(omp_get_thread_num() == 0 && i % 10000 == 0) {
				fprintf(stderr, "\rPhotonPass %5.2f%%", 100.0 * i / num_photon);
			}
			
			//double p = 100.0 * (i+1) / num_photon;
			//fprintf(stderr, "\rPhotonPass %5.2f%%",p);
			//for(int j = 0; j < 1000; j++) {
				double a = uniform_sampling_zeroone() * 4 - 2;
				double b = uniform_sampling_zeroone() * 4 - 2;
				Vec3 disturbance = Vec3(a, 0, b);
				//Vec3 disturbance = Vec3(0,0,0);
				//Vec3 dir = uniform_sampling_halfsphere(Vec3(0,-1,0));
				Vec3 dir = uniform_sampling_sphere();
				trace(lightorg + disturbance, dir, objs, Vec3(700,700,700)*(PI*4.0), Vec3(1,1,1), false, 	0, htable, 0, 0);
			//}
		}
	}

	//vector<Hitpoint> * ptr = htable.getPtr();
	for (int j = 0; j < htable.hashtable.size(); j++) {
		//vector<Hitpoint> * nptr = ptr + j;
		for (int i = 0; i < htable.hashtable[j].size(); i++) {
			Hitpoint hp = htable.hashtable[j][i];
			image[hp.h][hp.w] = image[hp.h][hp.w] + hp.flux * (1.0/(PI*hp.r2*num_photon*num_threads*num_of_samples));
		}
	}

	int debug_counter = 0;
	for (int i = 0; i < htable.hashtable.size(); i++) {
		debug_counter += htable.hashtable[i].size();
	}

	printf("\nhitpoints: %d\n", debug_counter);
}

int main(int argc, char *argv[]) {
	srand(time(0));
	//initialize the image data
	for (int i = 0; i < height; i++) {
		for (int j = 0; j < width; j++) {
			image[i][j] = Vec3(0,0,0);
		}
	}

	vector<Object *> objs;
	vector<Sphere> sphs;
	vector<Plane> plns;
	/*
	sphs.push_back(Sphere(Vec3(0.0, -10020, 0), 10000, Vec3(0.25, 0.25, 0.25), 0.0, 0.0));
	sphs.push_back(Sphere(Vec3(10020, 0.0, 0), 10000, Vec3(0.25, 0.75, 0.25), 0.0, 0.0));
	sphs.push_back(Sphere(Vec3(-10020, 0.0, 0), 10000, Vec3(0.75, 0.25, 0.25), 0.0, 0.0));
	sphs.push_back(Sphere(Vec3(0.0, 0.0, 10040), 10000, Vec3(0.25, 0.25, 0.25), 0.0, 0.0));
	sphs.push_back(Sphere(Vec3(0.0, 10020, 0), 10000, Vec3(0.25, 0.25, 0.25), 0.0, 0.0));
	*/
	//sphs.push_back(Sphere(Vec3(0.0, 0.0, -10015), 10000, Vec3(0, 0, 0), 0.0, 0.0));
	//sphs.push_back(Sphere(Vec3(-15.0, -20.0, 60), 10, Vec3(0.3, 0.3, 0.3), 0.0, 0.0));
	//sphs.push_back(Sphere(Vec3(10.0, -20.0, 60), 7, Vec3(1.0, 1.0, 1.0), 0.8, 0.0));
	//sphs.push_back(Sphere(Vec3(10.0, -20.0, 30), 7, Vec3(1.0, 1.0, 1.0), 0.8, 0.5));

	TriangleMesh tm1("model/dragon.txt", 1.5, Vec3(-5, -20, 30), Vec3(0.25, 0.25, 0.5), 0.0, 0.0, 1);
	//TriangleMesh tm("model/lowpolybunny.txt", 10, Vec3(0, -15, 40), Vec3(1.0, 1.0, 1.0), 0.8, 0.5);
	//TriangleMesh tm2("model/Mesh001.obj", 20, Vec3(0, -15, 30), Vec3(1.0, 1.0, 1.0), 0.8, 0.5, 2);
	//TriangleMesh tm2("model/water.txt", 7, Vec3(-20, -10, 40), Vec3(1.0, 1.0, 1.0), 0.8, 0.5, 2);

	//TriangleMesh tm2("model/dragon.txt", 1.5, Vec3(-10, -20, 30), Vec3(0.25, 0.25, 0.5), 0.0, 0.0, 1);

	int w, h, bpp;
    unsigned char * texture = stbi_load("texture/stone.jpg", &w, &h, &bpp, 3);
	vector<vector<Vec3> > tdata;
	int ctr = 0;
	for(int i = 0; i < h; i++) {
		vector<Vec3> v;
		for(int j = 0; j < w; j++) {
			Vec3 col = Vec3();
			col.x = (double)texture[ctr] / (double)256;
			ctr++;
			col.y = (double)texture[ctr] / (double)256;
			ctr++;
			col.z = (double)texture[ctr] / (double)256;
			ctr++;
			v.push_back(col);
		}
		tdata.push_back(v);
	}

	stbi_image_free(texture);

	Texture tex = Texture(tdata, Vec3(0,1,0), Vec3(-21, 0, 0), 42, 40, true);
	/*
	texture = stbi_load("texture/iiis.png", &w, &h, &bpp, 3);
	vector<vector<Vec3> > tdata2;
	ctr = 0;
	for(int i = 0; i < h; i++) {
		vector<Vec3> v;
		for(int j = 0; j < w; j++) {
			Vec3 col = Vec3();
			col.x = (double)texture[ctr] / (double)256;
			ctr++;
			col.y = (double)texture[ctr] / (double)256;
			ctr++;
			col.z = (double)texture[ctr] / (double)256;
			ctr++;
			v.push_back(col);
		}
		tdata2.push_back(v);
	}

	Texture tex2 = Texture(tdata2, Vec3(0,0,-1), Vec3(-10, -10, 40), 20, 10);

	stbi_image_free(texture);
	 */

	printf("testtesttest\n");

	//vector<Triangle> tris;
	plns.push_back(Plane(Vec3(0.0, -20, 0), Vec3(0,1,0), Vec3(0.15, 0.15, 0.15), 0.0, 0.0, tex));
	plns.push_back(Plane(Vec3(20, 0.0, 0), Vec3(-1,0,0), Vec3(0.15, 0.50, 0.15), 0.0, 0.0));
	plns.push_back(Plane(Vec3(-20, 0.0, 0), Vec3(1,0,0), Vec3(0.50, 0.15, 0.15), 0.0, 0.0));
	plns.push_back(Plane(Vec3(0.0, 0.0, 40), Vec3(0,0,-1), Vec3(0.15, 0.15, 0.15), 0.0, 0.0));
	//plns.push_back(Plane(Vec3(0.0, 0.0, 40), Vec3(0,0,-1), Vec3(0.15, 0.15, 0.15), 0.0, 0.0));
	plns.push_back(Plane(Vec3(0.0, 20, 0), Vec3(0,-1,0), Vec3(0.15, 0.15, 0.15), 0.0, 0.0));
	
	Object * obj;
	for (int i = 0; i < sphs.size(); i++) {
		obj = &sphs[i];
		objs.push_back(obj);
	}
	for (int i = 0; i < plns.size(); i++) {
		obj = &plns[i];
		objs.push_back(obj);
	}
	
	obj = &tm1;
	objs.push_back(obj);
	//obj = &tm2;
	//objs.push_back(obj);
	
	/*
	vector<Vec3> cp;
	cp.push_back(Vec3(0,-10,4));
	cp.push_back(Vec3(0,2,4));
	cp.push_back(Vec3(0,-2,0));
	cp.push_back(Vec3(0,10,2));
	Bezier b = Bezier(cp, Vec3(15,-10.1,35), Vec3(1.0, 1.0, 1.0), 0.5, 0.0);
	obj = &b;
	objs.push_back(obj);
	*/
	/*
	b.gradP(0.5).print();
	Vec3 org = Vec3(0,0,0);
	Vec3 dir = Vec3(0.1,0.1,1);
	Vec3 initial;
	initial.x = 22;
	initial.y = 0.6;
	initial.z = 2;
	initial.print();
	b.newtonMethod(initial, org, dir);
	*/
	/*
	b.funcValue(Vec3(0,0.5,0), Vec3(0,0,0), Vec3(0,0,0)).print();
	Vec3 resa, resb, resc;
	b.gradValue(Vec3(0.5,0.5,0.5), Vec3(0,0,0), Vec3(0,0,0), resa, resb, resc);
	resa.print();
	resb.print();
	resc.print();
	*/

	// render
	render(objs);
	// generate PNG
	int counter = 0;
	for (int i = 0; i < height; i++) {
		for (int j = 0; j < width; j++) {
			image_data[3 * counter] = gammaCorr(image[height - i - 1][j].x);
        	image_data[3 * counter + 1] = gammaCorr(image[height - i - 1][j].y);
        	image_data[3 * counter + 2] = gammaCorr(image[height - i - 1][j].z);
			counter++;
		}
	}
    stbi_write_png("test.png", width, height, 3, image_data, width * 3);
    return 0; 
}
