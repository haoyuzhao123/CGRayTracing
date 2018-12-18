// expanded smallppm (code is exactly the same as smallppm.cpp but with more comments)

#include <math.h>   // smallppm, Progressive Photon Mapping by T. Hachisuka
#include <stdlib.h> // originally smallpt, a path tracer by Kevin Beason, 2008
#include <stdio.h>  // Usage: ./smallppm 100000 && xv image.ppm
#define PI ((double)3.14159265358979) // ^^^^^^:number of photons emitted
#define ALPHA ((double)0.7) // the alpha parameter of PPM

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#include "headers/vec3.h"
#include "headers/ray.h"
#include "headers/util.h"

const int w=1024;
const int h=768;
char image_data[w * h * 3];

// Halton sequence with reverse permutation
int primes[61]={
	2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,
	83,89,97,101,103,107,109,113,127,131,137,139,149,151,157,163,167,173,179,181,
	191,193,197,199,211,223,227,229,233,239,241,251,257,263,269,271,277,281,283
};
inline int rev(const int i,const int p) {
	if (i==0) return i; else return p-i;
}
double hal(const int b, int j) {
	const int p = primes[b]; 
	double h = 0.0, f = 1.0 / (double)p, fct = f;
	while (j > 0) {
		h += rev(j % p, p) * fct; j /= p; fct *= f;
	}
	return h;
}

struct AABB {Vec3 min, max; // axis aligned bounding box
	inline void fit(const Vec3 &p)
	{
		if (p.x<min.x)min.x=p.x; // min
		if (p.y<min.y)min.y=p.y; // min
		if (p.z<min.z)min.z=p.z; // min
		max.x=MAX(p.x, max.x);
		max.y=MAX(p.y, max.y);
		max.z=MAX(p.z, max.z);
	}
	inline void reset() {
		min=Vec3(1e20,1e20,1e20); 
		max=Vec3(-1e20,-1e20,-1e20);
	}
};

struct HPoint {
	Vec3 f,pos,nrm,flux; 
	double r2; 
	unsigned int n; // n = N / ALPHA in the paper
	int pix;
};

struct List {HPoint *id; List *next;};
List* ListAdd(HPoint *i,List* h){
	List* p=new List;
	p->id=i;
	p->next=h;
	return p;
}

unsigned int num_hash, pixel_index, num_photon;
double hash_s; List **hash_grid; List *hitpoints = NULL; AABB hpbbox;

// spatial hash function
inline unsigned int hash(const int ix, const int iy, const int iz) {
	return (unsigned int)((ix*73856093)^(iy*19349663)^(iz*83492791))%num_hash;
}

void build_hash_grid(const int w, const int h) {
	// find the bounding box of all the measurement points
	hpbbox.reset();
	List *lst = hitpoints;
	while (lst != NULL) {
		HPoint *hp=lst->id; 
		lst=lst->next; 
		hpbbox.fit(hp->pos);
	}

	// heuristic for initial radius
	Vec3 ssize = hpbbox.max - hpbbox.min;
	double irad = ((ssize.x + ssize.y + ssize.z) / 3.0) / ((w + h) / 2.0) * 2.0;

	// determine hash table size
	// we now find the bounding box of all the measurement points inflated by the initial radius
	hpbbox.reset(); 
	lst = hitpoints; 
	int vphoton = 0; 
	while (lst != NULL) {
		HPoint *hp = lst->id; 
		lst = lst->next;
		hp->r2=irad *irad; 
		hp->n = 0; 
		hp->flux = Vec3();
		vphoton++; 
		hpbbox.fit(hp->pos-irad); 
		hpbbox.fit(hp->pos+irad);
	}

	// make each grid cell two times larger than the initial radius
	hash_s=1.0/(irad*2.0); 
	num_hash = vphoton; 

	// build the hash table
	hash_grid=new List*[num_hash];
	for (unsigned int i=0; i<num_hash;i++) hash_grid[i] = NULL;
	lst = hitpoints; 
	while (lst != NULL) 
	{ 
		HPoint *hp = lst->id; 
		lst = lst->next;
		Vec3 BMin = ((hp->pos - irad) - hpbbox.min) * hash_s;
		Vec3 BMax = ((hp->pos + irad) - hpbbox.min) * hash_s;
		for (int iz = abs(int(BMin.z)); iz <= abs(int(BMax.z)); iz++)
		{
			for (int iy = abs(int(BMin.y)); iy <= abs(int(BMax.y)); iy++)
			{
				for (int ix = abs(int(BMin.x)); ix <= abs(int(BMax.x)); ix++)
				{
					int hv=hash(ix,iy,iz); 
					hash_grid[hv]=ListAdd(hp,hash_grid[hv]);
				}
			}
		}
	}
}

enum Refl_t {DIFF, SPEC, REFR};  // material types, used in radiance()
struct Sphere {double rad; Vec3 p, c; Refl_t refl;
	Sphere(double r_,Vec3 p_,Vec3 c_,Refl_t re_) : rad(r_),p(p_),c(c_),refl(re_){}
	inline double intersect(const Ray &r) const { 
		// ray-sphere intersection returns distance
		Vec3 op=p-r.org; 
		double t, b=op.dot(r.dir), det=b*b-op.dot(op)+rad*rad;
		if (det < 0) {
			return 1e20; 
		}
		else {
			det = sqrt(det);
		}
		return (t=b-det) > 1e-4 ? t : ((t=b+det)>1e-4 ? t : 1e20);
	}
};

Sphere sph[] = { // Scene: radius, position, color, material
  Sphere(1e5, Vec3( 1e5+1,40.8,81.6), Vec3(.75,.25,.25),DIFF),//Left
  Sphere(1e5, Vec3(-1e5+99,40.8,81.6),Vec3(.25,.25,.75),DIFF),//Right
  Sphere(1e5, Vec3(50,40.8, 1e5),     Vec3(.75,.75,.75),DIFF),//Back
  Sphere(1e5, Vec3(50,40.8,-1e5+170), Vec3(),           DIFF),//Front
  Sphere(1e5, Vec3(50, 1e5, 81.6),    Vec3(.75,.75,.75),DIFF),//Bottomm
  Sphere(1e5, Vec3(50,-1e5+81.6,81.6),Vec3(.75,.75,.75),DIFF),//Top
  Sphere(16.5,Vec3(27,16.5,47),       Vec3(1,1,1)*.999, SPEC),//Mirror
  Sphere(16.5,Vec3(73,16.5,88),       Vec3(1,1,1)*.999, REFR),//Glass
  Sphere(8.5, Vec3(50,8.5,60),        Vec3(1,1,1)*.999, DIFF)};//Middle

// find the closet interection
inline bool intersect(const Ray &r,double &t,int &id){
	int n = sizeof(sph) / sizeof(Sphere); 
	double d, inf = 1e20; t = inf;
	for(int i=0;i<n;i++){
		d=sph[i].intersect(r);
		if(d<t){
			t=d;
			id=i;
		}
	}
	return t<inf;
}

// generate a photon ray from the point light source with QMC
void genp(Ray* pr, Vec3* f, int i) {
	*f = Vec3(2500,2500,2500)*(PI*4.0); // flux
	double p=2.*PI*hal(0,i),t=2.*acos(sqrt(1.-hal(1,i)));
	double st=sin(t);
	pr->dir=Vec3(cos(p)*st,cos(t),sin(p)*st);
	pr->org=Vec3(50,60,85);
}

void trace(const Ray &r,int dpt,bool m,const Vec3 &fl,const Vec3 &adj,int i) 
{
	double t;
	int id; 

	dpt++;
	if(!intersect(r,t,id)||(dpt>=20))return;
 
	int d3=dpt*3;
	const Sphere &obj = sph[id]; 
	Vec3 x=r.org+r.dir*t, n=(x-obj.p).normalize(), f=obj.c;
	Vec3 nl=n.dot(r.dir)<0?n:n*-1; 
	double p=f.x>f.y&&f.x>f.z?f.x:f.y>f.z?f.y:f.z;

	if (obj.refl == DIFF) {
		// Lambertian

		// use QMC to sample the next direction
		double r1=2.*PI*hal(d3-1,i),r2=hal(d3+0,i);
		double r2s=sqrt(r2);
		Vec3 w=nl,u=((fabs(w.x)>.1?Vec3(0,1):Vec3(1)).cross(w)).normalize();
		Vec3 v=w.cross(u), d = (u*cos(r1)*r2s + v*sin(r1)*r2s + w*sqrt(1-r2)).normalize();

		if (m) {
			// eye ray
			// store the measurment point
			HPoint* hp=new HPoint; 
			hp->f=f.mul(adj); 
			hp->pos=x;
			hp->nrm=n; 
			hp->pix = pixel_index; 
			hitpoints = ListAdd(hp, hitpoints);
		} 
		else 
		{
			// photon ray
			// find neighboring measurement points and accumulate flux via progressive density estimation
			Vec3 hh = (x-hpbbox.min) * hash_s;
			int ix = abs(int(hh.x)), iy = abs(int(hh.y)), iz = abs(int(hh.z));
			// strictly speaking, we should use #pragma omp critical here.
			// it usually works without an artifact due to the fact that photons are 
			// rarely accumulated to the same measurement points at the same time (especially with QMC).
			// it is also significantly faster.
			{
				List* hp = hash_grid[hash(ix, iy, iz)]; 
				while (hp != NULL) {
					HPoint *hitpoint = hp->id; 
					hp = hp->next; 
					Vec3 v = hitpoint->pos - x;
					// check normals to be closer than 90 degree (avoids some edge brightning)
					if ((hitpoint->nrm.dot(n) > 1e-3) && (v.dot(v) <= hitpoint->r2)) {
						// unlike N in the paper, hitpoint->n stores "N / ALPHA" to make it an integer value
						double g = (hitpoint->n*ALPHA+ALPHA) / (hitpoint->n*ALPHA+1.0);
						hitpoint->r2=hitpoint->r2*g; 
						hitpoint->n++;
						hitpoint->flux=(hitpoint->flux+hitpoint->f.mul(fl)*(1./PI))*g;
					}
				}
			}
			if (hal(d3+1,i)<p) trace(Ray(x,d),dpt,m,f.mul(fl)*(1./p),adj,i);
		}

	} else if (obj.refl == SPEC) {
		// mirror
		trace(Ray(x, r.dir-n*2.0*n.dot(r.dir)), dpt, m, f.mul(fl), f.mul(adj),i);

	} else {
		// glass
		Ray lr(x,r.dir-n*2.0*n.dot(r.dir)); 
		bool into = (n.dot(nl)>0.0);
		double nc = 1.0, nt=1.5, nnt = into?nc/nt:nt/nc, ddn = r.dir.dot(nl), cos2t;

		// total internal reflection
		if ((cos2t=1-nnt*nnt*(1-ddn*ddn))<0) return trace(lr,dpt,m,fl,adj,i);

		Vec3 td = (r.dir*nnt - n*((into?1:-1)*(ddn*nnt+sqrt(cos2t)))).normalize();
		double a=nt-nc, b=nt+nc, R0=a*a/(b*b), c = 1-(into?-ddn:td.dot(n));
		double Re=R0+(1-R0)*c*c*c*c*c,P=Re;Ray rr(x,td);Vec3 fa=f.mul(adj);
		if (m) {
			// eye ray (trace both rays)
			trace(lr,dpt,m,fl,fa*Re,i);
			trace(rr,dpt,m,fl,fa*(1.0-Re),i);
		} else {
			// photon ray (pick one via Russian roulette)
			(hal(d3-1,i)<P)?trace(lr,dpt,m,fl,fa,i):trace(rr,dpt,m,fl,fa,i);
		}
	}
}

int main(int argc, char *argv[]) {
	// samps * 1000 photon paths will be traced
	int samps = (argc==2) ? MAX(atoi(argv[1])/1000,1) : 1000;

	// trace eye rays and store measurement points
	Ray cam(Vec3(50,48,295.6), Vec3(0,-0.042612,-1).normalize());
	Vec3 cx=Vec3(w*.5135/h), cy=(cx.cross(cam.dir)).normalize()*.5135, *c=new Vec3[w*h], vw;
	for (int y=0; y<h; y++){
		fprintf(stderr, "\rHitPointPass %5.2f%%", 100.0*y/(h-1));
		for (int x=0; x<w; x++) {
			pixel_index = x + y * w;
			Vec3 d = cx * ((x + 0.5) / w - 0.5) + cy * (-(y + 0.5) / h + 0.5)+cam.dir;
			trace(Ray(cam.org + d * 140, d.normalize()), 0, true, Vec3(), Vec3(1, 1, 1),0);
		}
	}
	fprintf(stderr,"\n"); 
	
	// build the hash table over the measurement points
	build_hash_grid(w,h); 
	
	// trace photon rays with multi-threading
	num_photon=samps; 
	vw=Vec3(1,1,1);
	#pragma omp parallel for schedule(dynamic, 1)
	for(int i=0;i<num_photon;i++) {
		double p=100.*(i+1)/num_photon;
		fprintf(stderr,"\rPhotonPass %5.2f%%",p); 
		int m=1000*i; 
		Ray r; 
		Vec3 f;
		for(int j=0;j<1000;j++){
			genp(&r,&f,m+j); 
			trace(r,0,0>1,f,vw,m+j);
		}
	}

	// density estimation
	List* lst=hitpoints; 
	while (lst != NULL) {
		HPoint* hp=lst->id;
		lst=lst->next;
		int i=hp->pix;
		c[i]=c[i]+hp->flux*(1.0/(PI*hp->r2*num_photon*1000.0));
	}

	for(int i = 0; i< w * h; i++) {
        image_data[3*i] = gammaCorr(c[i].x);
        image_data[3*i+1] = gammaCorr(c[i].y);
        image_data[3*i+2] = gammaCorr(c[i].z); 
	}
    stbi_write_png("test.png", w, h, 3, image_data, w * 3);
    return 0;
}
