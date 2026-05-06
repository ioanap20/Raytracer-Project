#define _CRT_SECURE_NO_WARNINGS 1
#include <vector>
#include <cmath>
#include <random>
#include <omp.h>
#include <map>
#include <string>
#include <fstream>
#include <iostream>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#ifndef M_PI
#define M_PI 3.14159265358979323856
#endif

static std::default_random_engine engine[32];
static std::uniform_real_distribution<double> uniform(0, 1);

double sqr(double x) { return x * x; };

class Vector {
public:
	explicit Vector(double x = 0, double y = 0, double z = 0) {
		data[0] = x;
		data[1] = y;
		data[2] = z;
	}
	double norm2() const {
		return data[0] * data[0] + data[1] * data[1] + data[2] * data[2];
	}
	double norm() const {
		return sqrt(norm2());
	}
	void normalize() {
		double n = norm();
		data[0] /= n;
		data[1] /= n;
		data[2] /= n;
	}
	double operator[](int i) const { return data[i]; };
	double& operator[](int i) { return data[i]; };
	double data[3];
};

Vector operator+(const Vector& a, const Vector& b) {
	return Vector(a[0] + b[0], a[1] + b[1], a[2] + b[2]);
}
Vector operator-(const Vector& a, const Vector& b) {
	return Vector(a[0] - b[0], a[1] - b[1], a[2] - b[2]);
}
Vector operator*(const double a, const Vector& b) {
	return Vector(a*b[0], a*b[1], a*b[2]);
}
Vector operator*(const Vector& a, const double b) {
	return Vector(a[0]*b, a[1]*b, a[2]*b);
}
Vector operator/(const Vector& a, const double b) {
	return Vector(a[0] / b, a[1] / b, a[2] / b);
}
double dot(const Vector& a, const Vector& b) {
	return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}
Vector cross(const Vector& a, const Vector& b) {
	return Vector(a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]);
}

class Ray {
public:
	Ray(const Vector& origin, const Vector& unit_direction) : O(origin), u(unit_direction) {};
	Vector O, u;
};

class Object {
public:
	Object(const Vector& albedo, bool mirror = false, bool transparent = false) : albedo(albedo), mirror(mirror), transparent(transparent) {};

	virtual bool intersect(const Ray& ray, Vector& P, double& t, Vector& N) const = 0;

	Vector albedo;
	bool mirror, transparent;
};

class Sphere : public Object {
public:
	Sphere(const Vector& center, double radius, const Vector& albedo, bool mirror = false, bool transparent = false) : ::Object(albedo, mirror, transparent), C(center), R(radius) {};

	// returns true iif there is an intersection between the ray and the sphere
	// if there is an intersection, also computes the point of intersection P, 
	// t>=0 the distance between the ray origin and P (i.e., the parameter along the ray)
	// and the unit normal N
	bool intersect(const Ray& ray, Vector& P, double &t, Vector& N) const {
		// TODO (lab 1) : compute the intersection (just true/false at the begining of lab 1, then P, t and N as well)
		double R_squared = sqr(R);
		double delta;	
		delta = dot(ray.u, ray.O - C) * dot(ray.u, ray.O - C) - ((ray.O-C).norm2() - R_squared);	
		if (delta < 0){
			return false;
		}
		double t1 = dot(ray.u, C - ray.O) + sqrt(delta);
		double t2 = dot(ray.u, C - ray.O) - sqrt(delta);
		bool ok = false;
		if(t2 >= 0){
			t = t2;
			ok=true;
		}
		else if(t1>=0){
				t = t1;
			}
		else{return false;}
		
		P = ray.O + t * ray.u;
		N = P-C;
		N.normalize();
		return true;

	}

	double R;
	Vector C;
};


// Class only used in labs 3 and 4 
class TriangleIndices {
public:
	TriangleIndices(int vtxi = -1, int vtxj = -1, int vtxk = -1, int ni = -1, int nj = -1, int nk = -1, int uvi = -1, int uvj = -1, int uvk = -1, int group = -1) {
		vtx[0] = vtxi; vtx[1] = vtxj; vtx[2] = vtxk;
		uv[0] = uvi; uv[1] = uvj; uv[2] = uvk;
		n[0] = ni; n[1] = nj; n[2] = nk;
		this->group = group;
	};
	int vtx[3]; // indices within the vertex coordinates array
	int uv[3];  // indices within the uv coordinates array
	int n[3];   // indices within the normals array
	int group;  // face group
};

// I will provide you with an obj mesh loader (labs 3 and 4)
class TriangleMesh : public Object {
public:
	TriangleMesh(const Vector& albedo, bool mirror = true, bool transparent = false) : ::Object(albedo, mirror, transparent) {};

	// first scale and then translate the current object
	void scale_translate(double s, const Vector& t) {
		for (int i = 0; i < vertices.size(); i++) {
			vertices[i] = vertices[i] * s + t;
		}
		buildbox();
	}


	// read an .obj file
	void readOBJ(const char* obj) {
		std::ifstream f(obj);
		if (!f) return;

		std::map<std::string, int> mtls;
		int curGroup = -1, maxGroup = -1;

		// OBJ indices are 1-based and can be negative (relative), this normalizes them
		auto resolveIdx = [](int i, int size) {
			return i < 0 ? size + i : i - 1;
		};

		auto setFaceVerts = [&](TriangleIndices& t, int i0, int i1, int i2) {
			t.vtx[0] = resolveIdx(i0, vertices.size());
			t.vtx[1] = resolveIdx(i1, vertices.size());
			t.vtx[2] = resolveIdx(i2, vertices.size());
		};
		auto setFaceUVs = [&](TriangleIndices& t, int j0, int j1, int j2) {
			t.uv[0] = resolveIdx(j0, uvs.size());
			t.uv[1] = resolveIdx(j1, uvs.size());
			t.uv[2] = resolveIdx(j2, uvs.size());
		};
		auto setFaceNormals = [&](TriangleIndices& t, int k0, int k1, int k2) {
			t.n[0] = resolveIdx(k0, normals.size());
			t.n[1] = resolveIdx(k1, normals.size());
			t.n[2] = resolveIdx(k2, normals.size());
		};

		std::string line;
		while (std::getline(f, line)) {
			// Trim trailing whitespace
			line.erase(line.find_last_not_of(" \r\t\n") + 1);
			if (line.empty()) continue;

			const char* s = line.c_str();

			if (line.rfind("usemtl ", 0) == 0) {
				std::string matname = line.substr(7);
				auto result = mtls.emplace(matname, maxGroup + 1);
				if (result.second) {
					curGroup = ++maxGroup;
				} else {
					curGroup = result.first->second;
				}
			} else if (line.rfind("vn ", 0) == 0) {
				Vector v;
				sscanf(s, "vn %lf %lf %lf", &v[0], &v[1], &v[2]);
				normals.push_back(v);
			} else if (line.rfind("vt ", 0) == 0) {
				Vector v;
				sscanf(s, "vt %lf %lf", &v[0], &v[1]);
				uvs.push_back(v);
			} else if (line.rfind("v ", 0) == 0) {
				Vector pos, col;
				if (sscanf(s, "v %lf %lf %lf %lf %lf %lf", &pos[0], &pos[1], &pos[2], &col[0], &col[1], &col[2]) == 6) {
					for (int i = 0; i < 3; i++) col[i] = std::min(1.0, std::max(0.0, col[i]));
					vertexcolors.push_back(col);
				} else {
					sscanf(s, "v %lf %lf %lf", &pos[0], &pos[1], &pos[2]);
				}
				vertices.push_back(pos);
			}
			else if (line[0] == 'f') {
				int i[4], j[4], k[4], offset, nn;
				const char* cur = s + 1;
				TriangleIndices t;
				t.group = curGroup;

				// Try each face format: v/vt/vn, v/vt, v//vn, v
				if ((nn = sscanf(cur, "%d/%d/%d %d/%d/%d %d/%d/%d%n", &i[0], &j[0], &k[0], &i[1], &j[1], &k[1], &i[2], &j[2], &k[2], &offset)) == 9) {
					setFaceVerts(t, i[0], i[1], i[2]); 
					setFaceUVs(t, j[0], j[1], j[2]); 
					setFaceNormals(t, k[0], k[1], k[2]);
				} else if ((nn = sscanf(cur, "%d/%d %d/%d %d/%d%n", &i[0], &j[0], &i[1], &j[1], &i[2], &j[2], &offset)) == 6) {
					setFaceVerts(t, i[0], i[1], i[2]); 
					setFaceUVs(t, j[0], j[1], j[2]);
				} else if ((nn = sscanf(cur, "%d//%d %d//%d %d//%d%n", &i[0], &k[0], &i[1], &k[1], &i[2], &k[2], &offset)) == 6) {
					setFaceVerts(t, i[0], i[1], i[2]); 
					setFaceNormals(t, k[0], k[1], k[2]);
				} else if ((nn = sscanf(cur, "%d %d %d%n", &i[0], &i[1], &i[2], &offset)) == 3) {
					setFaceVerts(t, i[0], i[1], i[2]);
				}
				else continue;

				indices.push_back(t);
				cur += offset;

				// Fan triangulation for polygon faces (4+ vertices)
				while (*cur && *cur != '\n') {
					TriangleIndices t2;
					t2.group = curGroup;
					if ((nn = sscanf(cur, " %d/%d/%d%n", &i[3], &j[3], &k[3], &offset)) == 3) {
						setFaceVerts(t2, i[0], i[2], i[3]); 
						setFaceUVs(t2, j[0], j[2], j[3]); 
						setFaceNormals(t2, k[0], k[2], k[3]);
					} else if ((nn = sscanf(cur, " %d/%d%n", &i[3], &j[3], &offset)) == 2) {
						setFaceVerts(t2, i[0], i[2], i[3]); 
						setFaceUVs(t2, j[0], j[2], j[3]);
					} else if ((nn = sscanf(cur, " %d//%d%n", &i[3], &k[3], &offset)) == 2) {
						setFaceVerts(t2, i[0], i[2], i[3]); 
						setFaceNormals(t2, k[0], k[2], k[3]);
					} else if ((nn = sscanf(cur, " %d%n", &i[3], &offset)) == 1) {
						setFaceVerts(t2, i[0], i[2], i[3]);
					} else { 
						cur++; 
						continue; 
					}

					indices.push_back(t2);
					cur += offset;
					i[2] = i[3]; j[2] = j[3]; k[2] = k[3];
				}
			}
		}
		buildbox();
	}
	
	/*Vector compute_diag();
	int get_longest();
	Vector compute_barycenter();

	BBox compute_bbox(int starting_triangle, int end_triangle){

	}*/

	void buildbox(){
		bbox_min = vertices[0];
		bbox_max = vertices[0];

		for(int i=0; i < vertices.size(); i++){
			for(int j=0; j<3; j++){

				bbox_min[j] = std::min(bbox_min[j], vertices[i][j]);
				bbox_max[j] = std::max(bbox_max[j], vertices[i][j]);

			}
		}
	}

	bool intersect_box(const Ray& ray) const{
		double t_min = 0.0;
		double t_max = 1e30;

		for(int i=0; i<3; i++){
			double origin = ray.O[i];
			double direction = ray.u[i];

			double t_temp1 = (bbox_min[i] - origin) / direction;
			double t_temp2 = (bbox_max[i] - origin) / direction;

			if(t_temp1 > t_temp2){
				std::swap(t_temp1, t_temp2);
			}

			t_min = std::max(t_min, t_temp1);
			t_max = std::min(t_max, t_temp2);
			
			if(t_min > t_max){
				return false;
			}
		}
		return true;
	}
	
	// TODO ray-mesh intersection (labs 3 and 4)
	bool intersect(const Ray& ray, Vector& P, double& t, Vector& N) const {

		if(!intersect_box(ray)){
			return false;
		}

		bool hit = false;
		double t_min = 1e30;
		// lab 3 : for each triangle, compute the ray-triangle intersection with Moller-Trumbore algorithm
		for(int i=0; i < indices.size(); i++){
			TriangleIndices tri = indices[i];
			const Vector& A = vertices[tri.vtx[0]];
			const Vector& B = vertices[tri.vtx[1]];
			const Vector& C = vertices[tri.vtx[2]];

			Vector e1 = B - A;
			Vector e2 = C - A;
			Vector AO = A - ray.O;
			Vector N_geom = cross(e1, e2);
	
			double denom = dot(ray.u, N_geom);
			if(std::abs(denom) < 1e-12){
				continue;
			}
			double beta = dot(e2, cross(AO, ray.u)) / denom;
			double gamma = - dot(e1, cross(AO, ray.u)) / denom;
			double alpha = 1 - beta - gamma;
			if(beta < 0 || alpha < 0 || gamma < 0 || alpha > 1 || gamma > 1 || beta > 1){
				continue;
			}
			double t_temp = dot(AO, N_geom) / denom;
			if(t_temp < 1e-4){
				continue;
			}
			if(t_temp < t_min){
				t_min= t_temp;
				t = t_temp;

				P = ray.O + t * ray.u;
				hit = true;

				N = alpha * normals[tri.n[0]] + beta * normals[tri.n[1]] + gamma * normals[tri.n[2]];
				//N = N_geom;
				N.normalize();

			}
		}
		return hit;
		// lab 3 : once done, speed it up by first checking against the mesh bounding box
		// lab 4 : recursively apply the bounding-box test from a BVH datastructure
	}


	std::vector<TriangleIndices> indices;
	std::vector<Vector> vertices;
	std::vector<Vector> normals;
	std::vector<Vector> uvs;
	std::vector<Vector> vertexcolors;

	Vector bbox_min, bbox_max;
	bool bbox = false;

	struct BBox{
		int starting_triangle;
		int end_triangle;

		Vector bbmin;
		Vector bbmax;
	};

};



class Scene {
public:
	Scene() {};
	void addObject(const Object* obj) {
		objects.push_back(obj);
	}

	// returns true iif there is an intersection between the ray and any object in the scene
    // if there is an intersection, also computes the point of the *nearest* intersection P, 
    // t>=0 the distance between the ray origin and P (i.e., the parameter along the ray)
    // and the unit normal N. 
	// Also returns the index of the object within the std::vector objects in object_id
	bool intersect(const Ray& ray, Vector& P, double& t, Vector& N, int &object_id) const  {

		// TODO (lab 1): iterate through the objects and check the intersections with all of them, 
		// and keep the closest intersection, i.e., the one if smallest positive value of t
		double min = -1, t_tmp;
		bool is_intersection = false;
		for(int i=0; i < objects.size(); i++){
			Vector P_tmp, N_tmp;
			bool intersection = objects[i]->intersect(ray, P_tmp, t_tmp, N_tmp);
			if (intersection){
				if (min == -1 || t_tmp < min){
					is_intersection = true;
					min = t_tmp;
					P = P_tmp;
					N = N_tmp;
					t = t_tmp;
					object_id = i;
				}
			}
		}
		return is_intersection;
	}


Vector random_cos(const Vector N){
	double r1= uniform (engine[0]);
	double r2= uniform (engine[0]);
	double x = cos (2 * M_PI * r1) * sqrt(1 - r2);
	double y= sin (2 * M_PI * r1) * sqrt(1 - r2);
	double z = sqrt(r2);
	Vector T1(0, 0, 0);
	T1 = Vector(-N[1], N[0], 0);
	/*if(N[1] < N[2]){
		if(N[0] < N[1]){
			T1 = Vector(0, -N[2], N[1]);
		}
		else{
			T1 = Vector(-N[2],0 , N[0]);
		}
	}
	else if(N[0] < N[2]){
		T1 = Vector(0, -N[2], N[1]);
	}
	else{
		T1 = Vector(-N[1], N[0], 0);
	}*/

	T1.normalize();
	Vector T2 = cross(N, T1);
	T2.normalize();
	return x * T1 + y * T2 + z * N;
}

	// return the radiance (color) along ray
	Vector getColor(const Ray& ray, int recursion_depth) {

		if (recursion_depth >= max_light_bounce) return Vector(0, 0, 0);
		// TODO (lab 1) : if intersect with ray, use the returned information to compute the color ; otherwise black 
		// in lab 1, the color only includes direct lighting with shadows
		Vector P, N;
		double t;
		int object_id;
		Vector colour;
		if (intersect(ray, P, t, N, object_id)) {
			// Vector light = light_position - P;
			// dist = light.norm2();
			// light.normalize();
			// double attenuation = light_intensity/(4 * M_PI * dist);
			// Vector material = (objects[object_id]->albedo) / M_PI;
			// double solid_angle = dot(N, light);
			// colour = attenuation * material * solid_angle;
 
			if (objects[object_id]->mirror) {
				// return getColor in the reflected direction, with recursion_depth+1 (recursively)
				Vector reflected = ray.u - 2*dot(ray.u, N) * N;
				reflected.normalize();
				Ray reflected_ray(P + 1e-4 * N, reflected);
				return getColor(reflected_ray, recursion_depth+1);
			} // else

			if (objects[object_id]->transparent) { // optional

				// return getColor in the refraction direction, with recursion_depth+1 (recursively)
			} // else

			// test if there is a shadow by sending a new ray
		
			Vector light1 = light_position - P;
			double light1_norm = light1.norm();
			light1.normalize();

			Ray new_ray(P + 1e-4*N, light1);

			Vector P_tmp, N_tmp;
			double t_tmp;
			int object_id_shadow;
			bool shadow = false;

			bool is_intersection = intersect(new_ray, P_tmp, t_tmp, N_tmp, object_id_shadow);

			if (is_intersection){
				Vector P_prime = P_tmp;
				if(((P_prime - P).norm2()) <= (light1_norm * light1_norm)){
					shadow = true;
				}
			}// if there is no shadow, compute the formula with dot products etc.
			double vizibility = 1;
			// Sacha helped me with the shadow
			if(shadow){
				//vizibility = 0;
				Vector rand_vector = random_cos(N);
				Ray random_ray(P + 1e-4 * N, rand_vector);
				Vector albedo = objects[object_id]->albedo;
				Vector color1 = getColor (random_ray, recursion_depth + 1);
				Vector product(albedo[0] * color1[0], albedo[1] * color1[1], albedo[2] * color1[2]);
				return product;
				//return Vector(0, 0, 0);
			}

			else{
				Vector light2 = light_position - P;
				double light2_norm = light2.norm();
				double solid_angle = dot(N, light2 / light2_norm);
				if(solid_angle < 0){
					solid_angle = 0;
				}
				Vector Lo = Vector(0, 0, 0);
				//light2.normalize();
				double attenuation = light_intensity/(4 * M_PI * light2_norm * light2_norm);
				Vector material = (objects[object_id]->albedo) / M_PI;
				Lo = attenuation * material * solid_angle;
				//return colour;*/


			// TODO (lab 2) : add indirect lighting component with a recursive call
			

				Vector rand_vector = random_cos(N);
				Ray random_ray(P + 1e-4 * N, rand_vector);
				Vector albedo = objects[object_id]->albedo;
				Vector color1 = getColor (random_ray, recursion_depth + 1);
				Vector product(albedo[0] * color1[0], albedo[1] * color1[1], albedo[2] * color1[2]); 
				Lo = Lo + product;

				return Lo;
				}
		}

		return Vector(0, 0, 0);
	}

	std::vector<const Object*> objects;

	Vector camera_center, light_position;
	double fov, gamma, light_intensity;
	int max_light_bounce;
};


int main() {
	int W = 512;
	int H = 512;

	for (int i = 0; i<32; i++) {
		engine[i].seed(i);
	}

	Sphere center_sphere(Vector(0, 0, 0), 10., Vector(0.8, 0.8, 0.8), true);
	Sphere wall_left(Vector(-1000, 0, 0), 940, Vector(0.5, 0.8, 0.1));
	Sphere wall_right(Vector(1000, 0, 0), 940, Vector(0.9, 0.2, 0.3));
	Sphere wall_front(Vector(0, 0, -1000), 940, Vector(0.1, 0.6, 0.7));
	Sphere wall_behind(Vector(0, 0, 1000), 940, Vector(0.8, 0.2, 0.9));
	Sphere ceiling(Vector(0, 1000, 0), 940, Vector(0.3, 0.5, 0.3));
	Sphere floor(Vector(0, -1000, 0), 990, Vector(0.6, 0.5, 0.7));

	Scene scene;
	scene.camera_center = Vector(0, 0, 55);
	scene.light_position = Vector(-10,20,40);
	scene.light_intensity = 3E7;
	scene.fov = 60 * M_PI / 180.;
	scene.gamma = 2.2;    // TODO (lab 1) : play with gamma ; typically, gamma = 2.2
	scene.max_light_bounce = 2;

	//scene.addObject(&center_sphere);

	TriangleMesh mesh(Vector(0.8, 0.8, 0.8), false);
	mesh.readOBJ("cat.obj");
	mesh.scale_translate(0.6, Vector(0, -10, 0));

	scene.addObject(&mesh);
	scene.addObject(&wall_left);
	scene.addObject(&wall_right);
	scene.addObject(&wall_front);
	scene.addObject(&wall_behind);
	scene.addObject(&ceiling);
	scene.addObject(&floor);
	

	std::vector<unsigned char> image(W * H * 3, 0);

#pragma omp parallel for schedule(dynamic, 1)
	//double sum = 0;
	for (int i = 0; i < H; i++) {
		for (int j = 0; j < W; j++) {
			/*Vector color;

			// TODO (lab 1) : correct ray_direction so that it goes through each pixel (j, i)
			double X = j - W/2 + 0.5;
			double Y = H/2 - i - 0.5;
			double Z = -(W/(2*tan(scene.fov/2)));			
			Vector ray_direction(X,Y,Z);
			ray_direction.normalize();

			Ray ray(scene.camera_center, ray_direction);*/

			// TODO (lab 2) : add Monte Carlo / averaging of random ray contributions here
			Vector color = Vector(0, 0, 0);
			for(int k=0; k < 32; k ++){
				double r1= uniform (engine[0]);
				double r2= uniform (engine[0]);
				double x = sqrt (-2 * log ( r1 ) ) * cos (2 * M_PI * r2) * 0.3;
				double y= sqrt (-2 * log ( r1 ) ) * sin (2 * M_PI * r2 ) * 0.3;

				double X = j - W/2 + 0.5 + x;
				double Y = H/2 - i - 0.5 + y;
				double Z = -(W/(2*tan(scene.fov/2)));			
				Vector new_ray_direction(X,Y,Z);
				new_ray_direction.normalize();

				Ray new_ray(scene.camera_center, new_ray_direction);
				color = color + scene.getColor(new_ray, 0);

			}
			
			color = color / 32;

			// TODO (lab 2) : add antialiasing by altering the ray_direction here
			// TODO (lab 2) : add depth of field effect by altering the ray origin (and direction) here



			image[(i * W + j) * 3 + 0] = std::min(255., std::max(0., 255. * std::pow(color[0] / 255., 1. / scene.gamma)));
			image[(i * W + j) * 3 + 1] = std::min(255., std::max(0., 255. * std::pow(color[1] / 255., 1. / scene.gamma)));
			image[(i * W + j) * 3 + 2] = std::min(255., std::max(0., 255. * std::pow(color[2] / 255., 1. / scene.gamma)));
		}
	}


	stbi_write_png("image.png", W, H, 3, &image[0], 0);

	return 0;
}

/*clang++ -Xpreprocessor -fopenmp \
  -I/opt/homebrew/opt/libomp/include \
  main.cpp \
  -L/opt/homebrew/opt/libomp/lib -lomp*/