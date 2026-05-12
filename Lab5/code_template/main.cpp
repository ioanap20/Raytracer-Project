#include <vector>
#include <cmath>
#include <random>
#include <map>
#include <string>
#include <fstream>
#include <iostream>

#define _CRT_SECURE_NO_WARNINGS 1
#include <vector>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

//this is from Lab1
static std::default_random_engine engine[32];
static std::uniform_real_distribution<double> uniform(0, 1);

// got the code for Vector class from Lab1
class Vector {
public:
	explicit Vector(double x = 0, double y = 0, double z = 0) {
		data[0] = x;
		data[1] = y;
		data[2] = z;
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

Vector random_direction(){
	double r1= uniform (engine[0]);
	double r2= uniform (engine[0]);

	double x = cos (2 * M_PI * r1) * sqrt(1 - r2);
	double y= sin (2 * M_PI * r1) * sqrt(1 - r2);
	double z = 1 - 2 * r2;

	return Vector(x, y, z);
}

void sliced_optimal_transport(std:: vector <double>& I, std:: vector<double>& M, int n, int nbiter){
	std::vector<std::pair<double, int>> projI(n);
	std::vector<std::pair<double, int>> projM(n);
	for(int iter=0; iter < nbiter; iter++){
		Vector v = random_direction();

		for(int i=0; i<n; i++){

			double ri = I[i*3];
			double gi = I[i*3 + 1];
			double bi = I[i*3 + 2];

			double rm = M[i*3];
			double gm = M[i*3+1];
			double bm = M[i*3 +2];

			projI[i] = {dot(Vector(ri, gi, bi), v), i};
			projM[i] = {dot(Vector(rm, gm, bm), v), i};
		}

			std :: sort(projI.begin(), projI.end());
			std :: sort(projM.begin(), projM.end());

			for(int i=0; i<n; i++){
				int index1 = std::get<1>(projI[i]);
				int index2 = std::get<1>(projM[i]);

				double delta = (std::get<0>(projM[i]) - std::get<0>(projI[i]));

				I[index1 * 3] = I[index1 * 3] + delta * v[0];
				I[index1 * 3 + 1] = I[index1 * 3 + 1] + delta * v[1];
				I[index1 * 3 + 2] = I[index1 * 3 + 2] + delta * v[2];

			}

	}
}

int main() {

	int W, H, C;
	
	//stbi_set_flip_vertically_on_load(true);
	unsigned char *image = stbi_load("8733654151_b9422bb2ec_k.jpg",
                                 &W,
                                 &H,
                                 &C,
                                 STBI_rgb);

	int W_model, H_model, C_model;
	
	//stbi_set_flip_vertically_on_load(true);
	unsigned char *model = stbi_load("redim.jpg",
                                 &W_model,
                                 &H_model,
                                 &C_model,
                                 STBI_rgb);

	std::vector<double> image_double(W*H*3);
	std::vector<double> model_double(W_model*H_model*3);
	for (int i=0; i<W*H*3; i++)
		{
			image_double[i] = image[i];
			model_double[i] = model[i];
		}
	
	sliced_optimal_transport(image_double, model_double, W*H - 1, 100);
	
	std::vector<unsigned char> image_result(W*H * 3, 0);
	for (int i = 0; i < H; i++) {
		for (int j = 0; j < W; j++) {
			double r = image_double[(i*W+j)*3+0]*0.5;
			double g = image_double[(i*W+j)*3+1]*0.5;
			double b = image_double[(i*W+j)*3+2]*0.5;

			if(r<0){
				r=0;
			}
			if(r>255){
				r = 254;
			}

			if(g<0){
				g=0;
			}
			if(g>255){
				g = 254;
			}

			if(b<0){
				b=0;
			}
			if(b>255){
				b = 254;
			}

			image_result[(i*W + j) * 3 + 0] = r;
			image_result[(i*W + j) * 3 + 1] = g;
			image_result[(i*W + j) * 3 + 2] = b;
		}
	}
	stbi_write_png("image.png", W, H, 3, &image_result[0], 0);

	return 0;
}