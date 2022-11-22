#ifndef SRC_MATH_MATH_AUX_H_
#define SRC_MATH_MATH_AUX_H_

#ifdef BUILD_FOR_WINDOWS
#include <windows.h>
#include<time.h>
#endif

#include <iostream>
#include <cmath>
#include<vector>
#include "../lsqr/lsqrDense.h"
#include <random>

using namespace std;

#define SQRT2	   1.41421356237309504880168872420969807856967187537694807
#define PI 		   3.14159265358979323846264338327950288419716939937510582
#define TWOPI	   6.28318530717958647692528676655900576839433879875021164
#define PIOVERTWO  1.57079632679489661923132169163975144209858469968755291

namespace MATH_AUX {
	static std::mt19937 gen(time(0));
	static int legendre_flag = assoc_legendre(2,1,0.5)<0;
}


// class RandomDoer{
// 	public :
// 	RandomDoer();
// 	~RandomDoer();
// 	private:
// 	std::mt19937 gen;
// 	std::uniform_int_distribution<int>    *unidis_int;
// 	std::uniform_real_distribution<float> *unidis_01;
// }

template <typename T>
inline T sum(T* a){
	return a[0] + a[1] + a[2];
}

inline void verifyUnitRange(float *vec) {
	for (int i=0;i<3;i++){
		if (vec[i]> 1) vec[i]= 0.99999999999999999999999999999999999999;
		if (vec[i]<-1) vec[i]=-0.99999999999999999999999999999999999999;
	}
}

inline void normalize(float *vec) {
	float normi = 1.0/std::sqrt(vec[0]*vec[0]+vec[1]*vec[1]+vec[2]*vec[2]);
	vec[0] *= normi;
	vec[1] *= normi;
	vec[2] *= normi;
}

inline float norm(float *vec) {
	return std::sqrt(vec[0]*vec[0]+vec[1]*vec[1]+vec[2]*vec[2]);
}
template<typename T>
inline T dot(T* A, T* B) {
	return A[0]*B[0]+A[1]*B[1]+A[2]*B[2];
}

inline float dotsum(vector<float>& A, vector<float>& B){
	int n = A.size();
	if (B.size()!=n) {cout<<"incompatible size"<<endl;return 0 ;}
	float sum = 0;
	for (int i =0;i<n;i++) sum+= A[i]*B[i];
}
inline void dot(int rows, int cols, float** A, float**B, float** out){
	for (int i = 0;i<rows; i++){
		for (int j = 0;j<cols;j++){
			out[i][j] = A[i][j] * B[i][j];
		}
	}	
}

inline void cross(float *out, float* A, float *B) {
	out[0] = A[1]*B[2] - A[2]*B[1];
	out[1] = A[2]*B[0] - A[0]*B[2];
	out[2] = A[0]*B[1] - A[1]*B[0];
}

inline void rotate(float R[][3], float* direction, float theta){
	// vector<vector<float>> R(3, vector<float>(3,0));
	// vector<vector<float> > R(3, vector<float>(3));
	float ux = direction[0], uy = direction[1], uz = direction[2], cth = cos(theta), sth = sin(theta);
	/*
	R=[cth+(ux^2)*(1-cth) ux*uy*(1-cth)-uz*sth ux*uz*(1-cth)+uy*sth; 
	uy*ux*(1-cth)+uz*sth cth+(uy^2)*(1-cth) uy*uz*(1-cth)-ux*sth; 
	uz*ux*(1-cth)-uy*sth uz*uy*(1-cth)+ux*sth cth+(uz^2)*(1-cth)];
	*/
	R[0][0] = cth+ux*ux*(1-cth); R[0][1] = ux*uy*(1-cth)-uz*sth; R[0][2] = ux*uz*(1-cth)+uy*sth;
	R[1][0] = uy*ux*(1-cth)+uz*sth; R[1][1] = cth+(uy*uy)*(1-cth); R[1][2] = uy*uz*(1-cth)-ux*sth;
	R[2][0] = uz*ux*(1-cth)-uy*sth; R[2][1] = uz*uy*(1-cth)+ux*sth; R[2][2] = cth+(uz*uz)*(1-cth);
	// return R;
}

inline float** rotate(float direction[], float cth, float sth){
	// vector<vector<float>> R(3, vector<float>(3,0));
	float **R = new float*[3];
	R[0] = new float[3];
	R[1] = new float[3];
	R[2] = new float[3];
	float ux = direction[0], uy = direction[1], uz = direction[2];
	/*
	R=[cth+(ux^2)*(1-cth) ux*uy*(1-cth)-uz*sth ux*uz*(1-cth)+uy*sth; 
	uy*ux*(1-cth)+uz*sth cth+(uy^2)*(1-cth) uy*uz*(1-cth)-ux*sth; 
	uz*ux*(1-cth)-uy*sth uz*uy*(1-cth)+ux*sth cth+(uz^2)*(1-cth)];
	*/
	R[0][0] = cth+ux*ux*(1-cth); R[0][1] = ux*uy*(1-cth)-uz*sth; R[0][2] = ux*uz*(1-cth)+uy*sth;
	R[1][0] = uy*ux*(1-cth)+uz*sth; R[1][1] = cth+(uy*uy)*(1-cth); R[1][2] = uy*uz*(1-cth)-ux*sth;
	R[2][0] = uz*ux*(1-cth)-uy*sth; R[2][1] = uz*uy*(1-cth)+ux*sth; R[2][2] = cth+(uz*uz)*(1-cth);
	return R;
}

#endif
