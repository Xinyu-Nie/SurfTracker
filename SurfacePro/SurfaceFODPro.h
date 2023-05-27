#include <vector>
#include <string>
#include <algorithm>
#include <math.h>
#include <iostream>
#include <time.h>
#include <thread>
#include <complex>
#include "../header/load_mesh.h"
#include "../header/InterpFODontoTRI.h"
#include "../header/math_aux.h"
#include "../header/matrixB.h"
#include "../header/BuildMeshNbrhood.h"
#include "../header/Masktrimesh.h"
#include "Complex2real.h"
#include "Decomrotation.h"
#include "Real2complex.h"
#include "SpYrotation.h"
#include "SpZrotation.h"
#include "Decomrotation.h"
#include "DecomrotationPro.h"
// #include "/ifs/loni/faculty/shi/spectrum/Student_2020/jruan/include/eigen3/Eigen/Dense"
using namespace std;
// using namespace Eigen;
void SurfaceFODPro(string& FOD, string& Mesh, string& ROI, string& output, float shrinkdis);
void SurfaceThread(int t1,int t2,  vector<vector<float>>& FODProjected, vector<vector<float>>& coordS, 
vector<vector<float>>& trgS, vector<vector<float>>& x, vector<vector<complex<double>>>& C2R, 
vector<vector<complex<double>>>& C2R1, vector<vector<complex<double>>>& R2C, 
vector<vector<complex<double>>>& R2C1, vector<vector<float>>& UY90p, 
vector<vector<float>>& UY90n,vector<vector<float>>& BS, vector<float>& weight);

inline void mat_multi(vector<vector<complex<double>>>&res, vector<vector<complex<double>>>& mat1, vector<vector<complex<double>>>& mat2){
    int mat1_m = mat1.size(), mat1_n = mat1[0].size(), mat2_m = mat2.size(), mat2_n = mat2[0].size();
    if (mat1_n != mat2_m) {cout<<"matrix unmatch error\n"; exit(0);}
    
    for (int i = 0;i<mat1_m;i++){
        for (int j = 0;j<mat2_n;j++){
            for (int k = 0;k<mat1_n;k++){
                res[i][j] += mat1[i][k] * mat2[k][j];
            }
        }
    }
    // return res;
}

inline void mat_multi(vector<vector<complex<double>>>& res, vector<vector<float>>& mat1, vector<vector<complex<double>>>& mat2){
    int mat1_m = mat1.size(), mat1_n = mat1[0].size(), mat2_m = mat2.size(), mat2_n = mat2[0].size();
    if (mat1_n != mat2_m) {cout<<"matrix unmatch error\n"; exit(0);}
    // vector<vector<complex<double>>> res(mat1_m,vector<complex<double>>(mat2_n, 0));
    for (int i = 0;i<mat1_m;i++){
        for (int j = 0;j<mat2_n;j++){
            for (int k = 0;k<mat1_n;k++){
                res[i][j] += (double)mat1[i][k] * mat2[k][j];
            }
        }
    }
    // return res;
}

inline void mat_multi(vector<vector<complex<double>>>& res, vector<vector<complex<double>>>& mat1, vector<vector<float>>& mat2){
    int mat1_m = mat1.size(), mat1_n = mat1[0].size(), mat2_m = mat2.size(), mat2_n = mat2[0].size();
    if (mat1_n != mat2_m) {cout<<"matrix unmatch error\n"; exit(0);}
    // vector<vector<complex<double>>> res(mat1_m,vector<complex<double>>(mat2_n, 0));
    for (int i = 0;i<mat1_m;i++){
        for (int j = 0;j<mat2_n;j++){
            for (int k = 0;k<mat1_n;k++){
                res[i][j] += mat1[i][k] * (double)mat2[k][j];
            }
        }
    }
    // return res;
}
inline void save(string path, vector<vector<float>>&data){
    ofstream outFile;
	outFile.open(path, ios::out);
	outFile<<data.size()<<' '<<data[0].size()<<endl;
	for (int i = 0;i<data.size();i++){
		for (int j = 0;j<data[0].size();j++){
			if (j<data[0].size()-1) outFile<< data[i][j]<<' ';
			else outFile<< data[i][j]<<endl;

		}
	}
    outFile.close();
}
inline void save(string path, vector<float>&data){
    ofstream outFile;
	outFile.open(path, ios::out);
	outFile<<data.size()<<endl;
	for (int i = 0;i<data.size();i++){
		outFile<< data[i]<<endl;

    }
    outFile.close();
	
}

inline void binary_save(FILE* pfile, vector<vector<float>>& data){
    // FILE* pfile = fopen(path.c_str(), "wb");
    int m = data.size(), n = data[0].size();
    fwrite(&m,sizeof(int),1,pfile);
    fwrite(&n,sizeof(int),1,pfile);
    for (int i=0;i<m;i++){
        fwrite(&data[i][0], sizeof(float), n, pfile);
    }
    // fclose(pfile);
}
inline void binary_save(FILE* pfile, vector<float>& data){
    int n = data.size();
    fwrite(&n, sizeof(int), 1, pfile);
    fwrite(&data[0], sizeof(float), n, pfile);
}