#include<iostream>
#include<fstream>
#include<vector>
#include<stdlib.h>
#include<algorithm>
#include"math_aux.h"
using namespace std;
inline vector<float> randperm(int l, int r, int n){
    vector<float> randseed;
    for (int i=0;i<n;i++) {
        randseed.push_back(uniform_int_distribution<>(l,r)(MATH_AUX::gen));
    }
    
    return randseed;
}

inline void read_matrix(vector<vector<float> >& mat, string path){
    int m,n;
    FILE* pfile = fopen(path.c_str(), "r");
    fscanf(pfile, "%d %d", &m,&n);
    mat = vector<vector<float> >(m, vector<float>(n,0));
    for (int i =0;i<m;i++){
        for (int j = 0;j<n;j++){
            fscanf(pfile, "%f", &mat[i][j]);
        }
    }
    
}
inline void read_vector(vector<float>& mat, string path){
    int n;
    FILE* pfile = fopen(path.c_str(), "r");
    fscanf(pfile, "%d", &n);
    mat = vector<float>(n,0);
    for (int i = 0;i<n;i++){
        fscanf(pfile, "%f", &mat[i]);
    }

}
inline void binary_read(FILE* pfile, vector<vector<float>>& data){
    int m, n;
    fread(&m, sizeof(int), 1, pfile);
    fread(&n, sizeof(int), 1, pfile);
    data = vector<vector<float>>(m,vector<float>(n,0));
    for (int i=0;i<m;i++) fread(&data[i][0], sizeof(float), n, pfile);
}
inline void binary_read(FILE* pfile, vector<float>& data){
    int n;
    fread(&n, sizeof(int), 1, pfile);
    data = vector<float>(n,0);
    fread(&data[0], sizeof(float), n, pfile);
}

void Write2File(string file, vector<float> &LineInd, vector<vector<float>> &PtVec, vector<int>& dim, vector<float>& voxel_size);

void Ufibertracking(string output, string Mesh, string ROI, string Sulcuslabel,string Gyruslabel1, string Gyruslabel2, int iteration, vector<vector<float> >& coordS, vector<vector<float> >& trgS, vector<vector<float> >&TRIcen, vector<float>& THETA, vector<vector<float> >& FODProjected, vector<vector<float> >& trinerbor, vector<int>& dim, vector<float>& res, float minFOD, float trackdegree, int Gyrusflag);
