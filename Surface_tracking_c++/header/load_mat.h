#include<mat.h>
#include<iostream>
#include <fstream>
#include <sstream>
using namespace std;
inline void load_mat(string mat_path, vector<vector<float>>& FODProjected, vector<vector<float>>& trgS, vector<vector<float>>& coordS, vector<vector<float>>& TRIcen, vector<vector<float>>& trinerbor, vector<float>& THETA, vector<int>& dim, vector<float>& res){
    MATFile *pmatFile = NULL;
    mxArray *pMxArray = NULL;
    double *initA;
    pmatFile = matOpen(mat_path.c_str(), "r");
    ofstream os("/ifs/loni/faculty/shi/spectrum/Student_2020/jruan/SurfaceTrack_v1/saved_data/r.txt");
    // FODProjected
    pMxArray = matGetVariable(pmatFile, "FODProjected");
    initA = (double*) mxGetData(pMxArray);
    size_t M = mxGetM(pMxArray);
    size_t N = mxGetN(pMxArray);
    FODProjected = vector<vector<float>>(M,vector<float>(N,0));
    for (int i=0; i<M; i++)
    {
        for (int j=0; j<N; j++)
            {
                FODProjected[i][j] = (float)initA[M*j+i];
            }
    }
    mxFree(initA);
    //coordS
    pMxArray = matGetVariable(pmatFile, "coordS");
    initA = (double*) mxGetData(pMxArray);
    M = mxGetM(pMxArray);
    N = mxGetN(pMxArray);
    coordS = vector<vector<float>>(M,vector<float>(N,0));
    for (int i=0; i<M; i++)
    {
        for (int j=0; j<N; j++)
            {
                coordS[i][j] = (float)initA[M*j+i];
            }
    }
    mxFree(initA);
    // trgS
    
    pMxArray = matGetVariable(pmatFile, "trgS");
    initA = (double*) mxGetData(pMxArray);
    M = mxGetM(pMxArray);
    N = mxGetN(pMxArray);
    trgS = vector<vector<float>>(M,vector<float>(N,0));
    for (int i=0; i<M; i++)
    {
        for (int j=0; j<N; j++)
            {
                trgS[i][j] = (float)initA[M*j+i]-1;
            }
    }
    mxFree(initA);
    // dim
    pMxArray = matGetVariable(pmatFile, "dim");
    initA = (double*) mxGetData(pMxArray);
    M = mxGetM(pMxArray);
    N = mxGetN(pMxArray);
    for (int i=0; i<M; i++)
    {
        for (int j=0; j<N; j++)
            {
                dim[j] = (int)initA[M*j+i];
            }
    }
    mxFree(initA);
    // pixdim
    pMxArray = matGetVariable(pmatFile, "res");
    initA = (double*) mxGetData(pMxArray);
    M = mxGetM(pMxArray);
    N = mxGetN(pMxArray);
    for (int i=0; i<M; i++)
    {
        for (int j=0; j<N; j++)
            {
                res[j] = (float)initA[M*j+i];
            }
    }
    mxFree(initA);
    // THETA
    pMxArray = matGetVariable(pmatFile, "THETA");
    initA = (double*) mxGetData(pMxArray);
    M = mxGetM(pMxArray);
    N = mxGetN(pMxArray);
    THETA = vector<float>(N);
    for (int i=0; i<M; i++)
    {
        for (int j=0; j<N; j++)
            {
                THETA[j] = (float)initA[M*j+i];
            }
    }
    mxFree(initA);
    //TRIcen
    pMxArray = matGetVariable(pmatFile, "TRIcen");
    initA = (double*) mxGetData(pMxArray);
    M = mxGetM(pMxArray);
    N = mxGetN(pMxArray);
    TRIcen = vector<vector<float>>(M, vector<float>(N,0));
    for (int i=0; i<M; i++)
    {
        for (int j=0; j<N; j++)
            {
                TRIcen[i][j] = (float)initA[M*j+i];
            }
    }
    mxFree(initA);
    //trinerbor
    pMxArray = matGetVariable(pmatFile, "trinerbor");
    initA = (double*) mxGetData(pMxArray);
    M = mxGetM(pMxArray);
    N = mxGetN(pMxArray);
    trinerbor = vector<vector<float>>(M, vector<float>(N,-1));
    for (int i=0; i<M; i++)
    {
        for (int j=0; j<N; j++)
            {
                trinerbor[i][j] = (float)initA[M*j+i]-1;
                os<<trinerbor[i][j]<<' ';
            }
            os<<endl;
    }
    mxFree(initA);
    matClose(pmatFile);
}