// function [tract,flagendinge,flagendingp,biterminflag, TriCount] = trackingonsufaceparallel(coordS,trgS,tricenter,THETA,FOD,seedflag,trinerbor,Triendflag,Tribeginflag,minFOD,trackdegree)

#include<vector>
// #include "rejectsampling.h"
#include "math_aux.h"
#include "../lsqr/lsqrDense.h"
#include "projection.h"
#include "traveltrigparallel.h"
using namespace std;

void track(lsqrDense& lsqr, int a, int& cpedge, int& cnedge, vector<vector<float> >& tractline, vector<vector<float> >& tractbidirection, float* vd, vector<float>&v0, vector<float>& v1, vector<float>& v2, float seedflag);


void trackingonsufaceparallel(vector<vector<float> >& tract, int& flagendinge, int& flagendingp, vector<float>& biterminflag, vector<vector<float> >& coordS, vector<vector<float> >& trgS, vector<vector<float> >& tricenter, vector<float>& THETA, vector<vector<float> >& FOD, float seedflag, vector<vector<float> >& trinerbor,vector<float>& Triendflag,vector<float>& Tribeginflag, float minFOD, float trackdegree);
