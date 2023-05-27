#include <vector>
#include "../header/math_aux.h"
#include "../header/rejectsampling.h"
using namespace std;
void track(int a, lsqrDense& lsqr, float* vd, vector<float>&v0, vector<float>& vin1, vector<float>& vin2, vector<float>& outpoint, int &outedgef, int intri);
void traveltrigparallel(vector<float>& outpoint, int& outedgef, int& terminflag, vector<vector<float> >& coordS, vector<vector<float> >& trgS, vector<float>& inpoint, vector<float>& prepoint, vector<float>& THETA, vector<vector<float> >& FOD, int intri, int pretri, int inedge, float minFOD, float trackdegree);
