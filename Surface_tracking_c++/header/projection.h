#include <vector>
#include "math_aux.h"
using namespace std;

inline vector<float> projection(vector<vector<float> >& coordS, vector<vector<float> >& trgS, int index, float* Nv, float& q){
    auto v1 = coordS[trgS[index][0]];
    auto v2 = coordS[trgS[index][1]];
    auto v3 = coordS[trgS[index][2]];
    float v12[3], v13[3], v23[3];
    for (int i=0;i<3;i++){
        v12[i] = v2[i] - v1[i];
        v13[i] = v3[i] - v1[i];
        v23[i] = v3[i] - v2[i];
    }
    // Nv = new float[3];
    cross(Nv, v12,v13);
    q = norm(Nv);
    
}