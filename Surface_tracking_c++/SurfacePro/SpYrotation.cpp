#include "SpYrotation.h"

vector<vector<float>> SpYrotation(int maxOrder, float beta){ // U[74][66]-[69] check
    int N = (maxOrder+1)*(maxOrder+2)/2;
    vector<vector<float>> U(N,vector<float>(N,0));
    U[0][0] = 1;
    int t = 1;
    for (int order = 2;order<=maxOrder;order+=2){
        int len = 2*order+1;
        auto D = WignerDmatrix(order, beta);
        for (int i = 0;i<len;i++){
            for (int j = 0;j<len;j++){
                U[t+i][t+j] = D[i][j];
            }
        }
        t += len;
    }
    return U;
}