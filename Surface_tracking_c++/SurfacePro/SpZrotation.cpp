#include "SpZrotation.h"

vector<vector<complex<double>>> SpZrotation(int maxOrder, float alfa){
    int N = (maxOrder+1)*(maxOrder+2)/2;
    vector<vector<complex<double>>> U(N, vector<complex<double>>(N,0));
    U[0][0] = 1;
    int t = 1;
    for (int order = 2;order<=maxOrder;order+=2){
        int len = order*2+1;
        vector<vector<complex<double>>> A(len, vector<complex<double>>(len,0));
        A[order][order] = 1;
        for (int m = 1;m<=order;m++){
            A[order+m][order+m] = exp(c_i*(double)m*(double)alfa);
            A[order-m][order-m] = exp(-c_i*(double)m*(double)alfa);
        }
        for (int i = 0;i<len;i++){
            for (int j = 0;j<len;j++){
                U[t+i][t+j] = A[i][j];
            }
        }
        t+=len;
    }
    return U;
}