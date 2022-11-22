#include"Real2complex.h"

vector<vector<complex<double>>> Real2complex(int maxOrder){
    int N = (maxOrder+1)*(maxOrder+2)/2;
    vector<vector<complex<double>>> U(N,vector<complex<double>>(N,0));
    U[0][0] = 1;
    int t = 1;
    for (int order = 2; order<=maxOrder;order+=2){
        int len = order*2+1;
        vector<vector<complex<double>>> A(len, vector<complex<double>>(len, 0));
        A[order][order]=1;
        float k1=sqrt(0.5);
        complex<double> k2(0,k1);
        for (int m = 1;m<=order;m++){
            int flag;
            if (m%2) flag = -1;
            else flag = 1;
            A[m+order][order-m] = k2*complex<double>(flag);
            A[m+order][order+m] = complex<double>(k1)*complex<double>(flag);
            A[order-m][order-m] = -complex<double>(k2);
            A[order-m][order+m]= complex<double>(k1);
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