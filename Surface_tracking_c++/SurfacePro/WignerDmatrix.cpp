#include "WignerDmatrix.h"
double factorial(int n){
    if (n==0) return 1;
    if (n==1) return n;
    return (double)n*factorial(n-1);
}
vector<vector<float>> WignerDmatrix(int n, float beta){
    int len = 2*n+1;
    vector<vector<float>> D(len,vector<float>(len,0));
    for (int m = 0;m<=n;m++){
        for (int k = m;k<=n;k++){
    
            double x[1]={cos(beta)};
            auto tmp = j_polynomial(1,n-k,k-m,k+m,x);
            float Pbeta = tmp[n-k];
            delete tmp;
            D[n+m][n+k] = (1/(pow(2.0,k))) * sqrt(1.0*factorial(n-k)*factorial(n+k)/(factorial(n-m)*factorial(n+m))) * (pow(sin(beta),k-m)) * (pow(1+cos(beta),m)) * Pbeta;
            if (k>m){
                D[n+k][n+m] = D[n+m][n+k]*(pow(-1, k-m));
                D[n-k][n-m] = D[n+m][n+k];
            }
            D[n-m][n-k] = D[n+m][n+k]*pow(-1,k-m);
            
            if (m>0) { //m>0
                double x[1] = {cos(PI-beta)};
                // x[0] = cos(PI-beta);
                auto tmp2 = j_polynomial(1,n-k,k-m,k+m,x);
                float Pbetar = tmp2[n-k];
                delete tmp2;
                D[n+m][n-k] = (pow(-1,n-m) * (1/pow(2,k))) * sqrt(1.0*factorial(n-k)*factorial(n+k)/(factorial(n-m)*factorial(n+m))) * (pow(sin(PI-beta),k-m)) * (pow(1+cos(PI-beta),m)) * Pbetar;
                if (k>m){
                    D[n+k][n-m]=D[n+m][n-k];
                    D[n-k][n+m]=(pow(-1,-k-m))*D[n+m][n-k];
                }
                D[n-m][n+k] = (pow(-1,-k-m)) * D[n+m][n-k];
            }
        }
    }
    return D;
}

// vector<vector<double>> WignerDmatrix(int n, double beta){
//     int len = 2*n+1;
//     vector<vector<double>> D(len,vector<double>(len,0));
//     for (int m = 0;m<=n;m++){
//         for (int k = m;k<=n;k++){
    
//             double x[1]={cos(beta)};
//             auto tmp = j_polynomial(1,n-k,k-m,k+m,x);
//             double Pbeta = tmp[n-k];
//             delete tmp;
//             D[n+m][n+k] = (1/(pow(2.0,k))) * sqrt(1.0*factorial(n-k)*factorial(n+k)/(factorial(n-m)*factorial(n+m))) * (pow(sin(beta),k-m)) * (pow(1+cos(beta),m)) * Pbeta;
//             if (k>m){
//                 D[n+k][n+m] = D[n+m][n+k]*(pow(-1, k-m));
//                 D[n-k][n-m] = D[n+m][n+k];
//             }
//             D[n-m][n-k] = D[n+m][n+k]*pow(-1,k-m);
            
//             if (m>0) { //m>0
//                 double x[1] = {cos(PI-beta)};
//                 // x[0] = cos(PI-beta);
//                 auto tmp2 = j_polynomial(1,n-k,k-m,k+m,x);
//                 double Pbetar = tmp2[n-k];
//                 delete tmp2;
//                 D[n+m][n-k] = (pow(-1,n-m) * (1/pow(2,k))) * sqrt(1.0*factorial(n-k)*factorial(n+k)/(factorial(n-m)*factorial(n+m))) * (pow(sin(PI-beta),k-m)) * (pow(1+cos(PI-beta),m)) * Pbetar;
//                 if (k>m){
//                     D[n+k][n-m]=D[n+m][n-k];
//                     D[n-k][n+m]=(pow(-1,-k-m))*D[n+m][n-k];
//                 }
//                 D[n-m][n+k] = (pow(-1,-k-m)) * D[n+m][n-k];
//             }
//         }
//     }
//     return D;
// }