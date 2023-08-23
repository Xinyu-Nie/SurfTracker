
#include<vector>
#include<cmath>
#include"math_aux.h"

// #include <boost/math/special_functions/legendre.hpp>
// #define PI 		   3.1416
using namespace std;
inline float prod(int l, int r){
    float res = 1;
    while (l<=r){
        res*=l;
        l++;
    }
    return res;
}
inline void matrixB(vector<vector<float>>& B, vector<vector<float>>& A, int maxOrder){ // A size: 3*n
    /*
    %% Input: maxOrder = 0,2,4,6,8
%         A = matrix of size n x 3  
%% P = P^m_l(cos(theta))
% P = (P_1(cos(theta_1)).....P_R(cos(theta_1))
%      P_1(cos(theta_2)).....P_R(cos(theta_2))
%
%      .....................
%      P_1(cos(theta_n)).....P_R(cos(theta_n)))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% T = (...sin(2*phi_1) sin(phi_1) 1 sin(phi_1) sin(2*phi_1)...
%  ...sin(2*phi_2) sin(phi_2) 1 sin(phi_2) sin(2*phi_2)...
%      ..........................................
%  ...sin(2*phi_n) sin(phi_n) 1 sin(phi_n) sin(2*phi_n)...)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/
    int rows = A.size(), cols = A[0].size();
    float* thetavec = new float[cols];
    for (int i = 0;i<cols;i++){
        thetavec[i] = atan2(sqrt(A[0][i]*A[0][i] + A[1][i]*A[1][i]), A[2][i]);
    }
   
    float* phivec = new float[cols];
    for (int i = 0;i<cols;i++){
        phivec[i] = atan2(A[1][i], A[0][i]);
    }

    int* order = new int[maxOrder/2+1];
    int lengthOfOrder = maxOrder/2+1;
    for (int i = 0;i<lengthOfOrder;i++){
        order[i] = i*2;
    }
    
    int* R = new int[lengthOfOrder];
    for (int i =0;i<lengthOfOrder;i++){
        R[i] = (order[i]+1) * (order[i]+2)/2;
    }

    // float **P = new float*[cols], **T = new float*[cols];// n* degree
    vector<vector<float>> P(cols), T(cols);
    float* scale = new float[R[lengthOfOrder-1]];
    for (int i=0;i<cols;i++) {
        P[i] = vector<float>(R[lengthOfOrder-1],0);
        T[i] = vector<float>(R[lengthOfOrder-1],0);
        // P[i] = new float[R[lengthOfOrder-1]];
        // T[i] = new float[R[lengthOfOrder-1]];
        // P[i][0] = boost::math::legendre_p(0,0,cos(thetavec[i]));
#if defined(BUILD_FOR_WINDOWS) || defined(BUILD_FOR_MAC)
        P[i][0] = legendre_p(0,0,cos(thetavec[i]));
#else
        P[i][0] = assoc_legendre(0, 0, cos(thetavec[i]));
#endif
        T[i][0] = 1;
    }
    scale[0] = sqrt(1/(4.0*PI));

    int indx = 0;
    for (int i =1;i<lengthOfOrder;i++){
        // temp = legendre(order(i),cos(thetavec)); temp = temp';
        // float** tmp = new float*[cols];
        vector<vector<float>> tmp(cols);
        // for (int j=0;j<cols;j++) tmp[j] = new float[order[i]+1];
        for (int j=0;j<cols;j++) tmp[j] = vector<float>(order[i]+1,0);
        for (int j =0;j<cols;j++){
            for (int m = 0;m<=order[i];m++){
                // tmp[j][m] = boost::math::legendre_p(order[i],m,cos(thetavec[j]));
#if defined(BUILD_FOR_WINDOWS) || defined(BUILD_FOR_MAC)
                tmp[j][m] = legendre_p(order[i], m, cos(thetavec[j]));
#else
                tmp[j][m] = assoc_legendre(order[i], m, cos(thetavec[j]));
#endif
                if (!MATH_AUX::legendre_flag&& m%2) tmp[j][m] = -tmp[j][m];//check assoc_legendre
            }
        }
        //P(:,(R(i-1)+1):R(i)) = [temp(:,order(i)+1:-1:2) temp];
        // indx R[i-1] ==> R[i]-1
        for (int j = 0;j<cols;j++){
            indx = R[i-1];
            for (int m = order[i];m>=1;m--) {
                P[j][indx] = tmp[j][m];
                T[j][indx] = sin(phivec[j]*m);    
                indx++;
            }
            P[j][indx] = tmp[j][0];
            T[j][indx] = 1;
            indx ++;
            for (int m = 1;m<=order[i];m++){
                P[j][indx] = tmp[j][m];
                T[j][indx] = cos(phivec[j]*m);    
                indx++;
            }
        }
        float* tmpv = new float[order[i]+1];
        float numerator = (2*order[i]+1)/(2*PI);
        for (int m = 1;m<=order[i];m++) tmpv[m] = sqrt(numerator/prod(order[i]-m+1, order[i]+m));
        indx = R[i-1];
        for (int m = order[i];m>=1;m--) scale[indx++] = tmpv[m];
        scale[indx++] = sqrt((2*order[i]+1)/(4*PI));
        for (int m = 1;m<=order[i];m++) scale[indx++] = tmpv[m];
        delete tmpv;
    }
    B = vector<vector<float>>(cols, vector<float>(R[lengthOfOrder-1],0));
    // cout<<98<<endl;
    // for (int i = 0;i<cols;i++) B[i] = new float[R[lengthOfOrder-1]];
    //B = repmat(scale,rows,1).*P.*T;
    for (int i = 0; i<cols;i++){
        for (int j = 0;j<R[lengthOfOrder-1];j++){
            B[i][j] = scale[j] * P[i][j] * T[i][j];
        }
    }
    delete thetavec;
    delete phivec;
    delete order;
    delete R;
    // for (int i = 0;i<cols;i++){
    //     delete P[i];
    //     delete T[i];
    // }
    // delete[] P;
    // delete[] T; 
    delete scale;

    // return B;
    
}