#include "DecomrotationPro.h"
vector<float> DecomrotationPro(float R[][3]){
    
    float* xr = R[0];
    float* zr = R[2];
    normalize(xr);
    normalize(zr);
    float x[3] = {1,0,0}, y[3] = {0,1,0}, z[3] = {0,0,1};
    float n[3];
    cross(n,z,zr);
    float qn = norm(n);
    if (qn>0) {n[0]/=qn; n[1]/=qn; n[2]/=qn;}
    float alfa = acos(dot(n,y));
    if (n[0]>0) alfa = 2*PI - alfa;
    float Rz[3][3] = {{cos(alfa), -sin(alfa), 0}, {sin(alfa), cos(alfa), 0}, {0,0,1}};
    float x1[3];
    for (int i=0;i<3;i++) x1[i] = dot(Rz[i],x);
    float beta = acos(dot(z,zr));
    float Rn[3][3];
    rotate(Rn, n, beta);
    float zr1[3];
    for (int i =0;i<3;i++) zr1[i] = dot(Rn[i],z);
    float pz1 = abs(dot(zr1,zr));
    float Rnr[3][3];
    rotate(Rnr, n,2*PI-beta);
    float zr11[3];
    for (int i = 0;i<3;i++) zr11[i] = dot(Rnr[i],z);
    float pz2 = abs(dot(zr11,zr));
    if (pz1 < pz2){
        beta = 2*PI-beta;
    }
    float tmp[3][3];
    rotate(tmp,n,beta);
    float x2[3];
    for (int i =0;i<3;i++) x2[i] = dot(tmp[i],x1);
    float gamma = acos(dot(x2,xr));
    float Rzr[3][3];
    rotate(Rzr,zr,gamma);
    float xr1[3];
    for (int i=0;i<3;i++) xr1[i] = dot(Rzr[i],x2);
    float px1 = abs(dot(xr1,xr));
    float Rzrr[3][3];
    rotate(Rzrr, zr, 2*PI-gamma);
    float xr11[3];
    for (int i =0;i<3;i++) xr11[i] = dot(Rzrr[i], x2);
    float px2 = abs(dot(xr11,xr));
    if (px1<px2) gamma = 2*PI-gamma;
    return {alfa, beta, gamma};
}