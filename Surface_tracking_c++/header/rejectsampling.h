#include <vector>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <stdlib.h>
#include <time.h>
#include "math_aux.h"
using namespace std;
inline int rejectsampling(vector<float>& FOD, float threshold, float minfod){
    mt19937& gen = MATH_AUX::gen;
    int Lf = FOD.size();
    int direction;
    float sumfod=0;
    float fodmax = 0;
    for (int i =0;i<FOD.size();i++) {
        if (FOD[i]>fodmax) fodmax = FOD[i];
        sumfod+=FOD[i];
    }
    // float fodmax = *(max_element(FOD.begin(), FOD.end()));
    if (fodmax>0){
        // float sumfod = accumulate(FOD.begin(), FOD.end(), 0);
        vector<float> FODnorm(FOD.size(), 0);
        float fmax = -1e5;
        for (int i = 0;i<FOD.size();i++){
            FODnorm[i] = FOD[i]/sumfod;
            if (FODnorm[i] > fmax) fmax = FODnorm[i];
        }
        // float fmax = fodmax/sumfod;
        float fmin = minfod;
        float M = threshold * fmax;
        int Y = uniform_int_distribution<>(0,Lf-1)(gen);
        // int Y = rand() % Lf; // generate rand 0 - Lf-1
        float U = uniform_real_distribution<float>(0,1)(gen);
        // float U = (rand() % 100000) / 100000.0;
        float flag = FODnorm[Y]/M;
        int test = 1;
        vector<int> dirmax;
        for (int i =0;i<FODnorm.size();i++) {
            if (FODnorm[i] == fmax) dirmax.push_back(i);
        }
        direction = dirmax[0];
        while (U>flag || FOD[Y]<fmin){
            test++;
            if (test>1000){
                direction = dirmax[0];
                break;
            }
            Y = uniform_int_distribution<>(0,Lf-1)(gen);
            U = uniform_real_distribution<float>(0,1)(gen);
            flag = FODnorm[Y]/M;
            direction = Y;
        }
    }
    else direction = uniform_int_distribution<>(0,Lf-1)(gen);

    return direction;
}