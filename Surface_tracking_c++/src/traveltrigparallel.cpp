// function [outpoint,outedgef,terminflag] = traveltrigparallel(coordS,trgS,inpoint,prepoint,THETA,FOD,intri,pretri,inedge,minFOD,trackdegree)

#include "../header/traveltrigparallel.h"
void track(int a, lsqrDense& lsqr, float* vd, vector<float>&v0, vector<float>& vin1, vector<float>& vin2, vector<float>& outpoint, int &outedgef, int intri){
    
    double AA[6];
    double* A1[3];
    A1[0] = &AA[0], A1[1] = &AA[2], A1[2] = &AA[4];
    A1[0][0] = vd[0], A1[0][1] = vin1[0] - vin2[0];
    A1[1][0] = vd[1], A1[1][1] = vin1[1] - vin2[1];
    A1[2][0] = vd[2], A1[2][1] = vin1[2] - vin2[2];
    double b1[3] = {vin1[0]-v0[0], vin1[1]-v0[1], vin1[2]-v0[2]};
    // lsqr
    lsqr.SetMatrix(A1);
    double x1[2];
    lsqr.Solve(3,2,b1,x1);
   
   
    if (x1[1]>=0 && x1[1]<=1 && x1[0]>0){
            for (int j = 0;j<3;j++){
                outpoint[j] = x1[1]*vin2[j] + (1-x1[1])*vin1[j];
            }
            outedgef = a;
           
    }
        
}

void traveltrigparallel(vector<float>& outpoint, int& outedgef, int& terminflag, vector<vector<float> >& coordS, vector<vector<float> >& trgS, vector<float>& inpoint, vector<float>& prepoint, vector<float>& THETA, vector<vector<float> >& FOD, int intri, int pretri, int inedge, float minFOD, float trackdegree){
    float threshold = cos(trackdegree*PI/180);
    float toler = 1e-4;
    terminflag = 3;
    auto vin1 = coordS[trgS[intri][0]];
    auto vin2 = coordS[trgS[intri][1]];
    auto vin3 = coordS[trgS[intri][2]];
    float v12[3], v13[3], v23[3];
    for (int i =0;i<3;i++) {
        v12[i] = vin2[i] - vin1[i];
        v13[i] = vin3[i] - vin1[i];
        v23[i] = vin3[i] - vin2[i];
    }
    normalize(v12);
    normalize(v13);
    normalize(v23);
    float e[3][2] = 
    {
        {trgS[intri][0], trgS[intri][1]},
        {trgS[intri][2], trgS[intri][0]},
        {trgS[intri][1], trgS[intri][2]}
    };
    float Nvin[3];
    cross(Nvin, v12, v13);
    normalize(Nvin);
    auto v0 = inpoint;
    auto vpre1 = coordS[trgS[pretri][0]];
    auto vpre2 = coordS[trgS[pretri][1]];
    auto vpre3 = coordS[trgS[pretri][2]];
    float v12pre[3], v13pre[3], v23pre[3];
    for (int i =0;i<3;i++) {
        v12pre[i] = vpre2[i] - vpre1[i];
        v13pre[i] = vpre3[i] - vpre1[i];
        v23pre[i] = vpre3[i] - vpre2[i];
    }
    normalize(v12pre);
    normalize(v13pre);
    normalize(v23pre);
    float Nvpre[3];
    cross(Nvpre, v12pre, v13pre);
    normalize(Nvpre);
    float edge[2];
    if (inedge == 1){
        edge[0] = trgS[pretri][0]; edge[1]=trgS[pretri][1];
    }
    else if (inedge == 2){
        edge[0] = trgS[pretri][2]; edge[1]=trgS[pretri][0];
    }
    else if (inedge == 3){
        edge[0] = trgS[pretri][1]; edge[1]=trgS[pretri][2];
    }
    int flage[3] = {0};
    float* em;
    int m;
    for (m=0;m<3;m++){
        em = e[m];
        vector<float> u;
        if (em[0] == edge[0] || em[0] == edge[1]) u.push_back(em[0]);
        if (em[1] == edge[0] || em[1] == edge[1]) u.push_back(em[1]);
        if (u.size()>1.5) {flage[m] = 1;
            break;}
    }
    
    if (sum(flage) == 0) {
        cout<<"neighbor error";
        terminflag = 2;
        return;
    }
    if (em[0] == edge[0] && em[1] == edge[1]) {
        Nvpre[0] = -Nvpre[0]; Nvpre[1] = -Nvpre[1]; Nvpre[2] = -Nvpre[2]; }
        
    // Rotation matrix;
    float direction1[3];
    cross(direction1, Nvpre, Nvin);
    auto Len = norm(direction1);
    float ux = direction1[0]/Len, uy = direction1[1]/Len, uz = direction1[2]/Len;
    float cth = dot(Nvpre, Nvin) / sqrt(dot(Nvin, Nvin)*dot(Nvpre, Nvpre) );
    float sth = Len / sqrt(dot(Nvin, Nvin)*dot(Nvpre, Nvpre) );
    // auto R = rotate(direction1, cth, sth);
    float R[3][3];
    R[0][0] = cth+ux*ux*(1-cth); R[0][1] = ux*uy*(1-cth)-uz*sth; R[0][2] = ux*uz*(1-cth)+uy*sth;
	R[1][0] = uy*ux*(1-cth)+uz*sth; R[1][1] = cth+(uy*uy)*(1-cth); R[1][2] = uy*uz*(1-cth)-ux*sth;
	R[2][0] = uz*ux*(1-cth)-uy*sth; R[2][1] = uz*uy*(1-cth)+ux*sth; R[2][2] = cth+(uz*uz)*(1-cth);

    auto vfix = coordS[edge[0]];
    float rvector[3];
    for (int i = 0;i< 3;i++) rvector[i] = prepoint[i] - vfix[i];
    float invector[3];
    float rvector1[3];
    for (int i=0;i< 3;i++) rvector1[i] = dot(R[i], rvector);
    float pret[3];
    for (int i =0;i<3;i++) pret[i] = vfix[i] + rvector1[i];
    for (int i = 0;i<3;i++) invector[i] = inpoint[i] - pret[i];
    float qinv = norm(invector);
    if (qinv > 1e-10) {invector[0]/=qinv; invector[1]/=qinv; invector[2]/=qinv;}
    if (dot(invector, Nvin) -1e-3>0){
        // cout<<"invector rotation failed"<<endl;
        terminflag = 2;
        return;
    }

    int N = 1000;
    int tn = 0;
    outedgef = -1;
    outpoint = inpoint;
    while (terminflag == 3){
        tn++;
        if (tn>N){
            terminflag = 1;
            break;
        }
        vector<float> cFOD(FOD.size(), 0);
        for (int i = 0 ;i<FOD.size();i++) cFOD[i] = FOD[i][intri];
        auto alfa = rejectsampling(cFOD, 1, minFOD);
        // cout<<alfa<<endl;
        // int alfa = 25;
        float R12[3][3];
        rotate(R12, Nvin, THETA[alfa]);
        float vd[3];
        for (int i =0;i<3;i++) vd[i] = dot(R12[i], v12);
        float qd = norm(vd);
        if (qd>0) {vd[0]/=qd; vd[1]/=qd; vd[2]/=qd;}
        if (dot(vd, Nvin) > toler) {cout<<"vector rotation failed"<<endl; break;}
        float out[3];
        cross(out, vd, v12);
        float a1 = dot(out, out);
        cross(out, vd, v13);
        float a2 = dot(out, out);
        cross(out, vd, v23);
        float a3 = dot(out, out);
        lsqrDense lsqr;
        lsqr.SetMaximumNumberOfIterations( 20 );
        // float outpoint[3];
        // int outedgef=-1;
        double x1[2];
        if (a1 > toler && m!=0) {
            track(1, lsqr, vd, v0, vin1, vin2, outpoint, outedgef, intri);
        }
        if (a2 > toler && m!=1){
            track(2, lsqr, vd, v0, vin1, vin3, outpoint, outedgef, intri);
        }
        if (a3 > toler && m!=2){
            track(3, lsqr, vd, v0, vin2, vin3, outpoint, outedgef, intri);
        }
        if (outedgef!=-1){
            float vout[3];
            vout[0] = outpoint[0] - inpoint[0]; vout[1] = outpoint[1] - inpoint[1]; vout[2] = outpoint[2] - inpoint[2];
            normalize(vout);
            if (dot(vout,invector) > threshold){
                terminflag = 0;
                break;
            }
        }
    }

    if (outedgef == -1) outedgef = 0;
    

}

