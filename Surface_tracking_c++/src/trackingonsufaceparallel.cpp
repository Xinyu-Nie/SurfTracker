#include "../header/trackingonsufaceparallel.h"

void track(lsqrDense& lsqr, int a, int& cpedge, int& cnedge, vector<vector<float> >& tractline, vector<vector<float> >& tractbidirection, float* vd, vector<float>&v0, vector<float>& v1, vector<float>& v2, float seedflag){
    
            double AA[6];
            double* A1[3];
            A1[0] = &AA[0], A1[1] = &AA[2], A1[2] = &AA[4];
            A1[0][0] = vd[0], A1[0][1] = v1[0] - v2[0];
            A1[1][0] = vd[1], A1[1][1] = v1[1] - v2[1];
            A1[2][0] = vd[2], A1[2][1] = v1[2] - v2[2];
            double b1[3] = {v1[0]-v0[0], v1[1]-v0[1], v1[2]-v0[2]};
            // lsqr
            lsqr.SetMatrix(A1);
            double x1[2];
            lsqr.Solve(3,2,b1,x1);
            // delete(A1);
            if (x1[1]>=0 && x1[1]<=1){
                //cout<<x1[1]<<' ';
                if (x1[0]>0){
                    tractline.push_back(vector<float>(3,0));
                    for (int j = 0;j<3;j++){
                        tractline.back()[j] = x1[1]*v2[j] + (1-x1[1])*v1[j];
                    }
                    cpedge = a;
                    // BaryP = vector<float>(4);
                }
                else {
                    tractbidirection.push_back(vector<float>(3,0));
                    for (int j = 0;j<3;j++){
                        tractbidirection.back()[j] = x1[1]*v2[j] + (1-x1[1])*v1[j];
                    }
                    cnedge = a;
                }
            }
        
}

void trackingonsufaceparallel(vector<vector<float> >& tract, int& flagendinge, int& flagendingp, vector<float>& biterminflag, vector<vector<float> >& coordS, vector<vector<float> >& trgS, vector<vector<float> >& tricenter, vector<float>& THETA, vector<vector<float> >& FOD, float seedflag, vector<vector<float> >& trinerbor,vector<float>& Triendflag,vector<float>& Tribeginflag, float minFOD, float trackdegree){
    auto seedpoint = tricenter[(int)seedflag];
    // seedFOD=FOD(:,seedflag);?
    vector<float> seedFOD(FOD.size(),0);
    for (int i =0;i<FOD.size();i++) seedFOD[i] = FOD[i][(int) seedflag];
    // auto tractline = seedpoint;
    vector<vector<float> > tractline;
    tractline.push_back(seedpoint);
    vector<vector<float> > tractbidirection;
    tractbidirection.push_back(seedpoint);
    auto alfa = rejectsampling(seedFOD, 2, minFOD);
    // alfa = 66;
    // 135 15 67
    // 1 8 9
    // int alfa = 134;
    auto v1 = coordS[trgS[seedflag][0]];
    auto v2 = coordS[trgS[seedflag][1]];
    auto v3 = coordS[trgS[seedflag][2]];
    auto v0 = seedpoint;
    int cpedge = -1;
    int cnedge = -1;
    float v12[3], v13[3], v23[3];
    for (int i=0;i<3;i++){
        v12[i] = v2[i] - v1[i];
        v13[i] = v3[i] - v1[i];
        v23[i] = v3[i] - v2[i];
    }
    float Nv[3];
    cross(Nv, v12,v13);
    float q = norm(Nv);
    lsqrDense lsqr;
    lsqr.SetMaximumNumberOfIterations( 20 );
    if (q>0){
        for (int i =0;i<3;i++) Nv[i] = Nv[i]/q;
        normalize(v12);
        float R[3][3];
        rotate(R, Nv, THETA[alfa]);
        float vd[3];
        for (int i=0;i<3;i++) vd[i] = dot(R[i], v12);

        float out[3];
        cross(out, vd, v12);
        float a1 = dot(out, out);
        cross(out, vd, v13);
        float a2 = dot(out, out);
        cross(out, vd, v23);
        float a3 = dot(out, out);
        if (a1 > 0){
            track(lsqr, 1, cpedge, cnedge, tractline, tractbidirection, vd, v0, v1, v2, seedflag);
        }
        if (a2>0){
            track(lsqr, 2, cpedge, cnedge, tractline, tractbidirection, vd, v0, v1, v3, seedflag);
        }
        if (a3>0){
            track(lsqr, 3, cpedge, cnedge, tractline, tractbidirection, vd, v0, v2, v3, seedflag);
        }
        int flagpd = 0;
        int flagnd = 0;
        auto currentnegativetrig = seedflag;
        auto currentpositivetrig = seedflag;
        vector<float> Ptri, Ntri;
        Ptri.push_back(seedflag);
        Ntri.push_back(seedflag);
        biterminflag = vector<float>(2,0);
        int flap =1;
        flagendingp = 0;
        int terminflag = 0;
        if (cpedge!=-1) {
            while (flagpd==0 && cpedge!=0 && flap<1000 && terminflag ==0){
                auto Tric = trinerbor[currentpositivetrig];
                if (Tric[cpedge-1]>=0){
                    if (Triendflag[Tric[cpedge-1]]>0){
                        flagendingp = 1;
                        biterminflag[0] = 4;
                        break;
                    }
                    if (Tribeginflag[Tric[cpedge-1]]>0){
                        flagendingp = 2;
                        biterminflag[0] = 4.5;
                        break;
                    }
                    vector<float> outpoint; int outedgef;
                    traveltrigparallel(outpoint, outedgef, terminflag, coordS, trgS, tractline.back(),tractline[tractline.size()-2], THETA, FOD, (int)Tric[cpedge-1], (int)currentpositivetrig, cpedge, minFOD, trackdegree);
                    if (terminflag == 2) return;
                    currentpositivetrig = Tric[cpedge-1];
                    cpedge = outedgef;
                    if (outedgef > 0){
                        tractline.push_back(outpoint);
                        Ptri.push_back(currentpositivetrig);
                        flap++;
                    }
                }
                else{
                    flagpd = 1;
                    biterminflag[0] = 5;
                }
            }
            if (biterminflag[0]==0){
                if (terminflag>0){
                    biterminflag[0] = terminflag;
                }
                else if (flap==1000) biterminflag[0] = 6;
            }
        }
        else {
            biterminflag[0] = 7;
        }
        
        int flae = 1;
        flagendinge = 0;
        terminflag = 0;
        if (cnedge!=-1){
            while (flagnd==0 && cnedge!=0 && flae<1000 && terminflag ==0){
                auto Tric = trinerbor[currentnegativetrig];
                if (Tric[cnedge-1]>=0){
                    if (Triendflag[Tric[cnedge-1]]>0){
                        flagendinge = 1;
                        biterminflag[1] = 4;
                        break;
                    }
                    if (Tribeginflag[Tric[cnedge-1]]>0){
                        flagendinge = 2;
                        biterminflag[1] = 4.5;
                        break;
                    }
                    vector<float> outpoint; int outedgef;
                    traveltrigparallel(outpoint, outedgef, terminflag,coordS, trgS, tractbidirection.back(),tractbidirection[tractbidirection.size()-2], THETA, FOD, (int)Tric[cnedge-1], (int)currentnegativetrig, cnedge, minFOD, trackdegree);
                    if (terminflag == 2) return;
                    currentnegativetrig = Tric[cnedge-1];
                    cnedge = outedgef;
                    if (outedgef > 0){
                        tractbidirection.push_back(outpoint);
                        Ntri.push_back(currentnegativetrig);
                        flae++;
                    }
                }
                else{
                    flagnd = 1;
                    biterminflag[1] = 5;
                }
            }
            if (biterminflag[1]==0){
                if (terminflag>0){
                    biterminflag[1] = terminflag;
                }
                else if (flae==1000) biterminflag[1] = 6;
            }
        }
        else biterminflag[1] = 7;  

        for (int i=tractbidirection.size()-1;i>=0;i--) tract.push_back(tractbidirection[i]);
        int L = tractline.size();
        if (L>1){
            for (int i = 1;i<L;i++){
                tract.push_back(tractline[i]);
            }
        }
     
        
    }
    return;
    }
