#ifndef BUILD_FOR_WINDOWS
#include <unistd.h>
#include <sys/time.h>
#endif

#include"../header/Ufibertracking.h"
#include"../header/read_binary.h"
#include"../header/load_mesh.h"
#include"../header/Masktrimesh.h"
#include "../header/trackingonsufaceparallel.h"
#include <stdio.h>
#include <thread>
#include<algorithm>
#include<time.h>

/*
 * @param output output is the .trk file for tractogram
 * @param Mesh Mesh is the triangular mesh for half brain, the file is .obj
 * @param ROI ROI is the ROI to project the FOD data to track, it is a binary data
 * @param Sulcuslabel Sulcuslabel is a binary data to label the sulcus region for the Mesh, it
 * will intersect with the ROI binary data to generate the seed points
 * @param Gyruslabel1 Gyruslabel1 is a binary data to label the sulcus region for the Mesh, it
 * will intersect with the ROI binary data to generate the stopping points
 * @param Gyruslabel2 Gyruslabel2 is a binary data to label the sulcus region for the Mesh, it
 * will intersect with the ROI binary data to generate the stopping points
 * @param iteration iteration is the repeating times for each seed point
 * @param trackdegree Trackdegree is the upper bound degree for tracking 
 * @param Gyrusflag Gyrusflag is the flag for Gyruslabel, Gyrusflag=0, for Gyruslabel1 ~= Gyruslabel2;
 * Gyrusflag=0, for Gyruslabel1 ~= Gyruslabel2
 * @param Other_parameters Other input generated from the SurfaceFODpreprossing function
*/
void Ufibertracking(string output, string Mesh, string ROI, string Sulcuslabel,string Gyruslabel1, string Gyruslabel2, int iteration, vector<vector<float> >& coordS, vector<vector<float> >& trgS, vector<vector<float> >&TRIcen, vector<float>& THETA, vector<vector<float> >& FODProjected, vector<vector<float> >& trinerbor, vector<int>& dim, vector<float>& res, float minFOD, float trackdegree, int Gyrusflag){
    // srand(time(0));
    int numseed = 10000;
    vector<float> MTL;
    read_binary(ROI,MTL);
    vector<float> Gyrus1;
    read_binary(Gyruslabel1, Gyrus1);
    vector<float> Gyrus2;
    read_binary(Gyruslabel2, Gyrus2);
    vector<float> Suc;
    read_binary(Sulcuslabel, Suc);
 
    vector<vector<float> > coord, nS,trg;
    ReadObjShape(Mesh, coord, nS,trg);
    vector<float> SucdeepS;
    vector<vector<float> > tmp;
    vector<float> GY1;
    vector<float> GY2;
    vector<float> Faceori;
    Masktrimesh2(trg,Suc,MTL, SucdeepS, tmp, Faceori);
    Masktrimesh(trg,Gyrus1,MTL,GY1, tmp);
    Masktrimesh(trg,Gyrus2,MTL,GY2, tmp);
    int Ltri=trgS.size();
    vector<float> Gyrulabel1(Ltri,0);
    vector<float> Gyrulabel2 = Gyrulabel1;
    vector<float> SeedLabel(Ltri,0);
    vector<float> SeedP;
    for (int k = 0;k<Ltri;k++){
        auto triK = trgS[k];
        for (int m =0;m<3;m++){
            if (GY1[triK[m]] > 0){
                Gyrulabel1[k] = 1;
            }
            if (GY2[triK[m]] > 0){
                Gyrulabel2[k]=2;
            }
            if (SucdeepS[triK[m]]>0){
                SeedLabel[k] = 1;
            }

        }
    }

    for (int i = 0;i<SeedLabel.size();i++) if (SeedLabel[i]>0) SeedP.push_back(i);
    vector<float> seedflag;
    if (SeedP.size()<numseed){
        seedflag = SeedP;
    }
    else {
        auto flag = randperm(0,SeedP.size()-1, numseed);
        for (int i = 0;i<flag.size();i++) seedflag.push_back(SeedP[flag[i]]);
        sort(seedflag.begin(), seedflag.end());
    }
    int Nseed = seedflag.size();
  
    vector<vector<float> > Biflag(Nseed*iteration, vector<float>(2,0));
    vector<vector<vector<float> > > PtV(Nseed*iteration);
    // auto Tricoord = PtV;
    vector<vector<vector<float> > > Tricoord(Nseed*iteration);
    // auto PtV = cuvature;
    vector<float> Ltract(Nseed*iteration,0);
    auto flagend = Biflag;
    int badTrack = 0;
    thread task[100];
#ifdef BUILD_FOR_WINDOWS
    SYSTEM_INFO sysInfo;
    GetSystemInfo( &sysInfo );
    int num_cpu = sysInfo.dwNumberOfProcessors;
#else
    int num_cpu = sysconf(_SC_NPROCESSORS_ONLN);
#endif
    cout<<"cpu cores: "<< num_cpu<<endl;
    for (int m = 0;m<iteration;m++){
        int bm = m*Nseed;
        for (int t = 0;t<Nseed;t+=num_cpu){
            vector<vector<vector<float> > > tractline(num_cpu); vector<int> flagendinge(num_cpu), flagendingp(num_cpu); vector<vector<float> > biterminflag(num_cpu); vector<vector<float> > Tricount(num_cpu);
            for (int i=0;i<num_cpu && t+i<Nseed;i++){
                task[i] = thread(trackingonsufaceparallel, ref(tractline[i]), ref(flagendinge[i]), ref(flagendingp[i]), ref(biterminflag[i]), ref(coordS), ref(trgS), ref(TRIcen), ref(THETA), ref(FODProjected), seedflag[t+i], ref(trinerbor), ref(Gyrulabel2), ref(Gyrulabel1), minFOD, trackdegree);
            }
            for (int i=0;i<num_cpu && t+i<Nseed;i++) {
                task[i].join();
                Ltract[bm+t+i] = tractline[i].size();
                Biflag[bm+t+i][0] = biterminflag[i][0]; Biflag[bm+t+i][1] = biterminflag[i][1];
                PtV[bm+t+i] = tractline[i];
                flagend[bm+t+i] = {(float)flagendinge[i], (float)flagendingp[i]};
            }
            
        }
    }


    vector<float> LineInd;
    vector<vector<float>> PtVec;
    vector<int> Lineposition;
    if (Gyrusflag == 1){
        int flag = 0;
        int att = 0;// line position
        for (int k = 0;k<Ltract.size();k++){
            auto lu = Ltract[k];
            if (lu>20){
                if (Biflag[k][0] == 4 || Biflag[k][0]==4.5){
                    if (Biflag[k][1] == 4 || Biflag[k][1] == 4.5){
                        if (flag>0){
                            LineInd.push_back(LineInd.back()+lu);
                            PtVec.insert(PtVec.end(), PtV[k].begin(), PtV[k].end());
                            Lineposition.push_back(att);
                        }
                        else {
                            LineInd.push_back(lu);
                            PtVec = PtV[k];
                            flag = 1;
                            Lineposition.push_back(att);
                        }
                    }
                }
            }
            att++;
        }
    }
    else if (Gyrusflag == 0){
        int flag = 0;
        int att = 0;// line position
        for (int k = 0;k<Ltract.size();k++){
            auto lu = Ltract[k];
            if (lu>20){
                if (Biflag[k][0] == 4 || Biflag[k][0]==4.5){
                    if (Biflag[k][1] == 4 || Biflag[k][1] == 4.5){
                        if (Biflag[k][0]!=Biflag[k][1]){
                            if (flag>0){
                                LineInd.push_back(LineInd.back()+lu);
                                PtVec.insert(PtVec.end(), PtV[k].begin(), PtV[k].end());
                                Lineposition.push_back(att);
                            }
                            else {
                                LineInd.push_back(lu);
                                PtVec = PtV[k];
                                flag = 1;
                                Lineposition.push_back(att);
                            }
                        }
                    }
                }
            }
            att++;
        }
    }
    
    for (int i=0;i< PtVec.size();i++){
        for (int j =0;j<3;j++){
            PtVec[i][j] = (PtVec[i][j] - 1) * res[j];
        }
    }
    Write2File(output, LineInd, PtVec, dim, res);

}   
/*
%PtVec: points, assume in voxel xyz coordinate
% dim: 1x3 array, size of the image array
% voxel_size: 1x3 array, voxel scale
*/
void Write2File(string file, vector<float> &LineInd, vector<vector<float>> &PtVec, vector<int>& dim, vector<float>& voxel_size){
    FILE* pfile = fopen(file.c_str(), "wb");
    char id_string[] = {84, 82, 65, 67, 75, 0};
    fwrite(id_string, sizeof(char), sizeof(id_string), pfile);
    vector<int16_t> dim_16(dim.size());
    for (int i=0;i<dim.size();i++) dim_16[i] = (int16_t) dim[i];
    fwrite(&dim_16[0], sizeof(int16_t), 3, pfile);

    fwrite(&voxel_size[0], sizeof(float), voxel_size.size(), pfile);

    float origin[] = {0,0,0};
    fwrite(origin, sizeof(float), sizeof(origin)/sizeof(float), pfile);
    int16_t n_s = PtVec[0].size() - 3; // colors are stored as scalars
    fwrite(&n_s, sizeof(int16_t), 1, pfile);
    vector<char> sname(200,0);
    fwrite(&sname[0], sizeof(char), sname.size(), pfile);
    int16_t n_p = 0;
    fwrite(&n_p, sizeof(int16_t), 1, pfile);
    vector<char> pname(200,0);
    fwrite(&pname[0], sizeof(char), pname.size(), pfile);
    for (int i =0;i<4;i++){
        for (int j=0;j<4;j++){
            float x;
            if (i == j) x = 1;
            else x = 0;
            fwrite(&x, sizeof(float), 1, pfile);
        }
    }
    vector<char> reserved(444,0);
    fwrite(&reserved[0], sizeof(char), reserved.size(), pfile);
    char voxel_order[] = {'R', 'A', 'S', 0};
    fwrite(voxel_order, sizeof(char), sizeof(voxel_order), pfile);
    char pad2[] = {76, 80, 83, 0};
    fwrite(pad2, sizeof(char), sizeof(pad2), pfile);
    float image_orientation_patient[] = {1, 0, 0, 0, 1, 0};
    fwrite(image_orientation_patient, sizeof(float), sizeof(image_orientation_patient)/sizeof(float), pfile);
    char pad1[] = {0, 0};
    fwrite(pad1, sizeof(char), sizeof(pad1), pfile);
    vector<unsigned char> invert_xyz(3,0);
    fwrite(&invert_xyz[0], sizeof(unsigned char), invert_xyz.size(), pfile);
    vector<unsigned char> swap_xyz(3,0);
    fwrite(&swap_xyz[0], sizeof(unsigned char), swap_xyz.size(), pfile);
    int n_count = LineInd.size();
    fwrite(&n_count, sizeof(int32_t), 1, pfile);
    int version = 1;
    fwrite(&version, sizeof(int32_t), 1, pfile);
    int hdr_size = 1000;
    fwrite(&hdr_size, sizeof(int32_t), 1, pfile);

    // int count = 0;
    for (int i = 0;i<LineInd.size();i++){
        int num_pt;
        if (i==0) num_pt = LineInd[0];
        else num_pt = LineInd[i] - LineInd[i-1];
        fwrite(&num_pt, sizeof(int32_t), 1, pfile);
        for (int j = LineInd[i] - num_pt; j<LineInd[i];j++){
            fwrite(&PtVec[j][0], sizeof(float), PtVec[0].size(),pfile);
        }
    }
    fclose(pfile);
    cout<< "output file is "<<file<<endl;
}