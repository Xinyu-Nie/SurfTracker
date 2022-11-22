#include<iostream>
#include "SurfaceFODPro.h"
#include <unistd.h>
#include <sys/time.h>
using namespace std;
int SPHARM_Order;
int k1,k2;
float fai = PI/90;
float theta = PI/90;
/*
 * @param FOD the FOD data
 * @param Mesh the triangular mesh for half brain, the file is obj.\
 * @param ROI the ROI to project the FOD data to track, it is a binary data
 * @param output the folder saving computed FOD data on the triangular mesh
 * @param shrinkdis the distace to shrink the cortical surface into the WM
*/
void SurfaceFODPro(string& FOD, string& Mesh, string& ROI, string& output, float shrinkdis){

	// load nii
	NIImage nii(FOD);
	int N = nii.nim->nt;
	
	vector<vector<float>> coord, nS, trg;
 	ReadObjShape(Mesh, coord, nS, trg);
	/**************
	 warning trg index start from 1
	****************/
	// read ROI 
    // read raw binary ROI
	ifstream inFile(ROI,ios::in|ios::binary);
	vector<float> MTL;
	int s;
	while(inFile.read((char *)&s, sizeof(int))) { //一直读到文件结束
        MTL.push_back((float)s);
    }
	inFile.close();
	for (int i = 0;i<coord.size();i++){
		for (int j = 0;j<coord[0].size();j++) coord[i][j] -= shrinkdis*nS[i][j];
	}

	for (int j = 0;j<3;j++){
		for (int i = 0;i<coord.size();i++){
			coord[i][j] = coord[i][j]/nii.nim->pixdim[j+1] +1;
		}
	}
	vector<vector<float>> coordS, trgS;
	Masktrimesh(trg, coord, MTL, coordS, trgS);

	int L1 = trgS.size();

	vector<float> FAI;
	for (float i = PI/180;i<=PI;i+=fai){
		FAI.push_back(i);
	}
	vector<float> THETA;
	for (float i = PI/180; i<=2*PI;i+=theta){
		THETA.push_back(i);
	}
	k1 = FAI.size(); k2 = THETA.size();
	vector<float> weight(k1,0);
	for (int i = 0;i<k1;i++) weight[i] = sin(FAI[i]);
	int K = k1*k2;
	
	int m = trgS.size(), n = trgS[0].size();
	vector<vector<float>> TRIcen(m,vector<float>(n,0)), Normal(m,vector<float>(n,0));
	vector<float> flag(L1,1);
	for (int t = 0;t<L1;t++){
		for (int j = 0;j<TRIcen[0].size();j++){
			TRIcen[t][j] = (coordS[trgS[t][0]][j] + coordS[trgS[t][1]][j] + coordS[trgS[t][2]][j])/3;// correct
		}
	}
	vector<vector<float>> x = InterpFODontoTRI(nii, TRIcen); // correct x(33824 * 153)
	SPHARM_Order = (sqrt(8*N+1)-3)/2;
	vector<vector<float>> FODProjected(k2, vector<float>(L1,0));
    vector<vector<float>> FODcoordinate(3, vector<float>(K,0));
    for (int h2 = 0;h2<k2;h2++){
        float R2[3][3];
        float t[3] = {0,0,1};
        rotate(R2, t, THETA[h2]);
        float dr2[3];
		t[0]=1;t[1]=0;t[2]=0;
        for (int i = 0;i<3;i++){
            dr2[i] = dot(R2[i], t);
        }
		t[0]=0;t[1]=0;t[2]=1;
        for (int h1=0;h1<k1;h1++){
            float Nvv[3];
            cross(Nvv, t, dr2);
            normalize(Nvv);
            float R1[3][3];
            rotate(R1, Nvv, FAI[h1]);
            float dr1[3];
			t[0]=0;t[1]=0;t[2]=1;
            for (int i = 0;i<3;i++){
                dr1[i] = dot(R1[i], t);
            }
            normalize(dr1);
            for (int i = 0;i<FODcoordinate.size();i++){
                FODcoordinate[i][(h2)*k1+h1] = dr1[i];
            }
        }
    }
    vector<vector<float>> BS;
    matrixB(BS, FODcoordinate, SPHARM_Order);
    auto R2C=Real2complex((int)SPHARM_Order);
    auto C2R=Complex2real(SPHARM_Order);
	vector<vector<complex<double>>> C2R1(C2R[0].size(),vector<complex<double>>(C2R.size(), 0));
	vector<vector<complex<double>>> R2C1(R2C[0].size(),vector<complex<double>>(R2C.size(), 0));
	for (int i=0;i<C2R1.size();i++){
		for (int j=0;j<C2R1[0].size(); j++){
			C2R1[i][j] = C2R[j][i];
		}
	}
	for (int i =0;i<R2C1.size();i++){
		for (int j=0;j<R2C1[0].size();j++){
			R2C1[i][j] = R2C[j][i];
		}
	}
	float thetay=PI/2;
	auto UY90p=SpYrotation(SPHARM_Order, thetay);
	thetay = -PI/2;
	auto UY90n=SpYrotation(SPHARM_Order,thetay);
	// multi thread
	int num_cpu = sysconf(_SC_NPROCESSORS_ONLN);
	// int num_cpu = sysInfo.dwNumberOfProcessors;
	printf("system enable cpu num is %d\n", num_cpu);
	int n_per_core = FODProjected[0].size()/num_cpu + 1;
	thread task[100];
	
	// task[0] = thread(SurfaceThread, 0, 0, ref(FODProjected), ref(coordS), ref(trgS), ref(x), ref(C2R), ref(C2R1), ref(R2C), ref(R2C1), ref(UY90p), ref(UY90n), ref(BS), ref(weight));
	// task[0].join();
	timeval start, end;
	gettimeofday(&start, NULL);
	for (int i = 0;i<num_cpu;i++){
		if (i==num_cpu-1){
			task[i] = thread(SurfaceThread, i*n_per_core, FODProjected[0].size()-1, ref(FODProjected), ref(coordS), ref(trgS), ref(x), ref(C2R), ref(C2R1), ref(R2C), ref(R2C1), ref(UY90p), ref(UY90n), ref(BS), ref(weight));
		}
		else {
			task[i] = thread(SurfaceThread, i*n_per_core, (i+1)*n_per_core-1, ref(FODProjected), ref(coordS), ref(trgS), ref(x), ref(C2R), ref(C2R1), ref(R2C), ref(R2C1), ref(UY90p), ref(UY90n), ref(BS), ref(weight));
		}
	}
	for (int i=0;i<num_cpu;i++) task[i].join();
	gettimeofday(&end, NULL);
	double delta_t = (end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec)/1000000.0;
	// cout<<"time = "<<delta_t<<"s"<<endl;
	// save(output, FODProjected);
	vector<vector<float>> Vert2FaceMap;
	vector<vector<float>> TrgNbrMap;
	BuildMeshNbrhood(Vert2FaceMap, TrgNbrMap, coordS, trgS); // correct , TrgNbrMap trgS index from 0, -1 means nan
	auto trinerbor = TrgNbrMap;
	for (int i =0;i<trinerbor.size();i++){
		trinerbor[i][0] = TrgNbrMap[i][2];
		trinerbor[i][2] = TrgNbrMap[i][0];
	}
	gettimeofday(&start, NULL);
	/*
	save fodsphere data into binary file
	*/
	FILE* pfile = fopen(output.c_str(), "wb");
	binary_save(pfile, trinerbor);
	binary_save(pfile, trgS);
	binary_save(pfile,coordS);
	binary_save(pfile, TRIcen);
	binary_save(pfile,THETA);
	binary_save(pfile,FODProjected);
	fwrite(&nii.nim->dim[1], sizeof(int), 4,pfile);
	fwrite(&nii.nim->pixdim[1], sizeof(float), 3, pfile);
	fclose(pfile);
	cout<<"saved to "+ output<<endl;
	gettimeofday(&end, NULL);
	delta_t = (end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec)/1000000.0;
	// cout<<"time = "<<delta_t<<"s"<<endl;
}
void SurfaceThread(int t1, int t2, vector<vector<float>>& FODProjected, vector<vector<float>>& coordS, vector<vector<float>>& trgS, vector<vector<float>>& x, vector<vector<complex<double>>>& C2R, vector<vector<complex<double>>>& C2R1, vector<vector<complex<double>>>& R2C, vector<vector<complex<double>>>& R2C1, vector<vector<float>>& UY90p, vector<vector<float>>& UY90n,vector<vector<float>>& BS, vector<float>& weight){
	for (int t = t1;t<=t2;t++){
		float v1[3], v2[3];
		for (int j = 0;j<coordS[0].size();j++){
			v1[j] = coordS[trgS[t][1]][j] - coordS[trgS[t][0]][j];
			v2[j] = coordS[trgS[t][2]][j] - coordS[trgS[t][0]][j];
		}
		float Nv[3];
		cross(Nv, v1, v2);
		float q = norm(Nv);
		if (q>0){
			Nv[0] /= q;Nv[1] /= q;Nv[2] /= q;
			normalize(v1);
			float vy[3];
			cross(vy,Nv,v1);
			auto qy = norm(vy);
			if (qy>0){
				vy[0]/=qy;vy[1]/=qy;vy[2]/=qy;
			}
			float Rt[3][3] = {{v1[0],v1[1],v1[2]}, {vy[0],vy[1],vy[2]},{Nv[0],Nv[1],Nv[2]}};
			auto tmp = DecomrotationPro(Rt);
			float alfa = tmp[0], beta = tmp[1], gamma = tmp[2];
			auto Ylr = x[t];
			vector<complex<double>> Ylc(C2R[0].size(), 0); // Ylc = C2R1*Ylr, wrong

			for (int i=0;i<C2R1.size();i++){
				for (int j = 0;j<C2R1[0].size();j++){
					Ylc[i]+=C2R1[i][j]*(double)Ylr[j];
				}
			}
			auto Ualfa = SpZrotation((int)SPHARM_Order, alfa+PI/2);
			vector<complex<double>> Ylc1(Ualfa.size(),0); // Ylc1 wrong
			for (int i = 0;i<Ualfa.size();i++){
				for (int j = 0;j<Ualfa[0].size();j++){
					Ylc1[i]+=Ualfa[i][j]*Ylc[j];
				}
			}
			auto UZbeta = SpZrotation((int)SPHARM_Order, beta);
			int n = UY90n.size();
			// Ylc2=UY90n*UZbeta*UY90p*Ylc1;
			vector<complex<double>> tmp1(n,0), tmp2(n,0);
			for (int i = 0;i<n;i++) {
				for (int j=0;j<n;j++){
					tmp1[i]+=(double)UY90p[i][j]*Ylc1[j];
				}
			}
			for (int i = 0;i<n;i++) {
				for (int j=0;j<n;j++){
					tmp2[i]+=UZbeta[i][j]*tmp1[j];
				}
			}
			vector<complex<double>> Ylc2(n,0);
			for (int i = 0;i<n;i++) {
				for (int j=0;j<n;j++){
					Ylc2[i]+=(double)UY90n[i][j]*tmp2[j];
				}
			}
			auto Ugamma = SpZrotation((int)SPHARM_Order, gamma-PI/2);
			vector<complex<double>> Ylc3(Ugamma.size(),0);
			for (int i = 0;i<Ugamma.size();i++){
				for (int j = 0;j<Ugamma[0].size();j++){
					Ylc3[i]+=Ugamma[i][j]*Ylc2[j];
				}
			}
			vector<float> Ylr3(R2C1.size(),0);
			for (int i=0;i<R2C1.size();i++){
				complex<double> sum = 0;
				for (int j = 0;j<R2C1[0].size();j++){
					sum+=R2C1[i][j]*Ylc3[j];
				}
				Ylr3[i] = real(sum);
			}
			vector<float> FODPROJ(BS.size(),0);
			for (int i=0;i<BS.size();i++){
				for (int j = 0;j<BS[0].size();j++) FODPROJ[i]+=BS[i][j]*Ylr3[j];
				FODPROJ[i] = max(FODPROJ[i],(float)0);
			}
			for (int h = 0;h<k2;h++){
				vector<float> fodsphere(FODPROJ.begin()+h*k1, FODPROJ.begin()+(h+1)*k1);
				for (int i =0;i<weight.size();i++) FODProjected[h][t] += weight[i]*fodsphere[i];
				FODProjected[h][t]*=fai;
			}
		}
	
		
	}
}