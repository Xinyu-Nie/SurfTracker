#include "../header/InterpFODontoTRI.h"

vector<vector<float>> InterpFODontoTRI(NIImage& nii, vector<vector<float>>& coord){
	int N = nii.nim->nt;
    auto& img = nii.img;
	vector<vector<float>> FODonTri(coord.size(), vector<float>(N,0));
	for (int i = 0;i<coord.size();i++){
        vector<float>& qxyz = coord[i];// TRI coordinate
        int x = (int)qxyz[0], y = (int)qxyz[1], z = (int)qxyz[2];
        int c000[3] = {x,y,z}, c001[3] = {x,y,z+1}, c100[3] = {x+1,y,z}, c010[3] ={x,y+1,z},
            c110[3] ={x+1,y+1,z}, c101[3] = {x+1,y,z+1}, c011[3]={x,y+1,z+1}, c111[3] = {x+1,y+1,z+1};
        float dxyz[3];
        for (int iter = 0;iter<3;iter++){
            dxyz[iter] = (qxyz[iter] - c000[iter])/(c111[iter] - c000[iter]);
        }  
        // cout<<dxyz[0]<<' '<<dxyz[1]<<' '<<dxyz[2]<<endl;
        // cout<<img[c000[0]-1][c000[1]-1][c000[2]-1][0]<<endl;
        // vector<float> v000(N,0), v001(N,0), v010(N,0), v100(N,0), v110(N,0), v101(N,0), v011(N,0), v111(N,0);
        // if (c000[0]<1)
        for (int j=0;j<N;j++){
            FODonTri[i][j] = img[j][c000[2]-1][c000[1]-1][c000[0]-1] * (1-dxyz[0]) * (1-dxyz[1]) * (1-dxyz[2]) +
                             img[j][c100[2]-1][c100[1]-1][c100[0]-1] * (dxyz[0]) * (1-dxyz[1]) * (1-dxyz[2]) + 
                             img[j][c010[2]-1][c010[1]-1][c010[0]-1] * (1-dxyz[0]) * (dxyz[1]) * (1-dxyz[2]) + 
                             img[j][c001[2]-1][c001[1]-1][c001[0]-1] * (1-dxyz[0]) * (1-dxyz[1]) * (dxyz[2]) + 
                             img[j][c101[2]-1][c101[1]-1][c101[0]-1] * (dxyz[0]) * (1-dxyz[1]) * (dxyz[2]) + 
                             img[j][c110[2]-1][c110[1]-1][c110[0]-1] * (dxyz[0]) * (dxyz[1]) * (1-dxyz[2]) + 
                             img[j][c011[2]-1][c011[1]-1][c011[0]-1] * (1-dxyz[0]) * (dxyz[1]) * (dxyz[2]) + 
                             img[j][c111[2]-1][c111[1]-1][c111[0]-1] * (dxyz[0]) * (dxyz[1]) * (dxyz[2]);

            // cout<<"FODonTRI:" << FODonTri[i][j]<<endl;
        }
    }
    
	return FODonTri;
}
