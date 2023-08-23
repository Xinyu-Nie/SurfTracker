#include "../header/BuildMeshNbrhood.h"
void BuildMeshNbrhood(vector<vector<float>>& Vert2FaceMap, vector<vector<float>>& TrgNbrMap,vector<vector<float>>& coord, vector<vector<float>>& trg){
    Vert2FaceMap = vector<vector<float>>(coord.size(), vector<float>(20,0));
    for (int i = 0;i<coord.size();i++) Vert2FaceMap[i][0] = 1;
    for (int i =0;i<trg.size();i++){
        for (int j = 0;j<3;j++){
            Vert2FaceMap[trg[i][j]][0] ++;
            Vert2FaceMap[trg[i][j]][Vert2FaceMap[trg[i][j]][0]-1] = i;
        }
    }
    TrgNbrMap = vector<vector<float>>(trg.size(), vector<float>(3,-1));
    
    for (int i =0;i<trg.size();i++){
        // %faces sharing vert1 and 2
        for (int j =1;j<Vert2FaceMap[trg[i][0]][0];j++){
            for (int k = 1;k<Vert2FaceMap[trg[i][1]][0];k++){
                if (Vert2FaceMap[trg[i][0]][j]!=i && Vert2FaceMap[trg[i][0]][j] == Vert2FaceMap[trg[i][1]][k]){
                    TrgNbrMap[i][2] = Vert2FaceMap[trg[i][0]][j];
                    break;
                }
            }
            if (TrgNbrMap[i][2]>-1) break;
        }

        // %faces sharing vert1 and 3
        for (int j =1;j<Vert2FaceMap[trg[i][0]][0];j++){
            for (int k = 1;k<Vert2FaceMap[trg[i][2]][0];k++){
                if (Vert2FaceMap[trg[i][0]][j]!=i && Vert2FaceMap[trg[i][0]][j] == Vert2FaceMap[trg[i][2]][k]){
                    TrgNbrMap[i][1] = Vert2FaceMap[trg[i][0]][j];
                    break;
                }
            }
            if (TrgNbrMap[i][1]>-1) break;
        }

        // %faces sharing vert2 and 3
        for (int j =1;j<Vert2FaceMap[trg[i][1]][0];j++){
            for (int k = 1;k<Vert2FaceMap[trg[i][2]][0];k++){
                if (Vert2FaceMap[trg[i][1]][j]!=i && Vert2FaceMap[trg[i][1]][j] == Vert2FaceMap[trg[i][2]][k]){
                    TrgNbrMap[i][0] = Vert2FaceMap[trg[i][1]][j];
                    break;
                }
            }
            if (TrgNbrMap[i][0]>-1) break;
        }
    }
}