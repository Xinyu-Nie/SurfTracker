#include "../header/load_mesh.h"
void ReadObjShape(string fname, vector<vector<float> >& pt_coordinates, vector<vector<float> >& pt_normals, vector<vector<float> >& triangles){

    // char path[] = "C:\\Users\\SAM\\trekker\\data\\lh_white_DTISpace.obj";
    FILE* fid = fopen(fname.c_str(), "r");
    // char buf[sizeof(short)];
    char shape_type[10];
    // vector<vector<float>> pt_coordinates;
    // vector<vector<float>> pt_normals;
    // vector<vector<float>> pt_colors;
    // vector<vector<float>> triangles;
    fscanf(fid, "%s", shape_type);
    if (shape_type[0] == 'P'){
        vector<float> a = fscanf_n(fid,6);
        int num_pts = (int)a[5];
        for (int i = 0;i<num_pts;i++){
            pt_coordinates.push_back(fscanf_n(fid, 3));
        }
        for (int i = 0;i<num_pts;i++){
            pt_normals.push_back(fscanf_n(fid,3));
        }
        // get num of triangles;
        int num_triangles;
        float tmp;
        fscanf(fid, "%f", &tmp);
        num_triangles = (int)tmp;
    
        // get color type
        fscanf(fid, "%f", &tmp);
        int color_type = (int)tmp;
        vector<vector<float> > pt_colors;
        if (color_type == 2){
            for (int i = 0;i<num_pts;i++){
                pt_colors.push_back(fscanf_n(fid, 4));
            }
        }
        else if (color_type == 0){
            pt_colors.push_back(fscanf_n(fid, 4));
        }
        fscanf_n(fid,num_triangles);
        // load triangles
        for (int i = 0;i<num_triangles;i++){
            vector<float> tmp = fscanf_n(fid,3);
            // adjust for matlab index starting from 1 instead of 0
            tmp[0]+=1;
            tmp[1]+=1;
            tmp[2]+=1;
            triangles.push_back(tmp);
        }  
    }
    if (shape_type[0] == 'L'){
        vector<float> a = fscanf_n(fid, 2);
        int num_pts = (int)a[1];
        for (int i = 0;i<num_pts;i++){
            pt_coordinates.push_back(fscanf_n(fid,3));
        }
    }
    // cout<< "coordinates: "<<pt_coordinates.size()<<' '<<pt_coordinates[0].size()<<endl;
    // cout<< "normals "<<pt_normals.size()<<' '<<pt_normals[0].size()<<endl;
    // cout<< "triangles "<<triangles.size()<<' '<<triangles[0].size()<<endl;
    fclose(fid);
    // vector<vector<vector<float>>> res{pt_coordinates, pt_normals, triangles};
    // res.push_back(pt_coordinates);
    // res->push_back(pt_normals);
    // res->push_back(triangles);
 
    
}

vector<float> fscanf_n(FILE* fid, int n){ // scanf n times
    vector<float> data(n);
    for (int i = 0;i<n;i++){
        float tmp;
        fscanf(fid, "%f", &tmp);
        data[i] = tmp;
    }
    return data;

}

// int main(){
//     string Mesh = "/Users/samr/Documents/GitHub/trekker/data/lh_white_DTISpace.obj";
//     vector<vector<float> > coord, nS, trg;
//  	ReadObjShape(Mesh, coord, nS, trg);
//     return 1;
// }