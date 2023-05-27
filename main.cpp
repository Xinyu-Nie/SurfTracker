#ifdef BUILD_FOR_WINDOWS
#include<windows.h>
#include<io.h>
#else 
#include <unistd.h>
#include <sys/time.h>
#endif

#include <iostream>
#include "SurfacePro/SurfaceFODPro.h"
#include "header/Ufibertracking.h"
#include <unordered_map>
#include <sys/types.h>
#include <sys/stat.h>
//#include "header/load_mat.h"

using namespace std;
unordered_map<string, string> arg_map;
/*
* check file path is valid, if dir is not existed, then create one
*/
#ifdef BUILD_FOR_WINDOWS
bool createDirectoryIfNotExist(const std::string& path) {
    DWORD ftyp = GetFileAttributesA(path.c_str());
    if (ftyp == INVALID_FILE_ATTRIBUTES) {
        // �½�Ŀ¼
        if (CreateDirectoryA(path.c_str(), NULL)) {
            std::cout << "Directory created: " << path << std::endl;
            return true;
        }
        else {
            std::cerr << "Failed to create directory: " << path << std::endl;
            return false;
        }
    }

    if (ftyp & FILE_ATTRIBUTE_DIRECTORY) {
        // Ŀ¼�Ѵ���
        return true;
    }

    std::cerr << "Path exists, but not a directory: " << path << std::endl;
    return false;
}
#else
void createDirectoryIfNotExist(const std::string& path) {
    if (!path.size()) return;
    if (path.back() == '/' || path.back() == '\\') {
        if (access(path.c_str(), F_OK) == -1) {
            if (mkdir(path.c_str(), S_IRWXU) != 0)
            {
                cout << "mkdir failed" << endl;
                exit(0);
            }
        }
    }
    else {
        for (int i = path.size() - 1; i >= 0; i--) {
            if (path[i] == '/' || path[i] == '\\') {
                string dir(path.begin(), path.begin() + i + 1);
                if (access(dir.c_str(), F_OK) == -1) {
                    if (mkdir(dir.c_str(), S_IRWXU) != 0)
                    {
                        cout << "mkdir failed" << endl;
                        exit(0);
                    }
                }
            }
        }
    }
    return;
}
#endif
/*
* parse the commandline
*/
void parser(int argc, char* argv[]){
    vector<string> argS(argc);
    
    for (int i=1;i<argc;i++){
        argS[i] = string(argv[i]);
    }
    for (int i = 1;i<argc;i+=2){
        if (argS[i]=="-preprocess"){
            arg_map["preprocess"] = argS[i+1];
        }
        else if (argS[i]=="-track"){
            arg_map["track"] = argS[i+1];
        }
        else if (argS[i]=="-shrinkdis"){
            arg_map["shrinkdis"] = argS[i+1];
        }
        else if (argS[i]=="-minFOD"){
            arg_map["minFOD"] = argS[i+1];
        }
        else if (argS[i]=="-trackdegree"){
            arg_map["trackdegree"] = argS[i+1];
        }
        else if (argS[i]=="-iteration"){
            arg_map["iteration"] = argS[i+1];
        }
        else if (argS[i]=="-gyrusflag"){
            arg_map["gyrusflag"] = argS[i+1];
        }
        else if (argS[i]=="-fod"){
            arg_map["fod"] = argS[i+1];
        }
        else if (argS[i]=="-mesh"){
            arg_map["mesh"] = argS[i+1];
        }
        else if (argS[i]=="-roi"){
            arg_map["roi"] = argS[i+1];
        }
        else if (argS[i]=="-sulcus"){
            arg_map["sulcus"] = argS[i+1];
        }
        else if (argS[i]=="-gyrus1"){
            arg_map["gyrus1"] = argS[i+1];
        }
        else if (argS[i]=="-gyrus2"){
            arg_map["gyrus2"] = argS[i+1];
        }
        else if (argS[i]=="-output_trk"){
            arg_map["output"] = argS[i+1];
            //createDirectoryIfNotExist(arg_map["output"]);
        }
        else if (argS[i]=="-fodsphere"){
            arg_map["fodsphere"] = argS[i+1];
            //createDirectoryIfNotExist(arg_map["fodsphere"]);
        }
        else {
            cout<<argS[i]<<endl;
            cout<<"wrong input"<<endl;
            exit(0);
        }
    }
}
int main(int argc, char * argv[]) // 
{
    parser(argc, argv);
    string FOD, Mesh, ROI, SulcusLabel, GyrusLabel1, GyrusLabel2, output, saved_data_path;
    float minFOD = stof(arg_map["minFOD"]);
    float trackdegree = stof(arg_map["trackdegree"]);
    int iteration = stoi(arg_map["iteration"]);
    int GyrusFlag = stoi(arg_map["gyrusflag"]);
    float shrinkdis = stof(arg_map["shrinkdis"]);
    FOD = arg_map["fod"];
    Mesh = arg_map["mesh"];
    ROI = arg_map["roi"];
    SulcusLabel = arg_map["sulcus"];
    GyrusLabel1 = arg_map["gyrus1"];
    GyrusLabel2 = arg_map["gyrus2"];
    saved_data_path = arg_map["fodsphere"];
    
    if (arg_map["preprocess"] == "1") {
        cout<<"===============Start SurfacePreprocess================"<<endl;
        SurfaceFODPro(FOD, Mesh, ROI, saved_data_path, shrinkdis);
        cout<<"===============Complete SurfacePreprocess================"<<endl;
    }

    vector<vector<float> > FODProjected, trgS, coordS, TRIcen, trinerbor;
    vector<float> THETA;
    vector<int> dim(4);
    vector<float> res(3);
 
    if (arg_map["track"]=="1"){
        output = arg_map["output"];
        FILE* pfile = fopen(saved_data_path.c_str(), "rb");
        binary_read(pfile, trinerbor);
        binary_read(pfile, trgS);
        binary_read(pfile, coordS);
        binary_read(pfile, TRIcen);
        binary_read(pfile, THETA);
        binary_read(pfile, FODProjected);
        fread(&dim[0], sizeof(int), 4, pfile);
        fread(&res[0], sizeof(float), 3, pfile);
        fclose(pfile);
        cout<<"===============Start Tracking================"<<endl;
        Ufibertracking(output, Mesh, ROI,SulcusLabel,GyrusLabel1,GyrusLabel2,iteration, coordS,trgS,TRIcen, THETA,FODProjected, trinerbor, dim,res,minFOD,trackdegree, GyrusFlag);
        cout<<"===============Complete Traking================"<<endl;


    }
	return 1;
}

