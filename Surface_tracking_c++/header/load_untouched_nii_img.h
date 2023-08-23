#include "nifti1_io.h"
#include <iostream>
#include <string>
#include <vector>
#include <stdio.h>
#include <iostream>
#include <fstream>
using namespace std;

class NIImage{
    public:
    nifti_image* nim;
    string path;
    vector<vector<vector<vector<float>>>> img; //t z y x
    

    NIImage(string path);
    ~NIImage();
    void load_untouched_nii();
    nifti_image* get_nim();
    void printInfo();
    vector<vector<vector<vector<float>>>> get_img();
    void saveInfo(string path);

};
// void printInfo(nifti_image* nim);
// vector<vector<vector<vector<float>>>> load_untouched_nii(string FOD_path);

// namespace NII{
//     nifti_image* nim;

// };