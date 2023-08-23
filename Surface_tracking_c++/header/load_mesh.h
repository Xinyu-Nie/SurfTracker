#include <iostream>
#include <string>
#include <fstream>
#include <stdio.h>
#include <vector>
using namespace std;
void ReadObjShape(string fname, vector<vector<float> >& pt_coordinates, vector<vector<float> >& pt_normals, vector<vector<float> >& triangles);
vector<float> fscanf_n(FILE* fid, int n);