#include<vector>
#include<string>
using namespace std;
// Masktrimesh(trg, coord, MTL, coordS, trgS);
inline void Masktrimesh(vector<vector<float> >& faces0, vector<vector<float> >& vertices0, vector<float>& mask, vector<vector<float> >& vertices, vector<vector<float> >& faces){
	vector<int> Position;
	for (int i = 0;i<mask.size();i++){
		if (mask[i]>0) Position.push_back(i);
	}

	int Lv = Position.size();
	int L2 = vertices0[0].size();
	
	vertices = vector<vector<float> >(Lv, vector<float>(L2,0));

	for (int h = 0;h<Lv;h++){
		auto vert = vertices0[Position[h]];
		vertices[h] = vert;
	}
	int Lf1 = faces0.size(), Lf2 = faces0[0].size();
	// faces = vector<vector<float>> (Lf1, vector<float>(Lf2,0));
	// int t =0;
	
	for (int k = 0;k<Lf1;k++){
		auto v = faces0[k];
		v[0] = v[0]-1;
		v[1] = v[1]-1;
		v[2] = v[2]-1;
		if ((mask[(int)v[0]] + mask[(int)v[1]] + mask[(int)v[2]]) > 2.5) {
			float a1,a2,a3;
			for (int i = 0;i<Position.size();i++){
				if (Position[i] == (int)v[0]){
					a1 = i;
					break;
				}
			}
			for (int i = 0;i<Position.size();i++){
				if (Position[i] == (int)v[1]){
					a2 = i;
					break;
				}
			}
			for (int i = 0;i<Position.size();i++){
				if (Position[i] == (int)v[2]){
					a3 = i;
					break;
				}
			}
			faces.push_back({a1,a2,a3});
			
		}
	}
	// faces = 
	// return {vertices, vector<vector<float>>(faces.begin(), faces.begin()+t)};

}

inline void Masktrimesh(vector<vector<float> >& faces0, vector<float>& vertices0, vector<float>& mask, vector<float>& vertices, vector<vector<float> >& faces){
	vector<int> Position;
	for (int i = 0;i<mask.size();i++){
		if (mask[i]>0) Position.push_back(i);
	}

	int Lv = Position.size();
	// int L2 = vertices0[0].size();
	
	vertices = vector<float>(Lv,0);

	for (int h = 0;h<Lv;h++){
		auto vert = vertices0[Position[h]];
		vertices[h] = vert;
	}
	int Lf1 = faces0.size(), Lf2 = faces0[0].size();
	// faces = vector<vector<float>> (Lf1, vector<float>(Lf2,0));
	// int t =0;
	
	for (int k = 0;k<Lf1;k++){
		auto v = faces0[k];
		v[0] = v[0]-1;
		v[1] = v[1]-1;
		v[2] = v[2]-1;
		if ((mask[(int)v[0]] + mask[(int)v[1]] + mask[(int)v[2]]) > 2.5) {
			float a1,a2,a3;
			for (int i = 0;i<Position.size();i++){
				if (Position[i] == (int)v[0]){
					a1 = i;
					break;
				}
			}
			for (int i = 0;i<Position.size();i++){
				if (Position[i] == (int)v[1]){
					a2 = i;
					break;
				}
			}
			for (int i = 0;i<Position.size();i++){
				if (Position[i] == (int)v[2]){
					a3 = i;
					break;
				}
			}
			faces.push_back({a1,a2,a3});
			
		}
	}
	// faces = 
	// return {vertices, vector<vector<float>>(faces.begin(), faces.begin()+t)};

}

inline void Masktrimesh2(vector<vector<float> >& faces0, vector<float>& vertices0, vector<float>& mask, vector<float>& vertices, vector<vector<float> >& faces, vector<float>& faceOri){
	vector<int> Position;
	for (int i = 0;i<mask.size();i++){
		if (mask[i]>0) Position.push_back(i);
	}

	int Lv = Position.size();
	// int L2 = vertices0[0].size();
	
	vertices = vector<float>(Lv,0);

	for (int h = 0;h<Lv;h++){
		auto vert = vertices0[Position[h]];
		vertices[h] = vert;
	}
	int Lf1 = faces0.size(), Lf2 = faces0[0].size();
	// faces = vector<vector<float>> (Lf1, vector<float>(Lf2,0));
	// int t =0;
	
	for (int k = 0;k<Lf1;k++){
		auto v = faces0[k];
		v[0] = v[0]-1;
		v[1] = v[1]-1;
		v[2] = v[2]-1;
		if ((mask[(int)v[0]] + mask[(int)v[1]] + mask[(int)v[2]]) > 2.5) {
			float a1,a2,a3;
			for (int i = 0;i<Position.size();i++){
				if (Position[i] == (int)v[0]){
					a1 = i;
					break;
				}
			}
			for (int i = 0;i<Position.size();i++){
				if (Position[i] == (int)v[1]){
					a2 = i;
					break;
				}
			}
			for (int i = 0;i<Position.size();i++){
				if (Position[i] == (int)v[2]){
					a3 = i;
					break;
				}
			}
			faces.push_back({a1,a2,a3});
			faceOri.push_back(k);
		}
	}

}