#include<fstream>
#include<vector>
using namespace std;
inline void read_binary(string path, vector<float>& data){
    ifstream inFile(path,ios::in|ios::binary);
	int s;
	while(inFile.read((char *)&s, sizeof(int))) { //一直读到文件结束
        data.push_back(s);
    }
	inFile.close();
}

inline void read_binary(string path, vector<vector<float> >& data){
    // data.push_back({});
    ifstream inFile(path,ios::in|ios::binary);
	int s;
	while(inFile.read((char *)&s, sizeof(int))) { //一直读到文件结束
        data.push_back({(float)s});
    }
	inFile.close();
}