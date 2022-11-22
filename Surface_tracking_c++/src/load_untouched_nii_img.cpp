#include "../header/load_untouched_nii_img.h"
// using namespace NII;

NIImage::NIImage(string path){
	this->path = path;
	load_untouched_nii();
	// this->img = vector<vector<vector<vector<float>>>>(nx, vector<vector<vector<float>>>(ny, vector<vector<float>>(nz, vector<float>(nt,0))));
}
NIImage::~NIImage(){
	free(this->nim);
	this->img.clear();
}
void NIImage::printInfo() {

	std::cout << std::endl;
	std::cout << "Data description:       " << this->nim->descrip << std::endl;
	std::cout << "Nifti type:             ";
	switch (this->nim->nifti_type) {
	case 0:    std::cout<<"Analyze"                          << std::endl; break;
	case 1:    std::cout<<"NIFTI-1 (1 file)"                 << std::endl; break;
	case 2:    std::cout<<"NIFTI-2 (2 file)"                 << std::endl; break;
	case 3:    std::cout<<"NIFTI-ASCII (1 file)"             << std::endl; break;
	default:   std::cout<<"Not applicable"                   << std::endl; break;
	}

	std::cout << "Number of dimensions:   "   << this->nim->dim[0] << std::endl;
	std::cout << "Dimensions:             [";
	for (int i=1; i<(this->nim->dim[0]+1); i++) {
		std::cout << this->nim->dim[i];
		if (i!=this->nim->dim[0])
			std::cout << " x ";
	}
	std::cout << "]" << std::endl;

	std::cout << "Number of bytes/voxel:  "   << this->nim->nbyper << std::endl;

	std::cout << "Pixdim:                 [" << this->nim->pixdim[1] << " x " << this->nim->pixdim[2] << " x " << this->nim->pixdim[3] << "]" << std::endl;

	std::cout << "Unit of pixdim:         ";
	switch (this->nim->xyz_units) {
	case 0:    std::cout<<"unknown"                          << std::endl; break;
	case 1:    std::cout<<"meter"                            << std::endl; break;
	case 2:    std::cout<<"mm"                               << std::endl; break;
	case 3:    std::cout<<"micron"                           << std::endl; break;
	default:   std::cout<<"Not applicable"                   << std::endl; break;
	}

	std::cout << "Datatype:               " ;
	switch (this->nim->datatype) {
	case 2:    std::cout<<"uint8, unsigned char"             << std::endl; break;
	case 4:    std::cout<<"int16, signed short"              << std::endl; break;
	case 8:    std::cout<<"int32, signed int"                << std::endl; break;
	case 16:   std::cout<<"float32"                          << std::endl; break;
	case 64:   std::cout<<"float64, float"                  << std::endl; break;
	case 256:  std::cout<<"int8, signed char"                << std::endl; break;
	case 512:  std::cout<<"uint16, unsigned short"           << std::endl; break;
	case 768:  std::cout<<"uint32, unsigned int"             << std::endl; break;
	case 1024: std::cout<<"int64, long long int"             << std::endl; break;
	case 1280: std::cout<<"uint64, unsigned long long int"   << std::endl; break;
	case 1536: std::cout<<"float128, long float"            << std::endl; break;
	}

}

void NIImage::load_untouched_nii(){
    
	this->nim = nifti_image_read(this->path.c_str(),1);
	// img_FOD->printInfo();
	// printInfo();
	// load untouched nii
	// read the dataset header
	void* data = (float*)this->nim->data;
	int nx = this->nim->dim[1], ny = this->nim->dim[2], nz = this->nim->dim[3], nt = this->nim->dim[4];
	int img_size = nx*ny*nz*nt;
	this->img = vector<vector<vector<vector<float>>>>(nt, vector<vector<vector<float>>>(nz, vector<vector<float>>(ny, vector<float>(nx,0))));
	
	for (int t = 0;t<nt;t++){
		for (int z = 0;z<nz;z++){
			for (int y = 0;y<ny;y++){
				for (int x = 0;x<nx;x++){
					this->img[t][z][y][x] = *(float*)data;
					// if (this->img[x][y][z][t]!=0) cout<< this->img[x][y][z][t]<< ' ';
					data = (char*)data + sizeof(float);
				}
			}
		}
	}
	
}

vector<vector<vector<vector<float>>>> NIImage::get_img(){
	return this->img;
}

nifti_image* NIImage::get_nim(){
	return this->nim;
}

void NIImage::saveInfo(string path){
	FILE* fp = fopen(path.c_str(),"w");
	for (int i=1;i<=4;i++){
		fprintf(fp,"%d ",this->nim->dim[i]);
	}
	for (int i=1;i<=3;i++){
		fprintf(fp,"%f ",this->nim->pixdim[i]);
	}
	fclose(fp);
}

// int main(){
	
// 	string FOD = "C:\\Users\\SAM\\trekker\\data\\100206_FOD_Order16_WholeVol.nii";
	
// 	NIImage nii(FOD);
// 	auto img = nii.img;
// 	// printf("%.10f",img[59][53][110][0]);
// 	cout<<img[0][110][53][59]<<endl;
// 	return 1;
// }