main: 
	g++ main.cpp -pthread -std=c++17 -o ./main ./SurfacePro/*.cpp jacobiP/*.cpp ./src/*.cpp ./src/*.c ./lsqr/*.cpp -w -D HAVE_ZLIB -lz -g -I/usr/local/MATLAB/R2019b/extern/include -L/usr/local/MATLAB/R2019b/bin/glnxa64 -lmat -lmx -Wl,-rpath /usr/local/MATLAB/R2019b/bin/glnxa64

test:
	g++ test.cpp -pthread -std=c++17 -o test ./SurfacePro/*.cpp jacobiP/*.cpp ./src/*.cpp ./src/*.c ./lsqr/*.cpp -w -D HAVE_ZLIB -lz -g