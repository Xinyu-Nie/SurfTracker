main: 
	g++ main.cpp -pthread -std=c++17 -o ./main ./SurfacePro/*.cpp jacobiP/*.cpp ./src/*.cpp ./src/*.c ./lsqr/*.cpp -w -D HAVE_ZLIB -lz -g -I"C:\Program Files\MATLAB\R2022a\extern/include" "C:\Program Files\MATLAB\R2022a\extern\lib\win64\mingw64" -lmat -lmx -Wl,-rpath "C:\Program Files\MATLAB\R2022a\extern\lib\win64\mingw64" -D BUILD_FOR_WINDOWS

test:
	g++ test.cpp -pthread -std=c++17 -o test ./SurfacePro/*.cpp jacobiP/*.cpp ./src/*.cpp ./src/*.c ./lsqr/*.cpp -w -D HAVE_ZLIB -lz -g