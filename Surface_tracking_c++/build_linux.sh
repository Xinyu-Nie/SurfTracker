#!/bin/bash

cmakeExe="/usr/bin/cmake"
buildType="Release"

rm -rf build/Linux
mkdir build/Linux
${cmakeExe} . -B ./build/Linux -DCMAKE_BUILD_TYPE=${buildType}
cd build/Linux
make
cp ./main ../../main
