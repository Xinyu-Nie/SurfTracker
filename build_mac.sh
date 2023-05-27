#!/bin/bash

cmakeExe="/usr/local/bin/cmake"
buildType="Release"

rm -rf build/macOS
mkdir -p build/macOS
"$cmakeExe" . -B ./build/macOS -DCMAKE_BUILD_TYPE="$buildType"
cd ./build/macOS
make
cp ./main ../../bin/main_macos_x64