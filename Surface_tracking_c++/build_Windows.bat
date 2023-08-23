@SET cmakeExe="C:/Program Files/cmake-3.23.1-windows-x86_64/bin/cmake"
@SET buildType=Release

@rmdir build\Windows /s/q
@mkdir build\Windows
@%cmakeExe% . -B ./build/Windows -DCMAKE_BUILD_TYPE=Release
@cd ./build/Windows
@msbuild SurfTracker.sln /property:Configuration=Release
@copy .\Release\main.exe ..\..\