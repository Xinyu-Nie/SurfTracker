# set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/..)
set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g -w")
if (UNIX AND NOT APPLE)
	message(STATUS "Building for linux")
	cmake_minimum_required (VERSION 2.8)
	set(CMAKE_CXX_STANDARD 17)
	project (SurfTracker)
	aux_source_directory(./src SRC1)
	aux_source_directory(./SurfacePro SRC2)
	aux_source_directory(./jacobiP SRC3)
	aux_source_directory(./lsqr SRC4)
	add_executable(main main.cpp ${SRC1} ${SRC2} ${SRC3} ${SRC4})
	add_compile_definitions(HAVE_ZLIB)
	include_directories(./lsqr)
	include_directories(./jacobiP)
	include_directories(./header)
	target_link_libraries(main pthread z)
elseif (WIN32)
	message(STATUS "Building for Windows")
    	project (SurfTracker)
	aux_source_directory(./src SRC1)
	aux_source_directory(./SurfacePro SRC2)
	aux_source_directory(./jacobiP SRC3)
	aux_source_directory(./lsqr SRC4)
	add_compile_definitions(BUILD_FOR_WINDOWS)
	add_executable(main main.cpp ${SRC1} ${SRC2} ${SRC3} ${SRC4})
	include_directories(./boost_1_66_0)
	include_directories(./lsqr)
	include_directories(./jacobiP)
	include_directories(./header)
elseif (APPLE)
    message(STATUS "Building for macOS")
    cmake_minimum_required (VERSION 2.8)
    set(CMAKE_CXX_STANDARD 17)
    project (SurfTracker)
    aux_source_directory(./src SRC1)
    aux_source_directory(./SurfacePro SRC2)
    aux_source_directory(./jacobiP SRC3)
    aux_source_directory(./lsqr SRC4)
    add_executable(main main.cpp ${SRC1} ${SRC2} ${SRC3} ${SRC4})
    add_compile_definitions(BUILD_FOR_MAC)
    add_compile_definitions(HAVE_ZLIB)
    include_directories(./boost_1_66_0)
    include_directories(./lsqr)
    include_directories(./jacobiP)
    include_directories(./header)
    target_link_libraries(main pthread z)
endif()