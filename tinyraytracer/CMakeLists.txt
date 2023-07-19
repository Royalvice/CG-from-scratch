cmake_minimum_required(VERSION 3.15)
project(tinyraytracer)

set(CMAKE_CXX_STANDARD 14)
set(EGE_INC_DIR F:\\cpp_lib\\EGE20.08\\include)
set(EGE_LINK_DIR F:\\cpp_lib\\EGE20.08\\lib\\mingw64)

include_directories(${EGE_INC_DIR})
link_directories(${EGE_LINK_DIR})

add_executable(tinyraytracer raytrace.cpp raytrace.h main.cpp scene.h scene.h)

target_link_libraries(tinyraytracer libgraphics64.a
        libgdi32.a libimm32.a libmsimg32.a libole32.a
        liboleaut32.a libwinmm.a libuuid.a libgdiplus.a)
