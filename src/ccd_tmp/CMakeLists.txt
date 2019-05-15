cmake_minimum_required(VERSION 3.1)

add_library(ccd2d
		src/Curve.cpp
		src/Hash.cpp
		src/InterferenceVolume.cpp
)


target_include_directories(ccd2d PRIVATE ${EIGEN_INCLUDE_DIR})
target_include_directories(ccd2d PRIVATE src)
target_include_directories(ccd2d PUBLIC include)
