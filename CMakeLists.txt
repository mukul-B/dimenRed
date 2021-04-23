cmake_minimum_required(VERSION 3.17)
project(dimenRed)

set(CMAKE_CXX_STANDARD 14)

add_executable(dimenRed image.h image.cpp  ReadImage.cpp ReadImageHeader.cpp WriteImage.cpp Threshold.cpp jacobi.h jacobi.cpp lamdafunctions.h matrixfunctions.cpp matrixOperations.cpp matrixOperations.h matrixfunctions.h)