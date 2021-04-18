cmake_minimum_required(VERSION 3.17)
project(dimenRed)

set(CMAKE_CXX_STANDARD 14)

add_executable(dimenRed main.cpp worked.cpp image.h image.cpp  ReadImage.cpp ReadImageHeader.cpp WriteImage.cpp Threshold.cpp  MatricOP.h DiscriminantCases.h)