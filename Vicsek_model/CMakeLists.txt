cmake_minimum_required(VERSION 3.9)
project(Vicsek_model)

set(CMAKE_CXX_STANDARD 17)

add_executable(Vicsek_model main.cpp)
add_executable(Vicsek_model_box main_boxes.cpp)
add_executable(Potential main_potential.cpp)