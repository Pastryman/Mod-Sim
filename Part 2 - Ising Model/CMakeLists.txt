cmake_minimum_required(VERSION 3.9)
project(Part_2___Ising_Model)

set(CMAKE_CXX_STANDARD 11)

add_executable(Exc1.1 IsingModel.cpp)
add_executable(Exc1.2 IsingModel_Corr.cpp)
add_executable(Exc2.1 IsingModel_Wolff.cpp)