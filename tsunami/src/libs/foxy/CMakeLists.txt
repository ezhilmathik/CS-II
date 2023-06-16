cmake_minimum_required(VERSION 3.9)

project(foxy Fortran)

set(src_home ${CMAKE_CURRENT_LIST_DIR}) 
file(GLOB all_sources ${src_home}/*.F90 ${src_home}/*.f90)


add_library(${PROJECT_NAME} STATIC ${all_sources})


target_include_directories(${PROJECT_NAME} PUBLIC ${PROJECT_BINARY_DIR})


