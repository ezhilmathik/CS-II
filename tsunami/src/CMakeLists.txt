cmake_minimum_required(VERSION 3.9)

project(TsunAWI Fortran)

set(src_home ${CMAKE_CURRENT_LIST_DIR}) 
file(GLOB sources_all ${src_home}/*.F90)

add_executable(${PROJECT_NAME} ${sources_all})
set_target_properties(${PROJECT_NAME} PROPERTIES LINKER_LANGUAGE Fortran)

# Make it possible to switch off netcdf and/or foxy and/or double precision arithmetic
option(NETCDF "netcdf support" ON)
option(FOXY "XML support with libfoxy" OFF)

set(WITH_FLOATING_PRECISION "real64" CACHE STRING "floating point arithmetic, e.g. real64 (default: real32)")

if(NETCDF)
    include(${CMAKE_CURRENT_LIST_DIR}/../cmake/FindNETCDF.cmake)
    target_include_directories(${PROJECT_NAME} PRIVATE ${NETCDF_Fortran_INCLUDE_DIRECTORIES})
    target_link_libraries(${PROJECT_NAME} ${NETCDF_Fortran_LIBRARIES} ${NETCDF_C_LIBRARIES})	
else(NETCDF)
    add_definitions(-DNO_NETCDF)
endif(NETCDF)

