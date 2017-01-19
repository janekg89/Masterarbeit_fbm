cmake_minimum_required(VERSION 3.5)
project(lown_cython)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/")
find_package(FFTW3)
find_package(Boost COMPONENTS system REQUIRED)

add_executable(fractest main.cpp _generatefracincrements.cpp _generatefracincrements.h)
target_link_libraries(fractest ${Boost_LIBRARIES} ${FFTW3_LIBRARIES} gsl gslcblas)
target_include_directories(fractest PUBLIC ${Boost_INCLUDE_DIRS} ${FFTW_INCLUDE_DIRS})
