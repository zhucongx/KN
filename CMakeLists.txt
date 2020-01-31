cmake_minimum_required(VERSION 3.15)

# set the project name
project(KN)

set(CMAKE_CXX_STANDARD 17)

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin)
set(LIBRARY_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/lib)

set(KN_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/include)
set(KN_SOURCE_DIR ${PROJECT_SOURCE_DIR}/src)

include_directories(${KN_INCLUDE_DIR})
include_directories(${KN_SOURCE_DIR})

add_subdirectory(external/keras2cpp)

file(GLOB KN_SRC
        ${KN_SOURCE_DIR}/*.cpp
        ${KN_INCLUDE_DIR}/*)
add_executable(${PROJECT_NAME} ${KN_SRC})
#add_executable(${PROJECT_NAME} src/Atom.cpp src/Configuration.cpp)

include(FindPkgConfig)

find_package(armadillo REQUIRED)
find_package(keras2cpp REQUIRED)
find_package(mpi REQUIRED)
include_directories(${LIBARMADILLO_INCLUDE_DIR})
include_directories(${LIBKERAS2CPP_INCLUDE_DIR})
target_link_libraries(${PROJECT_NAME} ${LIBARMADILLO_LIBRARY})
target_link_libraries(${PROJECT_NAME} ${LIBKERAS2CPP_LIBRARY})
