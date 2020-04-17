cmake_minimum_required(VERSION 3.9)
set(CMAKE_C_COMPILER "mpicc")
set(CMAKE_CXX_COMPILER "mpicxx")

project(KN)
set(CMAKE_CXX_STANDARD 20)

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")

set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g3 -ggdb")

set(CMAKE_CXX_FLAGS_RELEASE "-O3 -ffast-math -fno-rtti")


file(GLOB KN_SRC
        ${PROJECT_SOURCE_DIR}/src/*.cpp
        src/include/*.h)

add_executable(${PROJECT_NAME} ${KN_SRC})
set_target_properties(${PROJECT_NAME} PROPERTIES OUTPUT_NAME "kn.exe")

#include(FindPkgConfig)
find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})
target_link_libraries(${PROJECT_NAME} ${MPI_LIBRARIES})
set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${MPI_COMPILE_FLAGS} ${MPI_LINK_FLAGS}")
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${MPI_COMPILE_FLAGS} ${MPI_LINK_FLAGS}")


find_package(OpenMP REQUIRED)
include_directories(${OpenMP_CXX_INCLUDE_DIRS})
include_directories(${OpenMP_C_INCLUDE_DIRS})
target_link_libraries(${PROJECT_NAME}  ${OpenMP_CXX_LIBRARIES})
target_link_libraries(${PROJECT_NAME}  ${OpenMP_C_LIBRARIES})
set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")


find_package(ARMADILLO REQUIRED)
include_directories(${LIBARMADILLO_INCLUDE_DIR})
target_link_libraries(${PROJECT_NAME} ${LIBARMADILLO_LIBRARY})

## Lib keras2cpp
#add_subdirectory(external/keras2cpp)
#include_directories(external/keras2cpp)
#target_link_libraries(${PROJECT_NAME} Keras2cpp)


