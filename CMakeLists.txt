cmake_minimum_required(VERSION 3.10)

#SET(CMAKE_C_COMPILER icc)
#SET(CMAKE_CXX_COMPILER icpc)

project(KN)
set(CMAKE_CXX_STANDARD 17)

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17"
        CACHE INTERNAL "")

if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    message("GUN flags")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}"
            CACHE INTERNAL ""
            )
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
    message("AppleClang flags")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}"
            CACHE INTERNAL ""
            )
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
    message("Intel flags")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}"
            CACHE INTERNAL ""
            )
endif ()



set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g3 -ggdb"
        CACHE INTERNAL "")

set(CMAKE_CXX_FLAGS_RELEASE "-O3"
        CACHE INTERNAL "")

file(GLOB KN_SRC
        ${PROJECT_SOURCE_DIR}/include/*.h
        ${PROJECT_SOURCE_DIR}/src/*.cpp)

add_executable(${PROJECT_NAME} ${KN_SRC})
set_target_properties(${PROJECT_NAME} PROPERTIES OUTPUT_NAME "kn.exe")

target_include_directories(${PROJECT_NAME} PRIVATE
        ${PROJECT_SOURCE_DIR}/include)

#include(FindPkgConfig)

find_package(MPI REQUIRED)
# include directories
target_include_directories(${PROJECT_NAME} INTERFACE ${MPI_C_INCLUDE_PATH} ${MPI_CXX_INCLUDE_PATH})
# link libraries
target_link_libraries(${PROJECT_NAME} INTERFACE ${MPI_C_LIBRARIES} ${MPI_CXX_LIBRARIES})
# set interface target compile options
target_compile_options(${PROJECT_NAME} INTERFACE
        "${MPI_C_COMPILE_DEFINITIONS} ${MPI_CXX_COMPILE_DEFINITIONS}")
target_compile_options(${PROJECT_NAME} INTERFACE
        "${MPI_C_COMPILE_OPTIONS} ${MPI_CXX_COMPILE_OPTIONS}")
set_target_properties(${PROJECT_NAME} PROPERTIES
        INTERFACE_LINK_OPTIONS "${MPI_C_LINK_FLAGS} ${MPI_CXX_LINK_FLAGS}")


#include_directories(${MPI_INCLUDE_PATH})
#target_link_libraries(${PROJECT_NAME} ${MPI_LIBRARIES})
#set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${MPI_C_COMPILE_FLAGS} ${MPI_C_LINK_FLAGS}")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${MPI_CXX_COMPILER_FLAGS} ${MPI_CXX_LINK_FLAGS}")

find_package(OpenMP REQUIRED)
# include directories
include_directories(${PROJECT_NAME} INTERFACE ${OpenMP_C_INCLUDE_DIRS} ${OpenMP_CXX_INCLUDE_DIRS})
# link libraries
target_link_libraries(${PROJECT_NAME} INTERFACE ${OpenMP_C_LIBRARIES} ${OpenMP_CXX_LIBRARIES})
# set interface target compile options
target_compile_options(${PROJECT_NAME} INTERFACE
        "${OpenMP_C_FLAGS} ${OpenMP_CXX_FLAGS}")

#include_directories(${OpenMP_CXX_INCLUDE_DIRS})
#include_directories(${OpenMP_C_INCLUDE_DIRS})
#target_link_libraries(${PROJECT_NAME} ${OpenMP_CXX_LIBRARIES})
#target_link_libraries(${PROJECT_NAME} ${OpenMP_C_LIBRARIES})
#set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")


find_package(ARMADILLO REQUIRED)
include_directories(${PROJECT_NAME} INTERFACE ${ARMADILLO_INCLUDE_DIRS})
target_link_libraries(${PROJECT_NAME} INTERFACE ${ARMADILLO_LIBRARIES})

## Lib keras2cpp
#add_subdirectory(external/keras2cpp)
#include_directories(external/keras2cpp)
#target_link_libraries(${PROJECT_NAME} Keras2cpp)


