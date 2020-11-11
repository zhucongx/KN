cmake_minimum_required(VERSION 3.10)

#set(CMAKE_C_COMPILER "mpicc") # {gcc, mpicc, icc}
#set(CMAKE_CXX_COMPILER "mpicxx") # {g++, mpicxx, icc}

project(KN)
set(CMAKE_CXX_STANDARD 20)

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin)
set(LIBRARY_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/lib)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} "
        CACHE INTERNAL "")
#set(CMAKE_VERBOSE_MAKEFILE ON)
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} \
-Wall -Wextra -Wshadow -Wnon-virtual-dtor -Wpointer-arith -Wold-style-cast -Wpedantic  \
-Wcast-qual -Wunused -Woverloaded-virtual -Winline -Wunreachable-code -Wcast-align  -Weffc++ \
-Wundef -Wmisleading-indentation -Wredundant-decls -pedantic-errors \
-Wnull-dereference  -Wdouble-promotion -Wconversion -Wfatal-errors  -Wextra-semi \
-O0 -g3 -finline-functions"
        CACHE INTERNAL "")
# -Wconversion
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -Ofast -Wall -Wextra"
        CACHE INTERNAL "")

message("-- Detecting ${CMAKE_CXX_COMPILER_ID} compiler")
if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE}"
            CACHE INTERNAL "")
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE}"
            CACHE INTERNAL "")
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE}"
            CACHE INTERNAL "")
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE}"
            CACHE INTERNAL "")
else ()
    message(WARNING "Please check which compiler you are using")
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -D_GLIBCXX_PARALLEL"
            CACHE INTERNAL "")
endif ()
##############################################################################################
#add_executable(${PROJECT_NAME} Kn_Main.cpp)
#add_subdirectory(kn/cfg)
#add_subdirectory(kn/gen)
#add_subdirectory(kn/ansys)
#add_subdirectory(kn/kmc)
#target_link_libraries(${PROJECT_NAME} PUBLIC cfg gen ansys kmc)
##############################################################################################
file(GLOB PROJECT_SRC
        ${PROJECT_SOURCE_DIR}/kn/cfg/src/*.cpp
        ${PROJECT_SOURCE_DIR}/kn/gen/src/*.cpp
        ${PROJECT_SOURCE_DIR}/kn/ansys/src/*.cpp
        ${PROJECT_SOURCE_DIR}/kn/kmc/src/*.cpp
        )
set(PROJECT_INCLUDE
        ${PROJECT_SOURCE_DIR}/kn/cfg/include
        ${PROJECT_SOURCE_DIR}/kn/gen/include
        ${PROJECT_SOURCE_DIR}/kn/ansys/include
        ${PROJECT_SOURCE_DIR}/kn/kmc/include
        )



add_executable(${PROJECT_NAME} Kn_Main.cpp ${PROJECT_SRC})
target_include_directories(${PROJECT_NAME} PUBLIC ${PROJECT_INCLUDE})
set_target_properties(${PROJECT_NAME} PROPERTIES OUTPUT_NAME "kn.exe")

#set(Boost_NO_SYSTEM_PATHS TRUE)
#if (Boost_NO_SYSTEM_PATHS)
#    set(BOOST_ROOT "/opt/apps/intel18/impi18_0/boost-mpi/1.68")
#    set(BOOST_INCLUDE_DIRS "${BOOST_ROOT}/include")
#    set(BOOST_LIBRARY_DIRS "${BOOST_ROOT}/lib")
#endif (Boost_NO_SYSTEM_PATHS)
find_package(Boost COMPONENTS mpi serialization REQUIRED)
target_include_directories(${PROJECT_NAME} PUBLIC ${Boost_INCLUDE_DIR})
target_link_libraries(${PROJECT_NAME} PUBLIC ${Boost_LIBRARIES})
##############################################################################################

#find_package(EIGEN3 REQUIRED)
#include_directories(${EIGEN3_INCLUDE_DIR})


#find_package(TensorFlow REQUIRED)
#target_include_directories(${PROJECT_NAME} PUBLIC ${TensorFlow_INCLUDE_DIRS})
#target_link_libraries(${PROJECT_NAME} PUBLIC ${TensorFlow_LIBRARIES})

#find_package(MPI REQUIRED)
## include directories
#target_include_directories(${PROJECT_NAME} INTERFACE ${MPI_C_INCLUDE_PATH} ${MPI_CXX_INCLUDE_PATH})
## link libraries
#target_link_libraries(${PROJECT_NAME} INTERFACE ${MPI_C_LIBRARIES} ${MPI_CXX_LIBRARIES})
## set interface target compile options
#target_compile_options(${PROJECT_NAME} INTERFACE
#        "${MPI_C_COMPILE_DEFINITIONS} ${MPI_CXX_COMPILE_DEFINITIONS}")
#target_compile_options(${PROJECT_NAME} INTERFACE
#        "${MPI_C_COMPILE_OPTIONS} ${MPI_CXX_COMPILE_OPTIONS}")
#set_target_properties(${PROJECT_NAME} PROPERTIES
#        INTERFACE_LINK_OPTIONS "${MPI_C_LINK_FLAGS} ${MPI_CXX_LINK_FLAGS}")

#find_package(OpenMP REQUIRED)
## include directories
#target_include_directories(${PROJECT_NAME} INTERFACE ${OpenMP_C_INCLUDE_DIRS} ${OpenMP_CXX_INCLUDE_DIRS})
## link libraries
#target_link_libraries(${PROJECT_NAME} INTERFACE ${OpenMP_C_LIBRARIES} ${OpenMP_CXX_LIBRARIES})
## set interface target compile options
#target_compile_options(${PROJECT_NAME} INTERFACE
#        "${OpenMP_C_FLAGS} ${OpenMP_CXX_FLAGS}")


#find_package(ARMADILLO REQUIRED)
#target_include_directories(${PROJECT_NAME} PUBLIC ${ARMADILLO_INCLUDE_DIRS})
#target_link_libraries(${PROJECT_NAME} INTERFACE ${ARMADILLO_LIBRARIES})
#
## Lib keras2cpp
#add_subdirectory(external/keras2cpp)
#target_include_directories(${PROJECT_NAME} INTERFACE external/keras2cpp)
#target_link_libraries(${PROJECT_NAME} INTERFACE keras2cpp)


