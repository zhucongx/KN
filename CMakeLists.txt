cmake_minimum_required(VERSION 3.11)
########################################## SET  COMPILERS ##########################################
#set(CMAKE_C_COMPILER "mpicc") # {gcc, clang, mpicc, icc}
#set(CMAKE_CXX_COMPILER "mpicxx") # {g++, clang++, mpicxx, icpc}
project(KN)
set(CMAKE_CXX_STANDARD 20)
############################################ SET  PATHS ############################################
list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin)
set(LIBRARY_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/lib)
############################################ SET  FLAGS ############################################
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
############################################ BUILD LIBS ############################################
#add_executable(${PROJECT_NAME} Kn_Main.cpp)
#add_subdirectory(kn/cfg)
#add_subdirectory(kn/gen)
#add_subdirectory(kn/ansys)
#add_subdirectory(kn/kmc)
#target_link_libraries(${PROJECT_NAME} PUBLIC cfg gen ansys kmc)
############################################ BUILD SRCS ############################################
file(GLOB PROJECT_SRC
        ${PROJECT_SOURCE_DIR}/kn/cfg/src/*.cpp
        ${PROJECT_SOURCE_DIR}/kn/gen/src/*.cpp
        ${PROJECT_SOURCE_DIR}/kn/ansys/src/*.cpp
        ${PROJECT_SOURCE_DIR}/kn/kmc/src/*.cpp
        ${PROJECT_SOURCE_DIR}/kn/api/src/*.cpp
        )
set(PROJECT_INCLUDE
        ${PROJECT_SOURCE_DIR}/kn/cfg/include
        ${PROJECT_SOURCE_DIR}/kn/gen/include
        ${PROJECT_SOURCE_DIR}/kn/ansys/include
        ${PROJECT_SOURCE_DIR}/kn/kmc/include
        ${PROJECT_SOURCE_DIR}/kn/api/include
        )
add_executable(${PROJECT_NAME} Kn_Main.cpp ${PROJECT_SRC})
target_include_directories(${PROJECT_NAME} PUBLIC ${PROJECT_INCLUDE})
set_target_properties(${PROJECT_NAME} PROPERTIES OUTPUT_NAME "kn.exe")
############################################## BOOST ###############################################
#set(Boost_NO_SYSTEM_PATHS TRUE)
#if (Boost_NO_SYSTEM_PATHS)
#    set(BOOST_ROOT "/opt/apps/intel18/impi18_0/boost-mpi/1.68")
#    set(BOOST_INCLUDE_DIRS "${BOOST_ROOT}/include")
#    set(BOOST_LIBRARY_DIRS "${BOOST_ROOT}/lib")
#endif (Boost_NO_SYSTEM_PATHS)
# Stampede2 doesn't have boost-mpi. Use self-compiled version
#SET(CMAKE_INCLUDE_PATH ${CMAKE_INCLUDE_PATH} "/home1/07075/tg863748/library/boost/include")
#SET(CMAKE_LIBRARY_PATH ${CMAKE_LIBRARY_PATH} "/home1/07075/tg863748/library/boost/lib")
find_package(Boost COMPONENTS mpi serialization REQUIRED)
target_include_directories(${PROJECT_NAME} PRIVATE ${Boost_INCLUDE_DIR})
target_link_libraries(${PROJECT_NAME} PRIVATE ${Boost_LIBRARIES})
############################################### JSON ###############################################
# Typically you don't care so much for a third party library's tests to be
# run from your own project's code.
#set(JSON_BuildTests OFF CACHE INTERNAL "")
# If you only include this third party in PRIVATE source files, you do not
# need to install it when your main project gets installed.
# set(JSON_Install OFF CACHE INTERNAL "")
#add_subdirectory(external/json)
include(FetchContent)
FetchContent_Declare(json
        GIT_REPOSITORY https://github.com/ArthurSonzogni/nlohmann_json_cmake_fetchcontent
        GIT_TAG v3.9.1)

FetchContent_GetProperties(json)
if(NOT json_POPULATED)
    FetchContent_Populate(json)
    add_subdirectory(${json_SOURCE_DIR} ${json_BINARY_DIR} EXCLUDE_FROM_ALL)
endif()
target_link_libraries(${PROJECT_NAME} PRIVATE nlohmann_json::nlohmann_json)
############################################### MPI ################################################
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
############################################## OPENMP ##############################################
#find_package(OpenMP REQUIRED)
## include directories
#target_include_directories(${PROJECT_NAME} INTERFACE ${OpenMP_C_INCLUDE_DIRS} ${OpenMP_CXX_INCLUDE_DIRS})
## link libraries
#target_link_libraries(${PROJECT_NAME} INTERFACE ${OpenMP_C_LIBRARIES} ${OpenMP_CXX_LIBRARIES})
## set interface target compile options
#target_compile_options(${PROJECT_NAME} INTERFACE
#        "${OpenMP_C_FLAGS} ${OpenMP_CXX_FLAGS}")
############################################ GOOGLETEST ############################################
#include(FetchContent)
#FetchContent_GetProperties(googletest)
#if (NOT googletest_POPULATED)
#    FetchContent_Populate(googletest)
#    add_subdirectory(${googletest_SOURCE_DIR} ${googletest_BINARY_DIR})
#endif ()
############################################ ARMADILLO #############################################
#find_package(ARMADILLO REQUIRED)
#target_include_directories(${PROJECT_NAME} PUBLIC ${ARMADILLO_INCLUDE_DIRS})
#target_link_libraries(${PROJECT_NAME} PUBLIC ${ARMADILLO_LIBRARIES})
############################################## EIGEN ###############################################
#find_package(EIGEN3 REQUIRED)
#target_include_directories(${PROJECT_NAME} PUBLIC ${EIGEN3_INCLUDE_DIRS})
############################################ TENSORFLOW ############################################
#find_package(TensorFlow REQUIRED)
#target_include_directories(${PROJECT_NAME} PUBLIC ${TensorFlow_INCLUDE_DIRS})
#target_link_libraries(${PROJECT_NAME} PUBLIC ${TensorFlow_LIBRARIES})



