project(ansys)

file(GLOB PROJECT_SRC ${PROJECT_SOURCE_DIR}/src/*.cpp)
set(PROJECT_INCLUDE "${PROJECT_SOURCE_DIR}/include")

add_library(${PROJECT_NAME} ${PROJECT_SRC})
target_include_directories(${PROJECT_NAME} PUBLIC ${PROJECT_INCLUDE})
set_target_properties(${PROJECT_NAME} PROPERTIES OUTPUT_NAME "ansys")

target_link_libraries(${PROJECT_NAME} PUBLIC cfg)

find_package(Boost COMPONENTS mpi serialization REQUIRED)
target_include_directories(${PROJECT_NAME} PUBLIC ${Boost_INCLUDE_DIR})
target_link_libraries(${PROJECT_NAME} PUBLIC ${Boost_LIBRARIES})