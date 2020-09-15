project(cfg)

file(GLOB PROJECT_SRC ${PROJECT_SOURCE_DIR}/src/*.cpp)
set(PROJECT_INCLUDE "${PROJECT_SOURCE_DIR}/include")

add_library(${PROJECT_NAME} ${PROJECT_SRC})
target_include_directories(${PROJECT_NAME} PUBLIC ${PROJECT_INCLUDE})
set_target_properties(${PROJECT_NAME} PROPERTIES OUTPUT_NAME "cfg")

find_package(Boost REQUIRED)
target_include_directories(${PROJECT_NAME} PUBLIC ${Boost_INCLUDE_DIR})



#[[
# test binary file
add_executable(test test.cpp ${PROJECT_SRC})
target_include_directories(test PUBLIC ${PROJECT_INCLUDE})
set_target_properties(test PROPERTIES OUTPUT_NAME "config_tools.exe")
target_include_directories(test PUBLIC ${Boost_INCLUDE_DIR})
]]
