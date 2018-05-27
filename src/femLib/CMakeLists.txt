cmake_minimum_required(VERSION 2.6)

project(femLib CXX)

set (CMAKE_CXX_STANDARD 11)
FILE(GLOB SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/*.h ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp)

include_directories(../optLib)
add_library(${PROJECT_NAME} ${SOURCES})
target_link_libraries(${PROJECT_NAME})
target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
