cmake_minimum_required(VERSION 2.6)

project(optLib CXX)

set (CMAKE_CXX_STANDARD 11)

find_package(OOQP REQUIRED)

include_directories(${OOQP_INCLUDE_DIRS})

FILE(GLOB SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/*.h ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp)

add_library(${PROJECT_NAME} ${SOURCES})
set_target_properties(${PROJECT_NAME} PROPERTIES LINKER_LANGUAGE CXX)
target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(
        ${PROJECT_NAME}
	    ${OOQP_LIBRARIES}
)
