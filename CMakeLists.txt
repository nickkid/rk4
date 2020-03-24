cmake_minimum_required(VERSION 2.8)

project(rk4)

aux_source_directory(. DIR_SRCS)
add_executable(rk4 ${DIR_SRCS}) 
target_link_libraries(rk4 m)