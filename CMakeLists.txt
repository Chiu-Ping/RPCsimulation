CMAKE_MINIMUM_REQUIRED(VERSION 3.10 FATAL_ERROR)
cmake_policy(SET CMP0048 NEW)
project(RPCsimulation VERSION 3.10 LANGUAGES C CXX)
add_executable(gasfile gasfile.C)
target_link_libraries(gasfile ${ROOT_LIBRARIES} Garfield)
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})

