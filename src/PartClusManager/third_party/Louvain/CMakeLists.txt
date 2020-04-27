
cmake_minimum_required (VERSION 3.1)

project (Louvain)

set(LIB_NAME "Louvain")

set(CMAKE_CXX_STANDARD 11)

###########################################################
# Source files
###########################################################

set(SRC_FILES
        ${CMAKE_CURRENT_SOURCE_DIR}/src/balmod.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/condora.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/devind.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/devuni.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/dp.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/goldberg.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/graph.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/graph_binary.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/louvain.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/main_convert.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/main_hierarchy.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/main_louvain.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/main_matrix.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/modularity.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/owzad.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/quality.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/shimalik.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/zahn.cpp
   )

###########################################################
# Library
###########################################################

add_definitions(-DLIBRARY_MODE)
add_LIBRARY(${LIB_NAME} STATIC ${SRC_FILES})

target_include_directories(${LIB_NAME} PUBLIC
                           ${CMAKE_CURRENT_SOURCE_DIR}/src)

target_link_libraries(${LIB_NAME})
