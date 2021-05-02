# -------------------------- DYNAMIC GRAPH REPRESENTATION LEARNING WITH METROPOLIS HASTINGS ----------------------------
cmake_minimum_required(VERSION 3.16)
project(dynamic_graph_representation_learning_with_metropolis_hastings)

# dependencies
include_directories(libs/parlay)
include_directories(libs/compressed_trees)

# internal includes
include_directories(graph)
include_directories(walks)
include_directories(config)
#include_directories(utils)

# set up google tests
enable_testing()
add_subdirectory(tests)

# experiments
#add_subdirectory(experiments)

#option(CLANG     "???" OFF)
#option(CILK      "???" OFF)
#option(OPENMP    "???" OFF)
#option(HOMEGROWN "???" OFF)
#option(SERIAL    "???" OFF)

## compile flags options
#set(OMPFLAGS     "-DOPENMP -fopenmp")
#set(CILKFLAGS    "-DCILK -fcilkplus")
#set(HGFLAGS      "-DHOMEGROWN -pthread")
#set(OPTIMIZATION "-O3 -g")
## https://stackoverflow.com/questions/18656433/the-assert-statement-does-not-work-in-c
#set(CFLAGS       "-DEDGELONG -mcx16 ${OPTIMIZATION} -ldl -std=c++17 -march=native -Wall -Wno-subobject-linkage -DUSEMALLOC ") # -DNDEBUG -DVERTEXLONG

#if(CLANG)
#    set(CMAKE_CXX_COMPILER "clang++")
#    set(PFLAGS ${CILKFLAGS})
#elseif(CILK)
#    set(CMAKE_CXX_COMPILER "g++")
#    set(PFLAGS ${CILKFLAGS})
#elseif(OPENMP)
#    set(CMAKE_CXX_COMPILER "g++")
#    set(PFLAGS ${OMPFLAGS})
#elseif(HOMEGROWN)
#    set(CMAKE_CXX_COMPILER "g++")
#    set(PFLAGS ${HGFLAGS})
#elseif(SERIAL)
#    set(CMAKE_CXX_COMPILER "g++")
#    set(PFLAGS "")
#else()
#    set(CMAKE_CXX_COMPILER "g++")
#    set(PFLAGS ${HGFLAGS})
#endif()

## set compile flags
#set(CMAKE_CXX_FLAGS "${CFLAGS} ${PFLAGS}")
#
## display build parameters
#message("Build parameters:")
#message("  -- Compiler: "                  ${CMAKE_CXX_COMPILER})
#message("  -- Compiler flags (base): "     ${CFLAGS})
#message("  -- Compiler flags (parallel): " ${PFLAGS})
# -------------------------- DYNAMIC GRAPH REPRESENTATION LEARNING WITH METROPOLIS HASTINGS ----------------------------