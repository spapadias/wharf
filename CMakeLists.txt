# -------------------------- Wharf ----------------------------
cmake_minimum_required(VERSION 3.13)
project(wharf)

# libraries
include_directories(libs/libcuckoo)
include_directories(libs/compressed_trees)
include_directories(libs/cds)

# include paths
include_directories(graph)
include_directories(walks)
include_directories(utils)
include_directories(config)

# set up google tests
#enable_testing()
#add_subdirectory(tests)

# experiments
add_subdirectory(experiments)
# -------------------------- Wharf ----------------------------