# -------------------------- DYNAMIC GRAPH REPRESENTATION LEARNING WITH METROPOLIS HASTINGS ----------------------------
cmake_minimum_required(VERSION 3.16)
project(dynamic_graph_representation_learning_with_metropolis_hastings)

# libraries
include_directories(libs/libcuckoo)
include_directories(libs/compressed_trees)

# include paths
include_directories(graph)
include_directories(walks)
include_directories(utils)
include_directories(config)

# set up google tests
enable_testing()
add_subdirectory(tests)

# experiments
add_subdirectory(experiments)
# -------------------------- DYNAMIC GRAPH REPRESENTATION LEARNING WITH METROPOLIS HASTINGS ----------------------------