#!/bin/bash

# script options
clean_build=False                           # cleans build folder after the execution

# execution options
walk_model="deepwalk"                       # deepwalk | node2vec
paramP=0.2                                  # node2vec's paramP
paramQ=0.7                                  # node2vec's paramQ
sampler_init_strategy="random"              # random | burnin | weight
vector_dimension=5                          # size of learned vectors
learning_strategy=2                         # 1: online | 2: mini-batch (default)
declare -a graphs=("aspen-paper-graph")     # array of graphs
declare -a walks_per_vertex=(2)             # walks per vertex to generate
declare -a walk_length=(10)                 # length of one walk

# 1. convert graphs in adjacency graph format if necessary
for graph in "${graphs[@]}"; do
  FILE=../data/"${graph}".adj
  if test -f "$FILE"; then
      echo 'Skipping conversion of a graph ' "${graph[@]}" ' to adjacency graph format!'
  else
      echo 'Converting graph ' "${graph[@]}" ' to adjacency graph format ...'
      ./../bin/SNAPtoAdj -s -f "../data/${graph}" "../data/${graph}.adj"
      echo 'Graph ' "${graph[@]}" ' converted to adjacency graph format!'
  fi
done

# 2. build the executable
mkdir -p ../../build;
cd ../../build;
cmake -DCMAKE_BUILD_TYPE=Release ..;
cd experiments;
make vertex-classification

# 3. execute experiments
for wpv in "${walks_per_vertex[@]}"; do
    for wl in "${walk_length[@]}"; do
        for graph in "${graphs[@]}"; do
            time ./vertex-classification -s -f "data/${graph}.adj" -w "${wpv}" -l "${wl}" -model "${walk_model}" -paramP "${paramP}" -paramQ "${paramQ}" -init "${sampler_init_strategy}" -d "${vector_dimension}" -le "${learning_strategy}"
        done
    done
done

# 4. clean build if necessary
if [ "$clean_build" = True ] ; then
    cd ../../;
    rm -rf build;
    rm experiments/data/*.adj
fi