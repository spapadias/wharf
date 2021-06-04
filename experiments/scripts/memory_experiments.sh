#!/bin/bash

# script options
clean_after_experiments=False

# execution options
walk_model="deepwalk" # deepwalk | node2vec
paramP=0.5  # node2vec paramP
paramQ=0.5  # node2vec paramQ
sampler_init_strategy="random" #  random | burnin | weight
declare -a graphs=("email-graph")
declare -a walks_per_node=(10)
declare -a walk_length=(80)

# convert graphs in adjacency graph format if necessary
for graph in "${graphs[@]}"; do
  FILE=../data/"${graph}".adj
  if test -f "$FILE"; then
      echo 'Skipping conversion of a graph ' "${graph[@]}" ' to adjacency graph format!'
  else
      echo 'Converting graph ' "${graph[@]}" ' to adjacency graph format ...'
      ./SNAPtoAdj -s -f "../data/${graph}" "../data/${graph}.adj"
      echo 'Graph ' "${graph[@]}" ' converted to adjacency graph format!'
  fi
done

# create the build directory
mkdir -p ../../cmake-build;
cd ../../cmake-build;
cmake ..;
cd experiments;

# build the learn embeddings experiment
make memory_experiments

# execute experiments
for wpv in "${walks_per_node[@]}"; do
    for wl in "${walk_length[@]}"; do
        for graph in "${graphs[@]}"; do
            printf "\n \n"
            echo "Executing experiment with parameters: " "${graph}" "${wpv}" "${wl}"
            time ./memory_experiments -s -f "data/${graph}.adj" -w "${wpv}" -l "${wl}" -model "${walk_model}" -paramP "${paramP}" -paramQ "${paramQ}" -init "${sampler_init_strategy}"
        done
    done
done

# clean build if necessary
if [ "$clean_after_experiments" = True ] ; then
    # delete the build
    cd ../../;
    rm -rf cmake-build;
    rm experiments/data/*.adj
fi