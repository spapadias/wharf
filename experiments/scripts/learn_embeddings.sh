#!/bin/bash

# script options
clean_after_experiments=True

# execution options
declare -a graphs=("aspen-paper-graph")
declare -a walks_per_node=(1)
declare -a walk_length=(10)

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
make learn_embeddings

## Execute experiments
#for wpv in "${walks_per_node[@]}"; do
#    for wl in "${walk_length[@]}"; do
#        for graph in "${graphs[@]}"; do
#            printf "\n \n \n"
#            echo "Executing experiment with parameters: " "${graph}" "${wpv}" "${wl}"
#            time ./learn_embeddings -s -f "data/${graph}.adj" -w "${wpv}" -l "${wl}"
#            #    numactl -i all ./memory_footprint -s -f "data/${graph}.adj" -w "${wpn}" -l "${wl}" todo: check numactl for parallel execution. read aspen readme
#        done
#    done
#done

# clean build if necessary
if [ "$clean_after_experiments" = True ] ; then
    # delete the build
    cd ../../;
    rm -rf cmake-build;
    rm experiments/data/*.adj
fi