#!/bin/bash

# script options
clean_build=True

# execution options
walk_model="deepwalk"                             # deepwalk | node2vec
paramP=4.0                                        # node2vec paramP
paramQ=1.0                                        # node2vec paramQ
sampler_init_strategy="weight"                    # random | burnin | weight
declare -a graphs=("flickr-graph")
declare -a walks_per_node=(10)
declare -a walk_length=(80)
determinism="true"

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
make memory-throughput-latency

# 3. execute experiments
for wpv in "${walks_per_node[@]}"; do
    for wl in "${walk_length[@]}"; do
        for graph in "${graphs[@]}"; do
            printf "${graph}"
            ./memory-throughput-latency -s -f "data/${graph}.adj" -w "${wpv}" -l "${wl}" -model "${walk_model}" -paramP "${paramP}" -paramQ "${paramQ}" -init "${sampler_init_strategy}" -det "${determinism}"
        done
    done
done

# 4. clean build if necessary
if [ "$clean_build" = True ] ; then
    cd ../../;
    rm -rf build;
#    rm experiments/data/*.adj
fi
