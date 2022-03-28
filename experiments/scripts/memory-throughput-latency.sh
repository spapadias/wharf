#!/bin/bash

# script options
clean_build=True

# execution options
walk_model="deepwalk"             # deepwalk | node2vec
paramP=0.5                        # node2vec paramP
paramQ=2.0                        # node2vec paramQ
sampler_init_strategy="weight"    # random | burnin | weight
declare -a graphs=("flickr-graph")
declare -a walks_per_node=(10)
declare -a walk_length=(80)
range_search="true"               # range search mode
determinism="true"                # determinism
num_of_batches=10                 # numbers of batches
half_of_batch_size=3500           # batch_size / 2
merge_wu_exec_mode="parallel"     # parallel | serial
merge_frequency=10                # every how many batches to merge

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
            printf "\n"
            printf "Graph: ${graph} \n"
            ./memory-throughput-latency -s -f "data/${graph}.adj" -w "${wpv}" -l "${wl}" -model "${walk_model}" -paramP "${paramP}" -paramQ "${paramQ}" -init "${sampler_init_strategy}" -rs "${range_search}" -d "${determinism}" -numbatch "${num_of_batches}" -sizebatch "${half_of_batch_size}" -mergefreq "${merge_frequency}$" -mergemode "${merge_wu_exec_mode}"
        done
    done
done

# 4. clean build if necessary
if [ "$clean_build" = True ] ; then
    cd ../../;
    rm -rf build;
#    rm experiments/data/*.adj
fi
