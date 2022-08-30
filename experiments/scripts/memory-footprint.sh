#!/bin/bash

# script options
clean_build=True               # cleans build folder after the execution

# execution options
model="deepwalk"                 # deepwalk | node2vec
p=0.5                            # node2vec's p parameter
q=2.0                            # node2vec's q parameter
init="random"                    # random | burnin | weight
declare -a graphs=("cora-graph") # array of graphs (can intake more than one dataset)
declare -a w=(10)                # walks per vertex to generate
declare -a l=(80)                # length of one walk
rs="true"                        # range search mode: ON | OFF
det="true"                       # determinism in random walking via predefined seeds (for debugging)
nb=1                             # numbers of batches
bs=1000                          # half the size of the batch size
merge_mode="parallel"            # mode for conducting the merge operation of the walk-trees: parallel | serial
merge_freq=1                     # every after "how many" batches to conduct the merging operation

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
make memory-footprint

# 3. execute experiments
for wpv in "${w[@]}"; do
    for wl in "${l[@]}"; do
        for graph in "${graphs[@]}"; do
            printf "\n"
            printf "Graph: ${graph}\n"
            ./memory_footprint -s -f "data/${graph}.adj" -w "${w}" -l "${l}" -model "${model}" -p "${p}" -q "${q}" -init "${init}" -rs "${rs}" -d "${det}" -nb "${nb}" -bs "${bs}" -merge_freq "${merge_freq}$" -merge_mode "${merge_mode}"
        done
    done
done

# 4. clean build if necessary
if [ "$clean_build" = True ] ; then
    cd ../../;
    rm -rf build;
     rm experiments/data/*.adj
fi
