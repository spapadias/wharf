#!/bin/bash

# script options
clean_build=True                    # removes build folder after the execution

# execution options
walk_model="deepwalk"                # deepwalk | node2vec
paramP=1.0                           # node2vec's paramP
paramQ=4.0                           # node2vec's paramQ
sampler_init_strategy="weight"       # random | weight | burnin
vector_dimension=128                 # size of learned vectors
learning_strategy=1                  # 1: online | 2: mini-batch (default)
edge_parition_size=500               # size of the edges parition
declare -a graphs=("cora-graph")     # array of graphs
declare -a walks_per_vertex=(10)     # walks per vertex to generate
declare -a walk_length=(80)          # length of one walk

# TEMP AREA
num_of_lines=$(< "../data/${graphs[0]}" wc -l)
head -n 2714 "../data/${graphs[0]}" > "../data/preload"
tail -n +2714 "../data/${graphs[0]}" > "../data/stream"

# 1. convert graphs in adjacency graph format if necessary
for graph in "${graphs[@]}"; do
  FILE=../data/preload.adj
  if test -f "$FILE"; then
      echo 'Skipping conversion of a graph' "${graph[@]}" 'to adjacency graph format'
  else
      echo 'Converting graph' "${graph[@]}" 'to adjacency graph format'
      ./../bin/SNAPtoAdj -s -f "../data/preload" "../data/preload.adj"
      echo 'Graph' "${graph[@]}" 'converted to adjacency graph format'
  fi
done

# 2. build the executable
mkdir -p ../../build;
cd ../../build;
cmake -DCMAKE_BUILD_TYPE=Release ..;
cd experiments;

# build python virtual environment
sudo pip3 install virtualenv
virtualenv venv
source venv/bin/activate
pip3 install pandas scikit-learn

# build executable
make vertex-classification

# 3. execute experiments
for wpv in "${walks_per_vertex[@]}"; do
    for wl in "${walk_length[@]}"; do
        for graph in "${graphs[@]}"; do
          printf "\n"
          printf "${graph}"
          ./vertex-classification -s -f "data/preload" -w "${wpv}" -l "${wl}" -model "${walk_model}" -paramP "${paramP}" -paramQ "${paramQ}" -init "${sampler_init_strategy}" -d "${vector_dimension}" -le "${learning_strategy}" -eps "${edge_parition_size}"
        done
    done
done

# 4. clean build if necessary
if [ "$clean_build" = True ] ; then
    cd ../../;
    rm -rf build;
    rm -rf yskip;
    rm experiments/data/*.adj
    rm experiments/data/preload*
    rm experiments/data/stream*
fi