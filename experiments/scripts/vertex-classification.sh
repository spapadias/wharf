#!/bin/bash

# script options
clean_build=False                    # removes build folder after the execution

# execution options
model="deepwalk"                     # deepwalk | node2vec
p=4.0                                # node2vec's p
q=1.0                                # node2vec's q
init="weight"                        # random | weight | burnin
vector_dimension=128                 # size of learned vectors
learning_strategy=1                  # 1: online | 2: mini-batch (default)
edge_parition_size=250               # size of the edges parition
declare -a graphs=("cora-graph")     # array of graphs
declare -a w=(10)                    # walks per vertex to generate
declare -a l=(80)                    # length of one walk

# 1. convert graphs in adjacency graph format if necessary
for graph in "${graphs[@]}"; do
  FILE=../data/"${graph}".adj
  if test -f "$FILE"; then
      echo 'Skipping conversion of a graph' "${graph[@]}" 'to adjacency graph format'
  else
      echo 'Converting graph' "${graph[@]}" 'to adjacency graph format'
      ./../bin/SNAPtoAdj -s -f "../data/${graph}" "../data/${graph}.adj"
      echo 'Graph' "${graph[@]}" 'converted to adjacency graph format'
  fi
done

# 2. build yskip
bash install-yskip.sh

# 3. build the executable
mkdir -p ../../build;
cd ../../build;
cmake -DCMAKE_BUILD_TYPE=Release ..;
cd experiments;

# build python virtual environment
pip3 install virtualenv
virtualenv venv
source venv/bin/activate
pip3 install pandas scikit-learn

# build executable
make vertex-classification

# 4. execute experiments
for wpv in "${w[@]}"; do
    for wl in "${l[@]}"; do
        for graph in "${graphs[@]}"; do
          printf "\n"
          printf "${graph}"
          ./vertex-classification -s -f "data/${graph}" -w "${w}" -l "${l}" -model "${model}" -p "${p}" -q "${q}" -init "${init}" -d "${vector_dimension}" -le "${learning_strategy}" -eps "${edge_parition_size}"
        done
    done
done

# 5. clean build if necessary
if [ "$clean_build" = True ] ; then
    cd ../../;
    rm -rf build;
    rm -rf yskip;
    rm experiments/data/*.adj
fi
