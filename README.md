# Wharf: Space-Efficient Random Walks on Streaming Graphs

This repository provides the code, scripts, and instructions for reproducing the experiments in our paper, Space-Efficient Random Walks on Streaming Graphs [1]. Our paper introduces Wharf, a main-memory streaming random walk system capable of updating, indexing, and compressing a maintained random walk corpus efficiently over streaming (i.e., highly dynamic) graphs that are evolving.

Wharf is designed for maintaining a dynamic graph that is subject to updates by a single writer, while supporting multiple concurrent readers by utilizing the C-tree structure introduced in [2]. Wharf represents the random walks in the form of triplets, which encodes as unique integers, and stores them in C-trees for achieving both high throughput and low memory footprint.

In the Getting Started Guide, we include functionality to reproduce the main results presented in our Wharf paper.

# Getting Started Guide

This Getting Started Guide gives an overview of

1. Using Wharf as a Random Walk Streaming System

2. Setting up Wharf

- Hardware and software dependencies
- Input formats
- Obtaining datasets

3. Experiment Workflow
- Memory Footprint in Wharf
- Throughput and Latency in Wharf
- Vertex Classification with Wharf's walks

4. Miscellaneous
- Wharf++
- References

# Using Wharf as a Streaming Random Walk System

We give a brief overview of the user-level API provided by Wharf for conducting and updating random walks on evolving graph workloads. Currently, we support DeepWalk and node2vec, yet we can support any random walking model.

```
/* Creates initial set of random walks */
generate_initial_random_walks() -> void

/* Inserts a batch of edges in the graph and updates the maintained walk corpus */
insert_edges_batch(size_t m, std::tuple<uintV, uintV>* edges, int batch_num, bool sorted = false, bool remove_dups = false, size_t nn = std::numeric_limits<size_t>::max(), bool apply_walk_updates = true, bool run_seq = false) -> pbbs::sequence<types::WalkID>

/* Deletes a batch of edges in the graph and updates the maintained walk corpus */
delete_edges_batch(size_t m, std::tuple<uintV, uintV>* edges, int batch_num, bool sorted = false, bool remove_dups = false, size_t nn = std::numeric_limits<size_t>::max(), bool apply_walk_updates = true, bool run_seq = false) -> pbbs::sequence<types::WalkID>
```

The pairings.h class has the following interface:

```
/* Maps pair of numbers to a unique number using the Szudzik pairing function */
template <class Type>
pair(const std::pair<Type, Type>& pair) -> Type

/* Converts a unique number produced by Szudzik function to a unique pair of numbers */
template <class Type>
unpair(const Type& encoded_value) -> std::pair<Type, Type>
```

The random_walk_model.h type has the following interface:

```
// State = the state is defined as a pair of two numbers, where the first represents the current vertex 
// and the second contains the extra information (e.g DeepWalk = current vertex, Node2Vec = previously visited vertex by the walker)
using State = std::pair<Vertex, Vertex>;

/* Determines an initial state of the walker */
initial_state(types::Vertex vertex) -> types::State

/* A new state of the walker */
new_state(const types::State& state, types::Vertex vertex) -> types::State
 
/* Calculates the edge weight based on the current state and the potentially proposed vertex */
weight(const types::State& state, types::Vertex vertex) -> float
 
/* Propose next vertex given current state */
propose_vertex(const types::State& state) -> types::Vertex
```

The metropolis_hastings_sampler.h type has the following interface:

```
/* MetropolisHastingsSampler constructor */
MetropolisHastingsSampler(types::State& current_state, dygrl::RandomWalkModel* model)

/* Sample new vertex */
sample(types::State& state, dygrl::RandomWalkModel* model) -> types::State
 
/* Previously sampled vertex */
previously_sampled_vertex() -> types::Vertex
 
/* Metropolis Hastings sampler initializer */
init(types::State& current_state, dygrl::RandomWalkModel* model, types::SamplerInitStartegy init_startegy) -> void
 
/* Accept candidate */
accept(float previous_weight, float new_weight) -> bool
```

Users can implement and run any random walking algorithm in Wharf due to the Metropolis-Hastings sampling approach utilized that can approximate and provably converge to any random walk model distribution. Currently, Wharf can generate DeepWalk random walks, i.e., first-order walks, and node2vec walks, i.e., second-order walks. The actual DeepWalk and node2vec instances can be found in the following files located in `walks/models/`:

```
deepwalk.h /* DeepWalk random walks */
node2vec.h /* node2vec random walks */
```

# Setting up Wharf

## Hardware Dependencies

Any modern x86-based multicore machine can be used. Wharf uses 128-bit CMPXCHG (requires -mcx16 compiler flag) but does not need hardware transactional memory (TSX). The experiments on all the graph datasets we used demand fewer than 1.5GB of main memory for storing both the graph and the maintain random walks corpus in memory. Note that the walk corpus is multiple times larger than the maintained graph, especially in Graph Representation Learnling (GRL) applications.

## Software Dependencies

Wharf requires g++ 7.5.0 or later versions supporting the Cilk Plus extensions. Note that all the experiments can by run using `numactl -i all` for better performance, yet we conducted our experiments directly without `numactl`.

## Input Formats

Wharf currently supports reading two formats: the `adjacency graph` and the `compressed graph` format. The former format is used by the Problem Based Benchmark suite, while the latter one was developed as part of the Graph Based Benchmark Suite.

The adjacency graph format starts with a sequence of offsets one for each vertex, followed by a sequence of directed edges ordered by their source vertex. The offset for a vertex i refers to the location of the start of a contiguous block of out edges for vertex i in the sequence of edges. The block continues until the offset of the next vertex, or the end if i is the last vertex. All vertices and offsets are 0 based and represented in decimal. The format we described above is depicted in the following:

```
Adjacency Graph Format 
<n>
<m>
<o0>
<o1>
...
<o(n-1)>
<e0>
<e1>
...
<e(m-1)>
```

This file is represented as plain text. The compressed format is the `bytePDA` format, which is similar to the `parallelByte` format of [Ligra+](https://github.com/jshun/ligra), extended with additional functionality. Note that for this prototype, we limit Wharf to processing symmetric, unweighted graph workloads. We plan to extend Wharf for a journal submission and enhance its functionality such that it can ultimately support both undirected and directed (weighted) graphs.

## Obtaining the Datasets

All the real graph workloads we used in our experiments can be obtained from the extremely helpful [graph dataset repository](https://renchi.ac.cn/datasets/) of Dr. [Renchi Yang](https://renchi.ac.cn/). We recommend using the soc-LiveJournal graph. Once one downloads a graph datsets, we provide the `SNAPtoAdj` executable, which is located at `experiments/bin/`, to symmetrize it, and store it in the the format(s) used by Wharf. The `SNAPToAdj` software is taken from the the [Ligra](https://github.com/jshun/ligra) repository, namely, from `ligra/utils/SNAPToAdj.C`.

For the synthetic datasets we used the [TrillionG](https://chan150.github.io/TrillionG/) which is the state of the art graph generator published in SIGMOD 2017. TrillionG uses the RMAT model to generate the a graph and has essentially four parameters a,b,c, and d. For producing each of the synthetic graphs please check the `Datasets` description in our paper.

Note that currently in our repo we only provide cora dataset (along with its label for the vertex classification experiment) for the sake of space.

# Experiment Workflow

This section assumes that you have downloaded (some of) the input graphs listed above to the `experiments/data/` folder.

The experiments presented in this guide share a common set of flags, which we list below:

```
-f filename     : provides a location to an input graph file
 
-s              : indicates that the input graph is symmetric
 
-c              : indicates that the input graph is compressed
 
-m              : indicates that the input graph should be mmap'd

-model          : choice for random walk model; currently, deepwalk and node2vec (others can be supported as well)
 
-w              : indicates the number of walks per vertex in the graph
 
-l              : indicates the length of each random walk
 
-p              : in case of node2vec, indicates the choice for the parameter p
 
-q              : in case of node2vec, indicates the choice for the parameter q

-init           : initialization strategy for the Metropolis-Hastings sampling

-nb             : indicates the number of batches to be inserted or deleted

-bs             : indicates the number of edges per batch

-rs             : range search mode: ON | OFF

-det            : determinism in random walking via predefined seeds (for debugging)

-merge_mode     : mode for conducting the merge operation of the walk-trees: parallel | serial

-merge_freq     : every after "how many" batches to conduct the merging operation
```

First, the user has to build `Wharf` using the following commands:

```
git clone https://github.com/spapadias/wharf.git
cd wharf/
mkdir build
cd build/
cmake ..
make
```
After the build, all the executables, namely, `memory-footprint`, `throughput-latency`, and `vertex-classification` are located in the folder `build/experiments/`. One may run the executables with a set of parameters or use the scripts we provide in the `experiments/scripts/`.

## Memory Footprint in Wharf

The memory usage of our codes can be measured using a tool called `memory-footprint.sh`, which loads a graph using the C-tree data structure and outputs the number of bytes used by the C-tree representation. Specifically, it separately measures the size for the vertex-tree, for the edge-tree, and for the walk-tree (in MBs and GBs). Furthermore, it calculates the size for storing the min and max id of the next vertex in each walk-tree of the hybrid-tree (Figure 4 in our paper), as well as the size to store all the Metropolis-Hastings samplers, which are both negligible compared to the size of the random walk corpus that is maintained and indexed in main memory. Note that in our experiments we set the probability of a node being selected as a head to p = 1/256, so the expected number of nodes in the edges tree is p*m.

The easiest way we can reproduce the memory footprint experimental results in Figure 8 of [1] is to open, select the appropriate parameters, and run the `memory-footprint.sh` script located in `experiments/scripts/`. However, the same results can be obtained after running the executable `memory-footprint` located in `build/experiments/` for the corresponding datasets and set of parameters. Bellow we demonstrate the output on `cora` dataset for a example set of parameters:

```
$ ./memory-footprint.sh  # run the script from the script/ folder
Graph: cora-graph
Running experiment with: 8 threads.
Parallel Merge and WU
Walks per vertex: 10
Walk length: 80
Walking model: DEEPWALK
Sampler strategy: RANDOM
n = 2708 m = 10556

Graph: 
        Vertices: 2708, Edges: 10556
Vertex Tree: 
        Heads: 2708, Head size: 88, Memory usage: 0.227264 MB = 0.000221938 GB
Edge Trees: 
        Heads: 68, Head size: 48, Lists memory: 0.0865784 MB = 8.45492e-05 GB, Total memory usage: 0.0896912 MB = 8.7589e-05 GB
Walks Trees: 
        Heads: 8408, Head size: 48, Lists memory: 13.3466 MB = 0.0130338 GB, Total memory usage: 13.7315 MB = 0.0134097 GB
Range Search (ON): 
        Total {min,max} memory usage: 0.0413208 MB = 4.03523e-05 GB
Samplers: 
        Total memory usage: 0.0413208 MB = 4.03523e-05 GB
Flat graph: 
        Total memory usage: 0.123962 MB = 0.000121057 GB
Total memory used: 
        14.1311 MB = 0.0137999 GB
```

## Throughput and Latency in Wharf

The throughput and latency of Wharf's random walk updates when batches of insertions and/or deletions arrive can be measured using a tool called `throughput-memory.sh`, which loads a graph into edge-trees, produces a walk corpus and loads it into the walk-trees of the hybrid-tree. Subsequently, it inserts a number of batches of certain number of edge insertions or deletions (we also have an experiment with mixed workload containing both insertions and deletions), and applies both the graph updates and the appropriate random walk updates such that the walk corpus always remains stastistically indistinguishable.

The easiest way we can reproduce the throughput and latency experimental results in Figure 6, and 9 of [1] is to open, select the appropriate parameters, and run the `throughput-latency.sh` script located in `experiments/scripts/`. However, the same results can be obtained after running the executable `throughput-latency` located in `build/experiments/` for the corresponding datasets and set of parameters. Bellow we demonstrate the output on `cora` dataset for a example set of parameters:

```
$ ./throughput-latency.sh  # run the script from the script/ folder
Graph: cora-graph 
Running experiment with 8 threads
Parallel Merge and WU
Walks per vertex: 10
Walk length: 80
Walking model: DEEPWALK
Sampler strategy: WEIGHT
n = 2708 m = 10556
GENERATE 27080 walks.
Total time to generate the walk corpus from scratch: 0.366158
Batch size = 1000

batch-1 and batch_seed-1
(merge) -- For batch-1 we are touching 2708 / 2708 vertices
(insert) -- For batch-1 we are touching 2577 / 2708 vertices
METRICS AT BATCH-1
Insert time (avg) = 0.637079
GUP (avg) = 0.055254
BWUP (avg, includes merge) = 0.565807, average walk affected = 25828, sampled vertices = 1931468
WUP (avg)   = 0.338497  (sampling= 0.29334, inserting= 0.045157)
MAV (avg)   = 0.0519266 MAV (min) = 0.0519266   MAV (max) = 0.0519266
Merge (avg,10 times) = 0.195651 Merge (min) = 0.195651  Merge (max) = 0.195651

batch-2 and batch_seed-2
(merge) -- For batch-2 we are touching 2708 / 2708 vertices
(insert) -- For batch-2 we are touching 2573 / 2708 vertices
METRICS AT BATCH-2
Insert time (avg) = 0.662097
GUP (avg) = 0.050795
BWUP (avg, includes merge) = 0.597204, average walk affected = 25914, sampled vertices = 3902997
WUP (avg)   = 0.356064  (sampling= 0.305324, inserting= 0.05074)
MAV (avg)   = 0.0489409 MAV (min) = 0.0459552   MAV (max) = 0.0519266
Merge (avg,10 times) = 0.213809 Merge (min) = 0.195651  Merge (max) = 0.231967

batch-3 and batch_seed-3
(merge) -- For batch-3 we are touching 2708 / 2708 vertices
(insert) -- For batch-3 we are touching 2584 / 2708 vertices
METRICS AT BATCH-3
Insert time (avg) = 0.681231
GUP (avg) = 0.0515549
BWUP (avg, includes merge) = 0.615161, average walk affected = 26028.7, sampled vertices = 5889032
WUP (avg)   = 0.360533  (sampling= 0.306352, inserting= 0.0541813)
MAV (avg)   = 0.0501107 MAV (min) = 0.0459552   MAV (max) = 0.104377
Merge (avg,10 times) = 0.228055 Merge (min) = 0.195651  Merge (max) = 0.452199

batch-4 and batch_seed-4
(merge) -- For batch-4 we are touching 2708 / 2708 vertices
(insert) -- For batch-4 we are touching 2570 / 2708 vertices
METRICS AT BATCH-4
Insert time (avg) = 0.698353
GUP (avg) = 0.0529832
BWUP (avg, includes merge) = 0.629816, average walk affected = 26094, sampled vertices = 7875110
WUP (avg)   = 0.370472  (sampling= 0.312895, inserting= 0.0575765)
MAV (avg)   = 0.0524014 MAV (min) = 0.0459552   MAV (max) = 0.105229
Merge (avg,10 times) = 0.234821 Merge (min) = 0.195651  Merge (max) = 0.487086

...

batch-9 and batch_seed-9
(merge) -- For batch-9 we are touching 2708 / 2708 vertices
(insert) -- For batch-9 we are touching 2574 / 2708 vertices
METRICS AT BATCH-9
Insert time (avg) = 0.691839
GUP (avg) = 0.0573472
BWUP (avg, includes merge) = 0.618374, average walk affected = 26253.1, sampled vertices = 17855585
WUP (avg)   = 0.364753  (sampling= 0.30611, inserting= 0.0586426)
MAV (avg)   = 0.0570162 MAV (min) = 0.0459552   MAV (max) = 0.299467
Merge (avg,10 times) = 0.233503 Merge (min) = 0.195651  Merge (max) = 1.13015

batch-10 and batch_seed-10
(merge) -- For batch-10 we are touching 2708 / 2708 vertices
(insert) -- For batch-10 we are touching 2569 / 2708 vertices
METRICS AT BATCH-10
Insert time (avg) = 0.682258
GUP (avg) = 0.0582243
BWUP (avg, includes merge) = 0.607583, average walk affected = 26275.2, sampled vertices = 19871658
WUP (avg)   = 0.358441  (sampling= 0.300859, inserting= 0.057582)
MAV (avg)   = 0.0578846 MAV (min) = 0.0459552   MAV (max) = 0.299467
Merge (avg,10 times) = 0.229404 Merge (min) = 0.195651  Merge (max) = 1.16389

(Last Merge) -- we are touching 2708 / 2708 vertices
Last merge (with MIN-MAX Ranges) time: 0.06975889

METRICS FOR ALL BATCHES
Insert time (avg) = 0.682258
GUP (avg) = 0.058224
BWUP (avg, includes merge) = 0.607583, average walk affected = 26275.200000, sampled vertices = 19871658
WUP (avg)   = 0.358441  (sampling= 0.300859, inserting= 0.057582)
MAV (avg)   = 0.057885  MAV (min) = 0.045955    MAV (max) = 0.299467
Merge (avg,10.000000 times) = 0.229404  Merge (min) = 0.195651  Merge (max) = 1.163887
Average walk insert latency = { 0.000022 0.000024 0.000046 0.000050 0.000069 0.000073 0.000092 0.000097 0.000114 0.000116 }
(1) throughput: 43245.41407922
(2) average latency: 0.00007036
```



## Vertex Classification with Wharf's walks

The vertex classification experiment validates whether the random walks produced and maintained by Wharf are capable of producing graph embeddings of high quality that can be effectively used in a downstream vertex classification task. The experiment can be run using a script called `vertex-classification.sh` located in `experiments/scripts/`. Note that for the vertex classification experiment the user should checkout to the `vertex_classification_exp` branch and use that code. 

The easiest way we can reproduce the vertex classification experimental results in Figure 11a of [1] is to open, select the appropriate parameters, and run the `vertex-classification.sh` script located in `experiments/scripts/` of the `vertex_classification_exp` branch. However, the same results can be obtained after running the executable `vertex-classification` located in `build/experiments/` for the corresponding datasets and set of parameters. 
Note that because the vertex classification is a supervised task, i.e., requires the true labels of the vertices to be classified, we should also set the proper label file of the dataset used in the script `experiments/bin/vertex-classification.py` (line 11). Bellow we demonstrate the output on `cora` dataset for a example set of parameters that we run on our laptop:

```
$ ./vertex-classification.sh  # run the script from the script/ folder
cora-graph - Running vertex classification with 8 threads.

Incremental Learning
Walks per vertex: 10
Walk length: 80
Vector dimension: 128
Learning strategy: ONLINE
Walking model: DEEPWALK
Sampler strategy: WEIGHT
Vertices: 2708 Edges: 10556
Learning initial embeddings
Initializing model... done
Training incremental SGNS
.. done (inf=27080/0 sent/sec)
Initializing model... done
...
Initializing model... done
Training incremental SGNS
. done (736.000000=14720/20 sent/sec)
IncrementalTimer: Total: 766.5944
Static Learning
Walks per vertex: 10
Walk length: 80
Vector dimension: 128
Learning strategy: BATCH
Walking model: DEEPWALK
Sampler strategy: WEIGHT
Vertices: 2708 Edges: 10556
Learning initial embeddings
Initializing model... done
Training incremental SGNS
.. done (inf=27080/0 sent/sec)
Initializing model... done
Training incremental SGNS
.. done (4513.333333=27080/6 sent/sec)
...
Initializing model... done
Training incremental SGNS
.. done (615.454545=27080/44 sent/sec)
Initializing model... done
Training incremental SGNS
.. done (530.980392=27080/51 sent/sec)
Initializing model... done
Training incremental SGNS
.. done (530.980392=27080/51 sent/sec)
StaticTimer: Total: 875.8127
```

# Wharf++

There are a few interesting directions to consider for extending Wharf further. If you are interested in helping us make the ship really "sail off into deep, cool waters and away from the wharf" drop us an email. ;-)

# References

[1]: Space-Efficient Random Walks on Streaming Graphs, Serafeim Papadias, Zoi Kaoudi, Jorge Quiane, and Volker Markl

[2]: Low-Latency Graph Streaming Using Compressed Purely-Functional Trees, Laxman Dhulipala, Guy Blelloch, and Julian Shun

# Cite Us

To cite our paper, add the following BibTeX snippet of the extended version to your bibliography:

```
@article{DBLP:journals/pvldb/PapadiasKQM22,
  author       = {Serafeim Papadias and
                  Zoi Kaoudi and
                  Jorge{-}Arnulfo Quian{\'{e}}{-}Ruiz and
                  Volker Markl},
  title        = {Space-Efficient Random Walks on Streaming Graphs},
  journal      = {Proc. {VLDB} Endow.},
  volume       = {16},
  number       = {2},
  pages        = {356--368},
  year         = {2022},
  url          = {https://www.vldb.org/pvldb/vol16/p356-papadias.pdf},
  timestamp    = {Mon, 26 Jun 2023 20:53:23 +0200},
  biburl       = {https://dblp.org/rec/journals/pvldb/PapadiasKQM22.bib},
  bibsource    = {dblp computer science bibliography, https://dblp.org}
}
```
