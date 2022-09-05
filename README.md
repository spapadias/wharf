# Wharf: Space-Efficient Random Walks on Streaming Graphs

This repository provides the code, unit tests, scripts, and instructions for reproducing the experiments in our paper, Space-Efficient Random Walks on Streaming Graphs [1]. Our paper introduces Wharf, a main-memory streaming random walk system capable of updating, indexing, and compressing a maintained random walk corpus efficiently over streaming (i.e., highly dynamic) graphs that are evolving.

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
- Unit Tests
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

Note that currently in our repo we only provide Cora dataset (along with its label for the vertex classification experiment) for the sake of space.

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

The easiest way we can reproduce the memory footprint experimental results in Figure 8 of [1] is to open, select the appropriate parameters, and run the `memory-footprint.sh` script located in `experiments/scripts/`. However, the same results can be obtained after running the executable `memory-footprint` located in `build/experiments/` for the corresponding datasets and set of parameters. Bellow we demonstrate the output on `Cora` dataset for a example set of parameters:

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

```
Graph: cora-graph 
Running experiment with 8 threads
Parallel Merge and WU
Walks per vertex: 10
Walk length: 80
Walking model: DEEPWALK
Sampler strategy: WEIGHT
n = 2708 m = 10556
GENERATE 27080 walks.

Total time to generate the walk corpus from scratch: 0.486084
Batch size = 7000 | batch-0 and batch_seed-0
6778 
^^^ merge at batch-1
all_to_delete size: 2708

(merge) -- For batch-1 we are touching 2708 / 2708 vertices

(insert) -- For batch-1 we are touching 2641 / 2708 vertices
METRICS AT BATCH-1
Insert time (avg) = 0.791444
GUP (avg) = 0.151061
BWUP (avg, includes merge) = 0.624453, average walk affected = 26430, sampled vertices = 2084937
WUP (avg)   = 0.406891  (sampling= 0.350641, inserting= 0.0562501)
MAV (avg)   = 0.139444  MAV (min) = 0.139444    MAV (max) = 0.139444
Merge (avg,10 times) = 0.206624 Merge (min) = 0.206624  Merge (max) = 0.206624
batch-1 and batch_seed-1
6784 
^^^ merge at batch-2
all_to_delete size: 2708

(merge) -- For batch-2 we are touching 2708 / 2708 vertices

(insert) -- For batch-2 we are touching 2638 / 2708 vertices
METRICS AT BATCH-2
Insert time (avg) = 0.834294
GUP (avg) = 0.13035
BWUP (avg, includes merge) = 0.690628, average walk affected = 26480, sampled vertices = 4182309
WUP (avg)   = 0.433448  (sampling= 0.367608, inserting= 0.06584)
MAV (avg)   = 0.12376   MAV (min) = 0.108076    MAV (max) = 0.139444
Merge (avg,10 times) = 0.241508 Merge (min) = 0.206624  Merge (max) = 0.276393
batch-2 and batch_seed-2
6776 
^^^ merge at batch-3
all_to_delete size: 2708

(merge) -- For batch-3 we are touching 2708 / 2708 vertices

(insert) -- For batch-3 we are touching 2656 / 2708 vertices
METRICS AT BATCH-3
Insert time (avg) = 0.869039
GUP (avg) = 0.154528
BWUP (avg, includes merge) = 0.695252, average walk affected = 26526.7, sampled vertices = 6287888
WUP (avg)   = 0.465676  (sampling= 0.409614, inserting= 0.056062)
MAV (avg)   = 0.145377  MAV (min) = 0.108076    MAV (max) = 0.328056
Merge (avg,10 times) = 0.205796 Merge (min) = 0.206624  Merge (max) = 0.340996
batch-3 and batch_seed-3
6762 
^^^ merge at batch-4
all_to_delete size: 2708

(merge) -- For batch-4 we are touching 2708 / 2708 vertices

(insert) -- For batch-4 we are touching 2639 / 2708 vertices
METRICS AT BATCH-4
Insert time (avg) = 0.798022
GUP (avg) = 0.138188
BWUP (avg, includes merge) = 0.643003, average walk affected = 26546, sampled vertices = 8390082
WUP (avg)   = 0.419004  (sampling= 0.364619, inserting= 0.0543845)
MAV (avg)   = 0.130901  MAV (min) = 0.108076    MAV (max) = 0.328056
Merge (avg,10 times) = 0.203914 Merge (min) = 0.206624  Merge (max) = 0.474662
batch-4 and batch_seed-4
6724 
^^^ merge at batch-5
all_to_delete size: 2708

(merge) -- For batch-5 we are touching 2708 / 2708 vertices

(insert) -- For batch-5 we are touching 2628 / 2708 vertices
METRICS AT BATCH-5
Insert time (avg) = 0.746375
GUP (avg) = 0.127611
BWUP (avg, includes merge) = 0.603363, average walk affected = 26557.6, sampled vertices = 10487758
WUP (avg)   = 0.39056   (sampling= 0.33812, inserting= 0.0524396)
MAV (avg)   = 0.122396  MAV (min) = 0.108076    MAV (max) = 0.416434
Merge (avg,10 times) = 0.193003 Merge (min) = 0.206624  Merge (max) = 0.490354
batch-5 and batch_seed-5
6726 
^^^ merge at batch-6
all_to_delete size: 2708

(merge) -- For batch-6 we are touching 2708 / 2708 vertices

(insert) -- For batch-6 we are touching 2629 / 2708 vertices
METRICS AT BATCH-6
Insert time (avg) = 0.718858
GUP (avg) = 0.119476
BWUP (avg, includes merge) = 0.585082, average walk affected = 26565.3, sampled vertices = 12585889
WUP (avg)   = 0.369326  (sampling= 0.31646, inserting= 0.0528667)
MAV (avg)   = 0.115124  MAV (min) = 0.108076    MAV (max) = 0.416434
Merge (avg,10 times) = 0.196194 Merge (min) = 0.206624  Merge (max) = 0.68681
batch-6 and batch_seed-6
6760 
^^^ merge at batch-7
all_to_delete size: 2708

(merge) -- For batch-7 we are touching 2708 / 2708 vertices

(insert) -- For batch-7 we are touching 2641 / 2708 vertices
METRICS AT BATCH-7
Insert time (avg) = 0.714115
GUP (avg) = 0.117085
BWUP (avg, includes merge) = 0.583147, average walk affected = 26570.9, sampled vertices = 14685039
WUP (avg)   = 0.361723  (sampling= 0.307634, inserting= 0.0540893)
MAV (avg)   = 0.112955  MAV (min) = 0.108076    MAV (max) = 0.516372
Merge (avg,10 times) = 0.203085 Merge (min) = 0.206624  Merge (max) = 0.734784
batch-7 and batch_seed-7
6754 
^^^ merge at batch-8
all_to_delete size: 2708

(merge) -- For batch-8 we are touching 2708 / 2708 vertices

(insert) -- For batch-8 we are touching 2631 / 2708 vertices
METRICS AT BATCH-8
Insert time (avg) = 0.694755
GUP (avg) = 0.112607
BWUP (avg, includes merge) = 0.567896, average walk affected = 26575, sampled vertices = 16783322
WUP (avg)   = 0.352945  (sampling= 0.299785, inserting= 0.0531594)
MAV (avg)   = 0.108348  MAV (min) = 0.108076    MAV (max) = 0.516372
Merge (avg,10 times) = 0.197533 Merge (min) = 0.206624  Merge (max) = 0.845478
batch-8 and batch_seed-8
6742 
^^^ merge at batch-9
all_to_delete size: 2708

(merge) -- For batch-9 we are touching 2708 / 2708 vertices

(insert) -- For batch-9 we are touching 2637 / 2708 vertices
METRICS AT BATCH-9
Insert time (avg) = 0.673956
GUP (avg) = 0.10743
BWUP (avg, includes merge) = 0.552404, average walk affected = 26578.2, sampled vertices = 18882937
WUP (avg)   = 0.343465  (sampling= 0.291623, inserting= 0.051841)
MAV (avg)   = 0.103379  MAV (min) = 0.108076    MAV (max) = 0.579997
Merge (avg,10 times) = 0.191671 Merge (min) = 0.206624  Merge (max) = 0.879559
batch-9 and batch_seed-9
6750 
^^^ merge at batch-10
all_to_delete size: 2708

(merge) -- For batch-10 we are touching 2708 / 2708 vertices

(insert) -- For batch-10 we are touching 2642 / 2708 vertices
METRICS AT BATCH-10
Insert time (avg) = 0.684484
GUP (avg) = 0.10867
BWUP (avg, includes merge) = 0.561878, average walk affected = 26582.4, sampled vertices = 20985342
WUP (avg)   = 0.344176  (sampling= 0.291681, inserting= 0.0524954)
MAV (avg)   = 0.105649  MAV (min) = 0.108076    MAV (max) = 0.579997
Merge (avg,10 times) = 0.200664 Merge (min) = 0.206624  Merge (max) = 1.12708

METRICS FOR ALL BATCHES
Insert time (avg) = 0.684484
GUP (avg) = 0.108670
BWUP (avg, includes merge) = 0.561878, average walk affected = 26582.400000, sampled vertices = 20985342
WUP (avg)   = 0.344176  (sampling= 0.291681, inserting= 0.052495)
MAV (avg)   = 0.105649  MAV (min) = 0.108076    MAV (max) = 0.579997
Merge (avg,10.000000 times) = 0.200664  Merge (min) = 0.206624  Merge (max) = 1.127080
Average walk insert latency = { 0.000024 0.000029 0.000050 0.000047 0.000067 0.000065 0.000088 0.000083 0.000104 0.000107 }
(1) throughput: 47309.92113572
(2) average latency: 0.00006627
FindPrev vertex in node2vec: 0.00000000
(11) SIMPLE SEARCH throughput: 47309.92113572
FindPrev SIMPLE in node2vec: 0.00000000
(111) RANGE SEARCH throughput: 47309.92113572
FindPrev vertex in node2vec: 0.00000000
all_to_delete size: 2708

(Last Merge) -- we are touching 2708 / 2708 vertices
Last merge (with MIN-MAX Ranges) time: 0.08695889

```



## Vertex Classification with Wharf's walks

The vertex classification experiment validates whether the random walks produced and maintained by Wharf are capable of producing graph embeddings of high quality that can be effectively used in a downstream vertex classification task. The experiment can be run using a script called `vertex-classification.sh` located in `experiments/scripts/`. Note that for the vertex classification experiment the user should checkout to the `vertex_classification_exp` branch and use that code. 

The easiest way we can reproduce the vertex classification experimental results in Figure 11a of [1] is to open, select the appropriate parameters, and run the `vertex-classification.sh` script located in `experiments/scripts/` of the `vertex_classification_exp` branch. However, the same results can be obtained after running the executable `vertex-classification` located in `build/experiments/` for the corresponding datasets and set of parameters. 
Note that because the vertex classification is a supervised task, i.e., requires the true labels of the vertices to be classified, we should also set the proper label file of the dataset used in the script `experiments/bin/vertex-classification.py` (line 11). Bellow we demonstrate the output on `Cora` dataset for a example set of parameters that we run on our laptop:

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

# Unit Tests

**Currently the units test are not included. Go to the main branch and fetch the clean unit tests in the final stage**

Note that the repo contains unit tests under the `tests` folder, which we wrote during the development of Wharf. The ones interested for extending Wharf may consider using our tests or creating their one in a similar manner. Currently, we are using the [latest version](https://github.com/google/googletest/releases/tag/release-1.12.1) of [GoogleTest](https://github.com/google/googletest).

# Wharf++

There are a few interesting directions to consider for extending Wharf further. If you are interested in helping us make the ship really "sail off into cool waters from the wharf" drop us an email. ;-)

# References

[1]: Space-Efficient Random Walks on Streaming Graphs, Serafeim Papadias, Zoi Kaoudi, Jorge Quiane, and Volker Markl

[2]: Low-Latency Graph Streaming Using Compressed Purely-Functional Trees, Laxman Dhulipala, Guy Blelloch, and Julian Shun
