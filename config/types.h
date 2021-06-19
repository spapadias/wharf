#ifndef DYNAMIC_GRAPH_REPRESENTATION_LEARNING_WITH_METROPOLIS_HASTINGS_TYPES_H
#define DYNAMIC_GRAPH_REPRESENTATION_LEARNING_WITH_METROPOLIS_HASTINGS_TYPES_H

namespace types
{
    // Vertex = graph vertex data type
    using Vertex               = uintV;

    // Degree = graph vertex degree data type
    using Degree               = uintE;

    // Neighbors = graph vertex neighbors
    using Neighbors            = std::tuple<Vertex*, Degree, bool>;

    // WalkID = unique walk id
    using WalkID               = uint32_t;

    // Position = the position of a vertex in the walk
    using Position             = uint8_t;

    // State = the state is defined as a pair of two numbers,
    // where the first represents the current vertex and the second contains the extra information
    // (e.g DeepWalk = current vertex, Node2Vec = previously visited vertex by the walker)
    using State                = std::pair<Vertex, Vertex>;

    // CompressedEdges = structure (augmented parallel binary tree) that stores compressed edges
    using CompressedEdges      = tree_plus::treeplus;
    using CompressedTreesLists = tree_plus::edge_list;

    // RandomWalkModelType = different walking models
    enum RandomWalkModelType   { DEEPWALK, NODE2VEC };

    // StartMode = edge sampler initialization strategy
    enum SamplerInitStartegy   { RANDOM, BURNIN, WEIGHT };

    // Global Map of Changes (MoC) = contains starting positions to crop the walk
    using MapOfChanges         = libcuckoo::cuckoohash_map<WalkID, std::pair<Position, Vertex>>;

    // ChangeAccumulator = accumulator of changes for walk updates
    using ChangeAccumulator    = libcuckoo::cuckoohash_map<Vertex, std::vector<std::pair<uintV, types::Vertex>>>;
}

#endif // DYNAMIC_GRAPH_REPRESENTATION_LEARNING_WITH_METROPOLIS_HASTINGS_TYPES_H
