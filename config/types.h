#ifndef DYNAMIC_GRAPH_REPRESENTATION_LEARNING_WITH_METROPOLIS_HASTINGS_TYPES_H
#define DYNAMIC_GRAPH_REPRESENTATION_LEARNING_WITH_METROPOLIS_HASTINGS_TYPES_H

namespace types
{
    // Vertex = vertex data type
    using Vertex               = uintV;

    // Degree = vertex degree
    using Degree               = uint16_t;

    // WalkID = unique walk id
    using WalkID               = uintV;

    // Position = the position of a vertex in the walk
    using Position             = uint8_t;

    // PairedTriplet = a triplet <walk id, position, next vertex> after encoded with the pairing function
    using PairedTriplet        = uintV;

    // State = the state is defined as a pair of two numbers, where the first represents the current vertex,
    // and the second contains the extra information
    // (e.g DeepWalk = current vertex, Node2Vec = previously visited vertex by the walker)
    using State                = std::pair<Vertex, Vertex>;

    // CompressedEdges = data structure (augmented parallel binary tree) that stores compressed edges & compressed walks
    using CompressedEdges      = tree_plus::treeplus;

    // CompressedTreesLists = data structure (augmented parallel binary tree) that stores compressed edges & compressed walks
    using CompressedTreesLists = tree_plus::edge_list;
}

#endif // DYNAMIC_GRAPH_REPRESENTATION_LEARNING_WITH_METROPOLIS_HASTINGS_TYPES_H
