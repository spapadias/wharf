#ifndef DYNAMIC_GRAPH_REPRESENTATION_LEARNING_WITH_METROPOLIS_HASTINGS_TYPES_H
#define DYNAMIC_GRAPH_REPRESENTATION_LEARNING_WITH_METROPOLIS_HASTINGS_TYPES_H

namespace types
{
    using Vertex               = uintV;

    using WalkID               = uintV;

    using Position             = unsigned char;

    using PairedTriplet        = uintV;

    using CompressedEdges      = tree_plus::treeplus;

    using CompressedTreesLists = tree_plus::edge_list;
}

#endif // DYNAMIC_GRAPH_REPRESENTATION_LEARNING_WITH_METROPOLIS_HASTINGS_TYPES_H
