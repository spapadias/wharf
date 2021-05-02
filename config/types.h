#ifndef DYNAMIC_GRAPH_REPRESENTATION_LEARNING_WITH_METROPOLIS_HASTINGS_TYPES_H
#define DYNAMIC_GRAPH_REPRESENTATION_LEARNING_WITH_METROPOLIS_HASTINGS_TYPES_H

namespace types
{
    // graph vertex type
    using Vertex = uintV;

    // class that stored edges in a compressed purely functional balanced binary tree
    using CompressedEdges = tree_plus::treeplus;

    // 
    using CompressedEdgesLists = tree_plus::edge_list;
}

#endif // DYNAMIC_GRAPH_REPRESENTATION_LEARNING_WITH_METROPOLIS_HASTINGS_TYPES_H
