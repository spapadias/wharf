#ifndef DYNAMIC_GRAPH_REPRESENTATION_LEARNING_WITH_METROPOLIS_HASTINGS_COMPRESSED_WALKS_H
#define DYNAMIC_GRAPH_REPRESENTATION_LEARNING_WITH_METROPOLIS_HASTINGS_COMPRESSED_WALKS_H

#include <graph/tree_plus/tree_plus.h>

namespace dynamic_graph_representation_learning_with_metropolis_hastings
{
    /**
    * CompressedWalks is a structure that stores graph random walks in a compressed format.
    * It essentially represents a compressed purely functional tree (C-Tree) that achieves compression
    * based on differential coding.
    */
    class CompressedWalks : public tree_plus::treeplus
    {
        public:
            /**
             * CompressedWalks default constructor.
             */
            CompressedWalks() : tree_plus::treeplus() {};

            /**
             * CompressedWalks constructor.
             *
             * @param plus
             * @param root
             */
            CompressedWalks(tree_plus::AT* plus, tree_plus::treeplus::Node* root)
                : tree_plus::treeplus(plus, root) {};

            /**
            * CompressedWalks constructor.
            *
            * @tparam Sequence
            *
            * @param sequence
            * @param source
            * @param flag
            */
            template<class Sequence>
            CompressedWalks(const Sequence &sequence, types::Vertex source, pbbs::flags flag = pbbs::no_flag)
                : tree_plus::treeplus(sequence, source, flag) {};
    };
}

#endif // DYNAMIC_GRAPH_REPRESENTATION_LEARNING_WITH_METROPOLIS_HASTINGS_COMPRESSED_WALKS_H
