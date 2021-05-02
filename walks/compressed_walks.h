#ifndef DYNAMIC_GRAPH_REPRESENTATION_LEARNING_WITH_METROPOLIS_HASTINGS_COMPRESSED_WALKS_H
#define DYNAMIC_GRAPH_REPRESENTATION_LEARNING_WITH_METROPOLIS_HASTINGS_COMPRESSED_WALKS_H

#include <graph/tree_plus/tree_plus.h>

namespace dynamic_graph_representation_learning_with_metropolis_hastings
{
    struct edge_entry
    {
        using key_t = uintV;
        using val_t = tree_plus::AT*;
        using aug_t = unsigned char;

        using entry_t = std::pair<key_t,val_t>;

        // Key x Key -> Key
        static bool comp(const key_t& a, const key_t& b)
        {
            return a < b;
        }

        // Key x Value -> Augmentation
        static aug_t from_entry(const key_t& k, const val_t& v)
        {
            return 0;
        }

        // Augmentation x Augmentation -> Augmentation
        static aug_t combine(const aug_t& a, const aug_t& b)
        {
            return a + b;
        }

        // Empty -> Augmentation (default augmentation)
        static aug_t get_empty()
        {
            return 0;
        }

        // Copy existing entry
        static entry_t copy_entry(const entry_t& e)
        {
            // TODO: Instead of copying, bump a ref-ct (note that copy_node and
            // deallocate can implement these semantics internally)
            return make_pair(e.first, lists::copy_node(e.second));
        }

        static void del(entry_t& e)
        {
            if (e.second)
            {
                // TODO: Should decrement ref-ct, free if only owner
                lists::deallocate(e.second);
            }
        }
    };

    using edge_list = aug_map<edge_entry>;
    using CompressedWalksLists = edge_list;

    /**
    * CompressedWalks represents a class that stores graph random walks in a compressed format.
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
    };
}

#endif // DYNAMIC_GRAPH_REPRESENTATION_LEARNING_WITH_METROPOLIS_HASTINGS_COMPRESSED_WALKS_H
