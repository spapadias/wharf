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
             * @brief CompressedWalks default constructor.
             */
            CompressedWalks() : tree_plus::treeplus() {};

            /**
             * @brief CompressedWalks constructor.
             *
             * @param plus - tree plus
             * @param root - tree root
             */
            CompressedWalks(tree_plus::AT* plus, tree_plus::treeplus::Node* root)
                : tree_plus::treeplus(plus, root) {};

            /**
            * @brief CompressedWalks constructor.
            *
            * @tparam Sequence
            *
            * @param sequence - sequence of walk parts
            * @param source   - graph vertex
            * @param flag
            */
            template<class Sequence>
            CompressedWalks(const Sequence &sequence, types::Vertex source, pbbs::flags flag = pbbs::no_flag)
                : tree_plus::treeplus(sequence, source, flag) {};

            /**
            * @brief Finds the next vertex in the walk given walk id and position
            *
            * @param walk_id  - unique walk id
            * @param position - position in the walk
            * @param source   - current walk vertex
            *
            * @return - next vertex in the walk
            */
            types::Vertex find_next(types::WalkID walk_id, types::Position position, types::Vertex source)
            {
                types::Vertex next_vertex = -1;

                bool result = this->iter_elms_cond(source, [&](auto value)
                {
                    auto pair = pairings::Szudzik<types::Vertex>::unpair(value);

                    auto this_walk_id  = pair.first / config::walk_length;
                    auto this_position = pair.first - (this_walk_id * config::walk_length);
                    next_vertex        = pair.second;

                    if (this_walk_id == walk_id && this_position == position)
                        return true;
                    else
                        return false;
                });

                #ifdef MALIN_DEBUG
                    if (!result || next_vertex == -1)
                    {
                        std::cerr << "Dock debug error! CompressedWalks::FindNext::walk_id = "
                                  << walk_id
                                  << ", position = "
                                  << (int) position
                                  << ", vertex = "
                                  << source
                                  << std::endl;

                        std::exit(1);
                    }
                #endif

                return next_vertex;
            }
    };
}

#endif // DYNAMIC_GRAPH_REPRESENTATION_LEARNING_WITH_METROPOLIS_HASTINGS_COMPRESSED_WALKS_H
