#ifndef DYNAMIC_GRAPH_REPRESENTATION_LEARNING_WITH_METROPOLIS_HASTINGS_VERTEX_H
#define DYNAMIC_GRAPH_REPRESENTATION_LEARNING_WITH_METROPOLIS_HASTINGS_VERTEX_H

namespace dynamic_graph_representation_learning_with_metropolis_hastings
{
    struct VertexEntry
    {
        types::CompressedEdges compressed_edges;

        /**
         * VertexEntry default constructor
         */
        VertexEntry()
        {
            this->compressed_edges = types::CompressedEdges();
        }

        /**
         * VertexEntry constructor
         *
         * @param compressed_edges
         */
        VertexEntry(const types::CompressedEdges& compressed_edges)
            : compressed_edges(compressed_edges) {};
    };

    struct Vertex
    {
        using key_t = uintV;        // key: vertex id
        using val_t = VertexEntry;  // value: compressed edges, compressed walks and metropolis hastings samplers
        using aug_t = uintE;        // augmentation: vertex degree

        using entry_t = std::pair<key_t, val_t>;  // vertex <vertex id, {compressed-edges, compressed-walks, MH samplers}>

        // Key x Key -> Key
        static bool comp(const key_t& keyX, const key_t& keyY)
        {
            return keyX < keyY;
        }

        // Key x Value -> Augmentation
        static aug_t from_entry(const key_t& key, const val_t& value)
        {
            return value.compressed_edges.size();
        }

        // Augmentation x Augmentation -> Augmentation
        static aug_t combine(const aug_t& augX, const aug_t& augY)
        {
            return augX + augY;
        }

        // Empty -> Augmentation (default augmentation)
        static aug_t get_empty()
        {
            return 0;
        }

        // Copy existing entry
        static entry_t copy_entry(const entry_t& entry)
        {
            auto ce_plus = lists::copy_node(entry.second.compressed_edges.plus);         // ce plus part; bumps ref-cnt
            auto ce_root = tree_plus::Tree_GC::inc(entry.second.compressed_edges.root);  // ce root part; bumps ref-cnt

//            auto cw_plus = lists::copy_node(entry.second.compressed_walks.plus);         // cw plus part; bumps ref-cnt
//            auto cw_root = tree_plus::Tree_GC::inc(entry.second.compressed_walks.root);  // cw root part; bumps ref-cnt

            return std::make_pair(entry.first, VertexEntry(types::CompressedEdges(ce_plus, ce_root)));
        }

        // Delete an entry
        static void del(entry_t& entry)
        {
            if (entry.second.compressed_edges.plus)
            {
                lists::deallocate(entry.second.compressed_edges.plus);
                entry.second.compressed_edges.plus = nullptr;
            }

            if (entry.second.compressed_edges.root)
            {
                auto T = tree_plus::edge_list();

                T.root = entry.second.compressed_edges.root;
                entry.second.compressed_edges.root = nullptr;
            }

//            if (entry.second.compressed_walks.plus)
//            {
//                lists::deallocate(entry.second.compressed_walks.plus);
//                entry.second.compressed_walks.plus = nullptr;
//            }
//
//            if (entry.second.compressed_walks.root)
//            {
//                auto T = tree_plus::edge_list();
//
//                T.root = entry.second.compressed_walks.root;
//                entry.second.compressed_walks.root = nullptr;
//            }
        }
    };
}

#endif // DYNAMIC_GRAPH_REPRESENTATION_LEARNING_WITH_METROPOLIS_HASTINGS_VERTEX_H
