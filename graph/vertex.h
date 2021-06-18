#ifndef DYNAMIC_GRAPH_REPRESENTATION_LEARNING_WITH_METROPOLIS_HASTINGS_VERTEX_H
#define DYNAMIC_GRAPH_REPRESENTATION_LEARNING_WITH_METROPOLIS_HASTINGS_VERTEX_H

#include <inverted_index.h>
#include <metropolis_hastings_sampler.h>

namespace dynamic_graph_representation_learning_with_metropolis_hastings
{
    // SamplerManager = parallel hash map that stores all samplers for one graph vertex
    using SamplerManager = libcuckoo::cuckoohash_map<types::Vertex, dygrl::MetropolisHastingsSampler>;

    /**
     * @brief VertexEntry represents a structure that contains the vertex data - compressed edges, inverted index, and sampler manager.
     */
    struct VertexEntry
    {
        types::CompressedEdges compressed_edges;
        dygrl::InvertedIndex   inverted_index;
        dygrl::SamplerManager* sampler_manager;

        /**
         * @brief VertexEntry default constructor.
         */
        VertexEntry()
        {
            this->compressed_edges = types::CompressedEdges();
            this->inverted_index   = dygrl::InvertedIndex();
            this->sampler_manager  = nullptr;
        }

        /**
         * @brief VertexEntry constructor.
         *
         * @param compressed_edges - compressed tree of edges
         * @param inverted_index   - walks inverted index
         * @param sampler_manager  - manager of MH samplers
         */
        VertexEntry(const types::CompressedEdges& compressed_edges, const dygrl::InvertedIndex& inverted_index, dygrl::SamplerManager* sampler_manager)
            : compressed_edges(compressed_edges), inverted_index(inverted_index), sampler_manager(sampler_manager) {};
    };

    /**
     * @brief Graph vertex structure.
     */
    struct Vertex
    {
        using key_t = types::Vertex;   // key: vertex id
        using val_t = VertexEntry;     // value: compressed edges, walks inverted index and metropolis hastings samplers
        using aug_t = types::Degree;   // augmentation: vertex degree

        using entry_t = std::pair<key_t, val_t>;  // vertex - <vertex id, {compressed-edges, walks index, MH samplers}>

        // key x key -> key
        static bool comp(const key_t& keyX, const key_t& keyY)
        {
            return keyX < keyY;
        }

        // key x value -> augmentation
        static aug_t from_entry(const key_t& key, const val_t& value)
        {
            return value.compressed_edges.size();
        }

        // augmentation x augmentation -> augmentation
        static aug_t combine(const aug_t& augX, const aug_t& augY)
        {
            return augX + augY;
        }

        // empty -> augmentation (default augmentation)
        static aug_t get_empty()
        {
            return 0;
        }

        // copy existing entry
        static entry_t copy_entry(const entry_t& entry)
        {
            // copy edges
            auto ce_plus = lists::copy_node(entry.second.compressed_edges.plus);         // ce plus part; bumps ref-cnt
            auto ce_root = tree_plus::Tree_GC::inc(entry.second.compressed_edges.root);  // ce root part; bumps ref-cnt

            // copy samplers
            auto sampler = new SamplerManager(entry.second.sampler_manager->size());
            for(auto& table_entry : entry.second.sampler_manager->lock_table())
            {
                sampler->insert(table_entry.first, table_entry.second);
            }

            return std::make_pair(entry.first,
                VertexEntry
                (
                    types::CompressedEdges(ce_plus, ce_root),
                    dygrl::InvertedIndex(entry.second.inverted_index),          // copy inverted index
                    sampler
                )
            );
        }

        // delete an entry
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

            entry.second.inverted_index.~InvertedIndex();

            delete entry.second.sampler_manager;
        }
    };
}

#endif // DYNAMIC_GRAPH_REPRESENTATION_LEARNING_WITH_METROPOLIS_HASTINGS_VERTEX_H
