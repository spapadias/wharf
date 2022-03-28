#ifndef VERTEX_H
#define VERTEX_H

#include <compressed_walks.h>
#include <metropolis_hastings_sampler.h>

namespace dynamic_graph_representation_learning_with_metropolis_hastings
{
    // SamplerManager = parallel hash map that stores all samplers for one graph vertex
    using SamplerManager = libcuckoo::cuckoohash_map<types::Vertex, dygrl::MetropolisHastingsSampler>;

    /**
     * @brief VertexEntry represents a structure that contains the vertex data - compressed edges, compressed walks, and sampler manager.
     */
    struct VertexEntry
    {
        types::CompressedEdges compressed_edges;
        std::vector<dygrl::CompressedWalks> compressed_walks;
        dygrl::SamplerManager* sampler_manager;

        /**
         * @brief VertexEntry default constructor.
         */
        VertexEntry()
        {
            this->compressed_edges = types::CompressedEdges();
            this->compressed_walks = vector<dygrl::CompressedWalks>();
            this->sampler_manager  = nullptr;
        }

        /**
         * @brief VertexEntry constructor.
         *
         * @param compressed_edges - compressed tree of edges
         * @param compressed_walks - compressed tree of walks
         * @param sampler_manager  - manager of MH samplers
         */
        VertexEntry (const types::CompressedEdges& compressed_edges, const std::vector<dygrl::CompressedWalks>& compressed_walks, dygrl::SamplerManager* sampler_manager)
            : compressed_edges(compressed_edges), compressed_walks(compressed_walks), sampler_manager(sampler_manager)
			{ };
    };

    /**
     * @brief Graph vertex structure.
     */
    struct Vertex
    {
        using key_t = types::Vertex;              // key: vertex id
        using val_t = VertexEntry;                // value: compressed edges, compressed walks and metropolis hastings samplers
        using aug_t = types::Degree;              // augmentation: vertex degree

        using entry_t = std::pair<key_t, val_t>;  // vertex - <vertex id, {compressed-edges, compressed-walks, MH samplers}>

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
            // copy compressed edges
            auto ce_plus = lists::copy_node(entry.second.compressed_edges.plus);         // ce plus part; bumps ref-cnt
            auto ce_root = edge_plus::Tree_GC::inc(entry.second.compressed_edges.root);  // ce root part; bumps ref-cnt

			std::vector<CompressedWalks> walk_trees_cpy;
			for (const auto& cw : entry.second.compressed_walks)
			{
				auto cw_plus = lists::copy_node(cw.plus);         // cw plus part; bumps ref-cnt
				auto cw_root = walk_plus::Tree_GC::inc(cw.root);  // cw root part; bumps ref-cnt

				// Copy the (min, max) range as well
				auto vnext_min = cw.vnext_min;
				auto vnext_max = cw.vnext_max;
				auto batch_num = cw.created_at_batch;

				walk_trees_cpy.push_back(dygrl::CompressedWalks(cw_plus, cw_root, vnext_min, vnext_max, batch_num));
			}

            // copy sampler manager
            auto sampler = new SamplerManager(entry.second.sampler_manager->size());
            for(auto& table_entry : entry.second.sampler_manager->lock_table())
            {
                sampler->insert(table_entry.first, table_entry.second);
            }

            return std::make_pair(entry.first, VertexEntry(types::CompressedEdges(ce_plus, ce_root), walk_trees_cpy, sampler));
        }

        // delete an entry
        static void del(entry_t& entry)
        {
            // delete compressed edges
            if (entry.second.compressed_edges.plus)
            {
                lists::deallocate(entry.second.compressed_edges.plus);
                entry.second.compressed_edges.plus = nullptr;
            }

            if (entry.second.compressed_edges.root)
            {
                auto T = edge_plus::edge_list();

                T.root = entry.second.compressed_edges.root;
                entry.second.compressed_edges.root = nullptr;
            }

	        for (auto& cw : entry.second.compressed_walks)
	        {
		        if (cw.plus)
		        {
			        lists::deallocate(cw.plus);
			        cw.plus = nullptr;
		        }

		        if (cw.root)
		        {
			        auto T = walk_plus::edge_list();

			        T.root = cw.root;
			        cw.root = nullptr;
		        }
			}

            // delete sampler manager
            entry.second.sampler_manager->clear();
            delete entry.second.sampler_manager;
        }
    };
}

#endif
