#ifndef DYNAMIC_GRAPH_REPRESENTATION_LEARNING_WITH_METROPOLIS_HASTINGS_DOCK_H
#define DYNAMIC_GRAPH_REPRESENTATION_LEARNING_WITH_METROPOLIS_HASTINGS_DOCK_H

#include <graph/api.h>

#include <utility.h>
#include <config.h>
#include <vertex.h>

namespace dynamic_graph_representation_learning_with_metropolis_hastings
{
    /**
     * Dock represents a structure that stores a graph as an augmented parallel balanced binary tree.
     * Keys in this tree are graph vertices and values are compressed edges, compressed walks, and metropolis hastings samplers.
     */
    class Dock
    {
        public:
            using Graph = aug_map<Vertex>;

            /**
             * Dock constructor.
             *
             * @param graph_vertices
             * @param graph_edges
             * @param offsets
             * @param edges
             * @param free_memory
             */
            Dock(long graph_vertices, long graph_edges, uintE* offsets, uintV* edges, bool free_memory = true)
            {
                #ifdef DOCK_TIMER
                    timer timer("Dock::Constructor");
                #endif

                // 1. Initialize memory pools
                Dock::init_memory_pools(graph_vertices, graph_edges);

                // 2. Create an empty vertex sequence
                using VertexStruct = std::pair<types::Vertex, VertexEntry>;
                auto vertices = pbbs::sequence<VertexStruct>(graph_vertices);

                // 3. In parallel construct graph vertices
                parallel_for(0, graph_vertices, [&](long index)
                {
                    size_t off = offsets[index];
                    size_t deg = ((index == (graph_vertices - 1)) ? graph_edges : offsets[index + 1]) - off;
                    auto S = pbbs::delayed_seq<uintV>(deg, [&](size_t j) { return edges[off + j]; });

                    if (deg > 0)
                        vertices[index] = std::make_pair(index, VertexEntry(types::CompressedEdges(S, index), dygrl::CompressedWalks()));
                    else
                        vertices[index] = std::make_pair(index, VertexEntry(types::CompressedEdges(), dygrl::CompressedWalks()));
                });

                // 4. Construct the graph
                auto replace = [](const VertexEntry& x, const VertexEntry& y) { return y; };
                this->graph_tree = Graph::Tree::multi_insert_sorted(nullptr, vertices.begin(), vertices.size(), replace, true);

                // 5. Determine random walk model
                // ...

                // 6. Memory cleanup
                if (free_memory)
                {
                    pbbs::free_array(offsets);
                    pbbs::free_array(edges);
                }

                vertices.clear();

                #ifdef DOCK_TIMER
                    timer.reportTotal("time(seconds)");
                #endif
            }

            /**
             * Number of vertices in a graph.
             *
             * @return the number of vertices in a graph
             */
            auto number_of_vertices() const
            {
                size_t n = this->graph_tree.size();
                auto last_vertex = this->graph_tree.select(n - 1);

                return n > 0 ? last_vertex.value.first + 1 : 0;
            }

            /**
             * Number of edges in a graph.
             *
             * @return the number of edges in a graph
             */
            auto number_of_edges() const
            {
                return this->graph_tree.aug_val();
            }

            /**
             * Flattens vertex tree to an array of vertex entries.
             *
             * @return the sequence of pointers to graph vertex entries
             */
            auto flatten_vertex_tree() const
            {
                #ifdef DOCK_TIMER
                    timer timer("Dock::FlattenVertexTree");
                #endif

                size_t n_vertices = this->number_of_vertices();
                auto flat_snapshot = pbbs::sequence<VertexEntry>(n_vertices);

                auto map_func = [&] (const Graph::E& entry, size_t ind)
                {
                    const uintV& key = entry.first;
                    const auto& value = entry.second;
                    flat_snapshot[key] = value;
                };

                this->map_vertices(map_func);

                #ifdef DOCK_TIMER
                    timer.reportTotal("time(seconds)");
                #endif

                return flat_snapshot;
            }

            /**
             * Traverses vertices and applies mapping function.
             *
             * @tparam F
             *
             * @param map_f
             * @param run_seq
             * @param granularity
             */
            template<class Function>
            void map_vertices(Function map_function, bool run_seq = false, size_t granularity = utils::node_limit) const
            {
                this->graph_tree.map_elms(map_function, run_seq, granularity);
            }

            /**
             * Destroys dock
             */
            void destroy()
            {
                this->graph_tree.~Graph();
                this->graph_tree.root = nullptr;
            }

            /**
             * Prints memory pool stats for the underlying lists.
             */
            void print_memory_pool_stats() const
            {
                // vertices memory pool stats
                std::cout << std::endl;
                std::cout << "Vertices tree memory pool stats: " << std::endl;
                Graph::print_stats();
                std::cout << std::endl;

                // edges memory pool stats
                std::cout << "Edges trees memory pool stats: " << std::endl;
                types::CompressedEdgesLists::print_stats();
                std::cout << std::endl;

                // walks memory pool stats
                std::cout << "Walks trees memory pool stats: " << std::endl;
                dygrl::CompressedWalksLists::print_stats();
                std::cout << std::endl;

                // compressed lists
                std::cout << "Compressed lists memory pool stats: " << std::endl;
                compressed_lists::print_stats();
                std::cout << std::endl;
            }

        private:
            Graph graph_tree;

            /**
            * Initializes memory pools for underlying lists.
            *
            * uses default sizeinit() s: 1 000 000 blocks, each block is put into a list of size 65 536 blocks
            * total allocated blocks = list_size(=65 5336) * (thread_count(=1..n) + ceil(allocated_blocks(=1M) / list_size(=65 536))
            * maximum blocks = ((3*RAM)/4)/size of one block.
            * @see list_allocator.h
            *
            * @param graph_vertices
            * @param graph_edges
            */
            static void init_memory_pools(size_t graph_vertices, size_t graph_edges)
            {
                types::CompressedEdgesLists::init();
                dygrl::CompressedWalksLists::init();

                compressed_lists::init(graph_vertices);

                Graph::init();
            }
    };
}

#endif // DYNAMIC_GRAPH_REPRESENTATION_LEARNING_WITH_METROPOLIS_HASTINGS_DOCK_H
