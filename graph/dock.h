#ifndef DYNAMIC_GRAPH_REPRESENTATION_LEARNING_WITH_METROPOLIS_HASTINGS_DOCK_H
#define DYNAMIC_GRAPH_REPRESENTATION_LEARNING_WITH_METROPOLIS_HASTINGS_DOCK_H

#include <graph/api.h>

#include <config.h>
#include <vertex.h>
#include <cuckoohash_map.hh>

#include <utility.h>
#include <pairings.h>

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
             * Destroys dock instance.
             */
            void destroy()
            {
                this->graph_tree.~Graph();
                this->graph_tree.root = nullptr;
            }

            // TODO: TBD
            void create_random_walks()
            {
                auto flat_graph     = this->flatten_vertex_tree();
                auto graph_vertices = this->number_of_vertices();

                auto walks  = graph_vertices * config::walks_per_vertex;
                auto cuckoo = libcuckoo::cuckoohash_map<types::Vertex, std::vector<types::Vertex>>();
                srand(time(nullptr));

                parallel_for(0, walks, [&](types::WalkID walk_id)
                {
                    auto current_vertex = walk_id % graph_vertices;
                    if (flat_graph[current_vertex].compressed_edges.degree() == 0) return;

                    for(types::Position position = 1; position <= config::walk_length; position++)
                    {
                        auto degree      = flat_graph[current_vertex].compressed_edges.degree();
                        auto neighbours  = flat_graph[current_vertex].compressed_edges.get_edges(current_vertex);
                        auto next_vertex = neighbours[rand() % degree];

                        cuckoo.upsert(current_vertex, [&](std::vector<types::Vertex>& walk_triplets)
                        {
                            types::PairedTriplet hash = pairings::Szudzik<types::Vertex>::pair({walk_id*config::walk_length + position, next_vertex});
                            walk_triplets.push_back(hash);
                        });

                        current_vertex = next_vertex;
                        pbbs::free_array(neighbours);
                    }
                });

                using VertexStruct = std::pair<types::Vertex, VertexEntry>;
                auto vertices = pbbs::sequence<VertexStruct>(graph_vertices);

                parallel_for(0, graph_vertices, [&](types::Vertex vertex)
                {
                    if (cuckoo.contains(vertex))
                    {
                        auto walk_parts = cuckoo.find(vertex);
                        auto sequence = pbbs::sequence<types::Vertex>(walk_parts.size());
                        for(auto i = 0; i < walk_parts.size(); i++)
                        {
                            sequence[i] = walk_parts[i];
                        }

                        sequence = pbbs::sample_sort(sequence, std::less<>());
                        vertices[vertex] = std::make_pair(vertex, VertexEntry(types::CompressedEdges(), dygrl::CompressedWalks(sequence, vertex)));
                    }
                    else
                    {
                        vertices[vertex] = std::make_pair(vertex, VertexEntry(types::CompressedEdges(), dygrl::CompressedWalks()));
                    }
                });

                auto replace = [&] (const uintV src, const VertexEntry& x, const VertexEntry& y)
                {
                    auto tree_plus = tree_plus::uniont(x.compressed_walks, y.compressed_walks, src);

                    // deallocate the memory
                    lists::deallocate(x.compressed_walks.plus);
                    tree_plus::Tree_GC::decrement_recursive(x.compressed_walks.root);
                    lists::deallocate(y.compressed_walks.plus);
                    tree_plus::Tree_GC::decrement_recursive(y.compressed_walks.root);

                    return VertexEntry(x.compressed_edges, CompressedWalks(tree_plus.plus, tree_plus.root));
                };

                this->graph_tree = Graph::Tree::multi_insert_sorted_with_values(this->graph_tree.root, vertices.begin(), vertices.size(), replace, true);
            }

            /**
             * Print memory footprint details.
             */
            void memory_footprint() const
            {
                std::cout << "----------------------------------------------------------------------------------"
                          << std::endl;

                size_t graph_vertices = this->number_of_vertices();
                size_t graph_edges    = this->number_of_edges();

                size_t vertex_node_size = Graph::node_size();
                size_t c_tree_node_size = types::CompressedTreesLists::node_size();

                size_t edges_heads = 0;
                size_t walks_heads = 0;
                size_t edges_bytes = 0;
                size_t walks_bytes = 0;
                size_t flat_graph_bytes = 0;
                auto flat_graph = this->flatten_vertex_tree();

                for (auto i = 0; i < flat_graph.size(); i++)
                {
                    flat_graph_bytes += sizeof(flat_graph[i]);

                    edges_heads += flat_graph[i].compressed_edges.edge_tree_nodes();
                    walks_heads += flat_graph[i].compressed_walks.edge_tree_nodes();

                    edges_bytes += flat_graph[i].compressed_edges.size_in_bytes(i);
                    walks_bytes += flat_graph[i].compressed_walks.size_in_bytes(i);
                }

                std::cout << "Total number of vertices: "
                          << graph_vertices
                          << " | Total number of edges: "
                          << graph_edges
                          << std::endl;

                std::cout << "Vertex tree: used nodes = "
                          << Graph::used_node()
                          << " -> node size = "
                          << vertex_node_size
                          << " -> memory used = "
                          << utility::MB(Graph::get_used_bytes())
                          << " MB"
                          << " = "
                          << utility::GB(Graph::get_used_bytes())
                          << " GB"
                          << std::endl;

                std::cout << "Edge compressed trees: heads = "
                          << edges_heads
                          << " -> node size = "
                          << c_tree_node_size
                          << " -> memory used = "
                          << utility::MB(edges_bytes + edges_heads*c_tree_node_size)
                          << " MB"
                          << " = "
                          << utility::GB(edges_bytes + edges_heads*c_tree_node_size)
                          << " GB"
                          << std::endl;

                std::cout << "Walk compressed trees: heads = "
                          << walks_heads
                          << " -> node size = "
                          << c_tree_node_size
                          << " -> memory used = "
                          << utility::MB(walks_bytes + walks_heads*c_tree_node_size)
                          << " MB"
                          << " = "
                          << utility::GB(walks_bytes + walks_heads*c_tree_node_size)
                          << " GB"
                          << std::endl;

                std::cout << "Flat snapshot size = "
                          << utility::MB(flat_graph_bytes)
                          << " MB = "
                          << utility::GB(flat_graph_bytes)
                          << " GB"
                          << std::endl << std::endl;

                size_t total_memory = Graph::get_used_bytes()
                        + walks_bytes + walks_heads*c_tree_node_size
                        + edges_bytes + edges_heads*c_tree_node_size;

                std::cout << "Total memory used: "
                          << utility::MB(total_memory)
                          << " MB = "
                          << utility::GB(total_memory)
                          << " GB"
                          << std::endl;

                std::cout << "----------------------------------------------------------------------------------"
                          << std::endl;
            }

            /**
             * Prints memory pool stats for the underlying lists.
             */
            void print_memory_pool_stats() const
            {
                // vertices memory pool stats
                std::cout << "----------------------------------------------------------------------------------"
                          << std::endl;
                std::cout << "Vertices tree memory lists: " << std::endl;
                std::cout << "\t"; Graph::print_stats();

                // edges and walks memory pool stats
                std::cout << "Edges & Walks trees memory lists: " << std::endl;
                std::cout << "\t"; types::CompressedTreesLists::print_stats();

                // compressed lists
                std::cout << "Pluses & Tails memory lists: " << std::endl;
                compressed_lists::print_stats();
                std::cout << "----------------------------------------------------------------------------------"
                          << std::endl;
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
                types::CompressedTreesLists::init();
                compressed_lists::init(graph_vertices);
                Graph::init();
            }
    };
}

#endif // DYNAMIC_GRAPH_REPRESENTATION_LEARNING_WITH_METROPOLIS_HASTINGS_DOCK_H
