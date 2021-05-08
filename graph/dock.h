#ifndef DYNAMIC_GRAPH_REPRESENTATION_LEARNING_WITH_METROPOLIS_HASTINGS_DOCK_H
#define DYNAMIC_GRAPH_REPRESENTATION_LEARNING_WITH_METROPOLIS_HASTINGS_DOCK_H

#include <graph/api.h>
#include <cuckoohash_map.hh>

#include <config.h>
#include <utility.h>
#include <pairings.h>
#include <vertex.h>
#include <snapshot.h>

#include <models/deepwalk.h>

namespace dynamic_graph_representation_learning_with_metropolis_hastings
{
    /**
     * @brief Dock represents a structure that stores a graph as an augmented parallel balanced binary tree.
     * Keys in this tree are graph vertices and values are compressed edges, compressed walks, and metropolis hastings samplers.
     */
    class Dock
    {
        public:
            using Graph = aug_map<dygrl::Vertex>;

            /**
             * @brief Dock constructor.
             *
             * @param graph_vertices - total vertices in a graph
             * @param graph_edges    - total edges in a graph
             * @param offsets        - vertex offsets for its neighbors
             * @param edges          - edges
             * @param free_memory    - free memory excess after graph is loaded
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
                        vertices[index] = std::make_pair(index, VertexEntry(types::CompressedEdges(S, index), dygrl::CompressedWalks(), dygrl::SamplerManager(0)));
                    else
                        vertices[index] = std::make_pair(index, VertexEntry(types::CompressedEdges(), dygrl::CompressedWalks(), dygrl::SamplerManager(0)));
                });

                // 4. Construct the graph
                auto replace = [](const VertexEntry& x, const VertexEntry& y) { return y; };
                this->graph_tree = Graph::Tree::multi_insert_sorted(nullptr, vertices.begin(), vertices.size(), replace, true);

                // 5. Memory cleanup
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
             * @brief Number of vertices in a graph.
             *
             * @return - the number of vertices in a graph
             */
            [[nodiscard]] auto number_of_vertices() const
            {
                size_t n = this->graph_tree.size();
                auto last_vertex = this->graph_tree.select(n - 1);

                return n > 0 ? last_vertex.value.first + 1 : 0;
            }

            /**
             * @brief Number of edges in a graph.
             *
             * @return - the number of edges in a graph
             */
            [[nodiscard]] auto number_of_edges() const
            {
                return this->graph_tree.aug_val();
            }

            /**
             * @brief Flattens the vertex tree to an array of vertex entries.
             *
             * @return - the sequence of pointers to graph vertex entries
             */
            [[nodiscard]] FlatVertexTree flatten_vertex_tree() const
            {
                #ifdef DOCK_TIMER
                    timer timer("Dock::FlattenVertexTree");
                #endif

                types::Vertex n_vertices = this->number_of_vertices();
                auto flat_vertex_tree    = FlatVertexTree(n_vertices);

                auto map_func = [&] (const Graph::E& entry, size_t ind)
                {
                    const types::Vertex& key = entry.first;
                    const auto& value = entry.second;
                    flat_vertex_tree[key] = value;
                };

                this->map_vertices(map_func);

                #ifdef DOCK_TIMER
                    timer.reportTotal("time(seconds)");

                    std::cout << "Flat vertex tree memory footprint: "
                              << utility::MB(flat_vertex_tree.size_in_bytes())
                              << " MB = " << utility::GB(flat_vertex_tree.size_in_bytes())
                              << " GB" << std::endl;
                #endif

                return flat_vertex_tree;
            }

            /**
            * @brief Flattens the graph to an array of vertices, their degrees, neighbors, and sampler managers.
            *
            * @return - the sequence of vertices, their degrees, neighbors, and sampler managers
            */
            [[nodiscard]] FlatGraph flatten_graph() const
            {
                #ifdef DOCK_TIMER
                    timer timer("Dock::FlattenGraph");
                #endif

                size_t n_vertices = this->number_of_vertices();
                auto flat_graph   = FlatGraph(n_vertices);

                auto map_func = [&] (const Graph::E& entry, size_t ind)
                {
                    const uintV& key  = entry.first;
                    const auto& value = entry.second;

                    flat_graph[key].neighbors = entry.second.compressed_edges.get_edges(key);
                    flat_graph[key].degrees   = entry.second.compressed_edges.degree();
                    flat_graph[key].samplers  = entry.second.sampler_manager;
                };

                this->map_vertices(map_func);

                #ifdef DOCK_TIMER
                    timer.reportTotal("time(seconds)");

                    auto size = flat_graph.size_in_bytes();

                    std::cout << "Dock::FlattenGraph: Flat vertex tree memory footprint: "
                              << utility::MB(size)
                              << " MB = " << utility::GB(size)
                              << " GB" << std::endl;
                #endif

                return flat_graph;
            }

            /**
             * @brief Traverses vertices and applies mapping function.
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
             * @brief Destroys dock instance.
             */
            void destroy()
            {
                this->graph_tree.~Graph();
                this->graph_tree.root = nullptr;
            }

            // TODO : working on it
            void create_random_walks()
            {
                auto graph          = this->flatten_graph();
                auto total_vertices = this->number_of_vertices();
                auto walks          = total_vertices * config::walks_per_vertex;
                auto cuckoo         = libcuckoo::cuckoohash_map<types::Vertex, std::vector<types::Vertex>>(total_vertices);

                using VertexStruct  = std::pair<types::Vertex, VertexEntry>;
                auto vertices       = pbbs::sequence<VertexStruct>(total_vertices);
                RandomWalkModel* model;

                switch (config::random_walk_model)
                {
                    case types::DEEPWALK:
                        model = new DeepWalk(&graph);
                        break;
                    case types::NODE2VEC:
                        std::cerr << "Unrecognized random walking model" << std::endl;
                        std::exit(1);
                    default:
                        std::cerr << "Unrecognized random walking model" << std::endl;
                        std::exit(1);
                }

                parallel_for(0, walks, [&](types::WalkID walk_id)
                {
                    types::State state  = model->initial_state(walk_id % total_vertices);
                    if (graph[state.first].degrees == 0) return;

                    for(types::Position position = 0; position < config::walk_length; position++)
                    {
                        std::cout << state.first << " ";

                        if (!graph[state.first].samplers->contains(state.second))
                        {
                            graph[state.first].samplers->insert(state.second, MetropolisHastingsSampler(state, model));
                        }

                        auto new_state = graph[state.first].samplers->find(state.second).sample(state, model);

//                        if (!cuckoo.contains(state.first))
//                            cuckoo.insert(state.first, std::vector<types::Vertex>());
//
//                        cuckoo.update_fn(state.first, [&](auto& vector)
//                        {
//                            types::PairedTriplet hash = pairings::Szudzik<types::Vertex>::pair({walk_id*config::walk_length + position, new_state.first});
//                            vector.push_back(hash);
//                        });

                        state = new_state;
                    }

                    std::cout << std::endl;
                });

//                parallel_for(0, total_vertices, [&](types::Vertex vertex)
//                {
//                    if (cuckoo.contains(vertex))
//                    {
//                        auto triplets = cuckoo.find(vertex);
//                        auto sequence = pbbs::sequence<types::Vertex>(triplets.size());
//
//                        for(auto index = 0; index < triplets.size(); index++)
//                            sequence[index] = triplets[index];
//
//                        pbbs::sample_sort_inplace(pbbs::make_range(sequence.begin(), sequence.end()), std::less<>());
//                        vertices[vertex] = std::make_pair(vertex, VertexEntry(types::CompressedEdges(), dygrl::CompressedWalks(sequence, vertex), SamplerManager(0)));
//                    }
//                    else
//                    {
//                        vertices[vertex] = std::make_pair(vertex, VertexEntry(types::CompressedEdges(), dygrl::CompressedWalks(), SamplerManager(0)));
//                    }
//                });

//                auto replace = [&] (const uintV src, const VertexEntry& x, const VertexEntry& y)
//                {
//                    auto tree_plus = tree_plus::uniont(x.compressed_walks, y.compressed_walks, src);
//
//                    // deallocate the memory
//                    lists::deallocate(x.compressed_walks.plus);
//                    tree_plus::Tree_GC::decrement_recursive(x.compressed_walks.root);
//                    lists::deallocate(y.compressed_walks.plus);
//                    tree_plus::Tree_GC::decrement_recursive(y.compressed_walks.root);
//
//                    return VertexEntry(x.compressed_edges, CompressedWalks(tree_plus.plus, tree_plus.root), *x.sampler_manager);
//                };
//
//                this->graph_tree = Graph::Tree::multi_insert_sorted_with_values(this->graph_tree.root, vertices.begin(), vertices.size(), replace, true);

                delete model;
            }

            /**
             * @brief Walks through the walk given walk id.
             *
             * @param walk_id - unique walk ID
             *
             * @return - walk string representation
             */
            std::string rewalk(types::WalkID walk_id)
            {
                // 1. Grab the first vertex in the walk
                types::Vertex current_vertex = walk_id % this->number_of_vertices();

                std::stringstream string_stream;
                string_stream << "WalkID: " << walk_id << " | ";

                // 2. Rewalk
                for (types::Position position = 0; position < config::walk_length; position++)
                {
                    position + 1 == config::walk_length ? string_stream << current_vertex : string_stream << current_vertex << ", ";

                    auto tree_node = this->graph_tree.find(current_vertex);

                    #ifdef DOCK_DEBUG
                        if (!tree_node.valid)
                        {
                            std::cerr << "Dock debug error! Dock::Rewalk::Vertex="
                                      << current_vertex << " is not found in the vertex tree!"
                                      << std::endl;

                            std::exit(1);
                        }
                    #endif

                    if (tree_node.value.compressed_edges.degree() == 0) break;
                    current_vertex = tree_node.value.compressed_walks.find_next(walk_id, position, current_vertex);
                }

                return string_stream.str();
            }

            /**
             * @brief Prints memory footprint details.
             */
            void memory_footprint() const
            {
                std::cout << std::endl;

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

                std::cout << "Graph: \n\t" << "Vertices: " << graph_vertices << ", Edges: " << graph_edges << std::endl;

                std::cout << "Vertex Tree: \n\t"
                          << "Heads: " << Graph::used_node()
                          << ", Head size: " << vertex_node_size
                          << ", Memory usage: " << utility::MB(Graph::get_used_bytes()) << " MB"
                          << " = " << utility::GB(Graph::get_used_bytes()) << " GB" << std::endl;

                std::cout << "Edge Trees: \n\t"
                          << "Heads: " << edges_heads
                          << ", Head size: " << c_tree_node_size
                          << ", Lists memory: " << utility::MB(edges_bytes) << " MB"
                          << " = " << utility::GB(edges_bytes) << " GB"
                          << ", Total memory usage: " << utility::MB(edges_bytes + edges_heads*c_tree_node_size)
                          << " MB = " << utility::GB(edges_bytes + edges_heads*c_tree_node_size)
                          << " GB" << std::endl;

                std::cout << "Walks Trees: \n\t"
                          << "Heads: " << walks_heads
                          << ", Head size: " << c_tree_node_size
                          << ", Lists memory: " << utility::MB(walks_bytes) << " MB"
                          << " = " << utility::GB(walks_bytes) << " GB"
                          << ", Total memory usage: " << utility::MB(walks_bytes + walks_heads*c_tree_node_size)
                          << " MB = " << utility::GB(walks_bytes + walks_heads*c_tree_node_size)
                          << " GB" << std::endl;

                std::cout << "Flat graph: \n\t"
                          << "Total memory usage: " << utility::MB(flat_graph_bytes)
                          << " MB = " << utility::GB(flat_graph_bytes)
                          << " GB" << std::endl;

                size_t total_memory = Graph::get_used_bytes()
                        + walks_bytes + walks_heads*c_tree_node_size
                        + edges_bytes + edges_heads*c_tree_node_size;

                std::cout << edges_bytes << " " << walks_bytes << " " << edges_bytes + walks_bytes  << std::endl;

                std::cout << "Total memory used: \n\t" << utility::MB(total_memory) << " MB = "
                          << utility::GB(total_memory) << " GB" << std::endl;

                std::cout << std::endl;
            }

            /**
             * @brief Prints memory pool stats for the underlying lists.
             */
            void print_memory_pool_stats() const
            {
                std::cout << std::endl;

                // vertices memory pool stats
                std::cout << "Vertices tree memory lists: \n\t";
                Graph::print_stats();

                // edges and walks memory pool stats
                std::cout << "Edges & Walks trees memory lists: \n\t";
                types::CompressedTreesLists::print_stats();

                // compressed lists
                std::cout << "Pluses & Tails memory lists: \n";
                compressed_lists::print_stats();

                std::cout << std::endl;
            }

        private:
            Graph graph_tree;

            /**
            * Initializes memory pools for underlying lists.
            *
            * uses default size init(): 1 000 000 blocks, each block is put into a list of size 65 536 blocks
            * total allocated blocks = list_size(=65 5336) * (thread_count(=1..n) + ceil(allocated_blocks(=1M) / list_size(=65 536))
            * maximum blocks = ((3*RAM)/4)/size of one block.
            * @see list_allocator.h
            *
            * @param graph_vertices - total graph vertices
            * @param graph_edges    - total graph edges
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
