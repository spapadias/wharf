#ifndef DYNAMIC_GRAPH_REPRESENTATION_LEARNING_WITH_METROPOLIS_HASTINGS_MALIN_H
#define DYNAMIC_GRAPH_REPRESENTATION_LEARNING_WITH_METROPOLIS_HASTINGS_MALIN_H

#include <graph/api.h>
#include <cuckoohash_map.hh>

#include <config.h>
#include <pairings.h>
#include <vertex.h>
#include <snapshot.h>

#include <models/deepwalk.h>
#include <models/node2vec.h>

namespace dynamic_graph_representation_learning_with_metropolis_hastings
{
    /**
     * @brief Malin represents a structure that stores a graph as an augmented parallel balanced binary tree.
     * Keys in this tree are graph vertices and values are compressed edges, compressed walks, and metropolis hastings samplers.
     */
    class Malin
    {
        public:
            using Graph = aug_map<dygrl::Vertex>;

            /**
             * @brief Malin constructor.
             *
             * @param graph_vertices - total vertices in a graph
             * @param graph_edges    - total edges in a graph
             * @param offsets        - vertex offsets for its neighbors
             * @param edges          - edges
             * @param free_memory    - free memory excess after graph is loaded
             */
            Malin(long graph_vertices, long graph_edges, uintE* offsets, uintV* edges, bool free_memory = true)
            {
                #ifdef MALIN_TIMER
                    timer timer("Malin::Constructor");
                #endif

                // 1. Initialize memory pools
                Malin::init_memory_pools(graph_vertices, graph_edges);

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
                        vertices[index] = std::make_pair(index, VertexEntry
                        (
            types::CompressedEdges(S, index),
                           dygrl::CompressedWalks(),
                           new dygrl::SamplerManager(0)
                        ));
                    else
                        vertices[index] = std::make_pair(index, VertexEntry
                        (
                            types::CompressedEdges(),
                            dygrl::CompressedWalks(),
                            new dygrl::SamplerManager(0)
                        ));
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

                #ifdef MALIN_TIMER
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
                #ifdef MALIN_TIMER
                    timer timer("Malin::FlattenVertexTree");
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

                #ifdef MALIN_TIMER
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
                #ifdef MALIN_TIMER
                    timer timer("Malin::FlattenGraph");
                #endif

                size_t n_vertices = this->number_of_vertices();
                auto flat_graph   = FlatGraph(n_vertices);

                auto map_func = [&] (const Graph::E& entry, size_t ind)
                {
                    const uintV& key  = entry.first;
                    const auto& value = entry.second;

                    flat_graph[key].neighbors = entry.second.compressed_edges.get_edges(key);
                    flat_graph[key].degree    = entry.second.compressed_edges.degree();
                    flat_graph[key].samplers  = entry.second.sampler_manager;
                };

                this->map_vertices(map_func);

                #ifdef MALIN_TIMER
                    timer.reportTotal("time(seconds)");

                    auto size = flat_graph.size_in_bytes();

                    std::cout << "Malin::FlattenGraph: Flat graph memory footprint: "
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
             * @param map_f   - map function
             * @param run_seq - determines whether to run part of the code sequantially
             * @param granularity
             */
            template<class Function>
            void map_vertices(Function map_function, bool run_seq = false, size_t granularity = utils::node_limit) const
            {
                this->graph_tree.map_elms(map_function, run_seq, granularity);
            }

            /**
             * @brief Destroys malin instance.
             */
            void destroy()
            {
                this->graph_tree.~Graph();
                this->graph_tree.root = nullptr;
            }

            /**
            * @brief Creates initial set of random walks.
            */
            void generate_initial_random_walks()
            {
                auto graph             = this->flatten_graph();
                auto total_vertices    = this->number_of_vertices();
                auto walks_to_generate = total_vertices * config::walks_per_vertex;
                auto cuckoo            = libcuckoo::cuckoohash_map<types::Vertex, std::vector<types::Vertex>>(total_vertices);

                using VertexStruct  = std::pair<types::Vertex, VertexEntry>;
                auto vertices       = pbbs::sequence<VertexStruct>(total_vertices);

                RandomWalkModel* model;
                switch (config::random_walk_model)
                {
                    case types::DEEPWALK:
                        model = new DeepWalk(&graph);
                        break;
                    case types::NODE2VEC:
                        model = new Node2Vec(&graph, config::paramP, config::paramQ);
                        break;
                    default:
                        std::cerr << "Unrecognized random walking model" << std::endl;
                        std::exit(1);
                }

                // 1. walk in parallel
                parallel_for(0, walks_to_generate, [&](types::WalkID walk_id)
                {
                    if (graph[walk_id % total_vertices].degree == 0)
                    {
                        types::PairedTriplet hash = pairings::Szudzik<types::Vertex>::pair({walk_id*config::walk_length, std::numeric_limits<uint32_t>::max() - 1});

                        cuckoo.insert(walk_id % total_vertices, std::vector<types::Vertex>());
                        cuckoo.update_fn(walk_id % total_vertices, [&](auto& vector)
                        {
                            vector.push_back(hash);
                        });

                        return;
                    }

                    types::State state  = model->initial_state(walk_id % total_vertices);

                    for(types::Position position = 0; position < config::walk_length; position++)
                    {
                        if (!graph[state.first].samplers->contains(state.second))
                        {
                            graph[state.first].samplers->insert(state.second, MetropolisHastingsSampler(state, model));
                        }

                        auto new_state = graph[state.first].samplers->find(state.second).sample(state, model);

                        if (!cuckoo.contains(state.first))
                            cuckoo.insert(state.first, std::vector<types::Vertex>());

                        types::PairedTriplet hash = (position != config::walk_length - 1) ?
                                pairings::Szudzik<types::Vertex>::pair({walk_id*config::walk_length + position, new_state.first}) :
                                pairings::Szudzik<types::Vertex>::pair({walk_id*config::walk_length + position, std::numeric_limits<uint32_t>::max() - 1});

                        cuckoo.update_fn(state.first, [&](auto& vector)
                        {
                            vector.push_back(hash);
                        });

                        state = new_state;
                    }
                });

                // 2. build vertices
                parallel_for(0, total_vertices, [&](types::Vertex vertex)
                {
                    if (cuckoo.contains(vertex))
                    {
                        auto triplets = cuckoo.find(vertex);
                        auto sequence = pbbs::sequence<types::Vertex>(triplets.size());

                        for(auto index = 0; index < triplets.size(); index++)
                            sequence[index] = triplets[index];

                        pbbs::sample_sort_inplace(pbbs::make_range(sequence.begin(), sequence.end()), std::less<>());
                        vertices[vertex] = std::make_pair(vertex, VertexEntry(types::CompressedEdges(), dygrl::CompressedWalks(sequence, vertex), new dygrl::SamplerManager(0)));
                    }
                    else
                    {
                        vertices[vertex] = std::make_pair(vertex, VertexEntry(types::CompressedEdges(), dygrl::CompressedWalks(), new dygrl::SamplerManager(0)));
                    }
                });

//                for(auto& e : vertices)
//                {
//                    std::cout << e.first << " ";
//                }
//
//                std::cout << std::endl;

                auto replace = [&] (const uintV& src, const VertexEntry& x, const VertexEntry& y)
                {
                    auto tree_plus = tree_plus::uniont(x.compressed_walks, y.compressed_walks, src);

                    // deallocate the memory
                    lists::deallocate(x.compressed_walks.plus);
                    tree_plus::Tree_GC::decrement_recursive(x.compressed_walks.root);
                    lists::deallocate(y.compressed_walks.plus);
                    tree_plus::Tree_GC::decrement_recursive(y.compressed_walks.root);

                    return VertexEntry(x.compressed_edges, CompressedWalks(tree_plus.plus, tree_plus.root), x.sampler_manager);
                };

                this->graph_tree = Graph::Tree::multi_insert_sorted_with_values(this->graph_tree.root, vertices.begin(), vertices.size(), replace, true);
                delete model;

            }

            /**
             * @brief Walks through the walk given walk id.
             *
             * @param walk_id - unique walk ID
             *
             * @return - walk string representation
             */
            std::string walk(types::WalkID walk_id)
            {
                // 1. Grab the first vertex in the walk
                types::Vertex current_vertex = walk_id % this->number_of_vertices();
                std::stringstream string_stream;
                types::Position position = 0;

                // 2. Walk
                while (current_vertex != std::numeric_limits<uint32_t>::max() - 1)
                {
                    string_stream << current_vertex << " ";

                    auto tree_node = this->graph_tree.find(current_vertex);

                    #ifdef MALIN_DEBUG
                        if (!tree_node.valid)
                        {
                            std::cerr << "Malin debug error! Malin::Walk::Vertex="
                                      << current_vertex << " is not found in the vertex tree!"
                                      << std::endl;

                            std::exit(1);
                        }
                    #endif

                    current_vertex = tree_node.value.compressed_walks.find_next(walk_id, position++, current_vertex);
                }

                return string_stream.str();
            }

            types::Vertex vertex_at_walk(types::WalkID walk_id, types::Position position)
            {
                types::Vertex current_vertex = walk_id % this->number_of_vertices();

                for (types::Position pos = 0; pos < position; pos++)
                {
                    auto tree_node = this->graph_tree.find(current_vertex);

                    #ifdef MALIN_DEBUG
                        if (!tree_node.valid)
                        {
                            std::cerr << "Malin debug error! Malin::Walk::Vertex="
                                      << current_vertex << " is not found in the vertex tree!"
                                      << std::endl;

                            std::exit(1);
                        }
                    #endif

                    current_vertex = tree_node.value.compressed_walks.find_next(walk_id, pos, current_vertex);
                }

                return current_vertex;
            }

            /**
            * @brief Inserts a batch of edges in the graph.
            *
            * @param m                  - size of the batch
            * @param edges              - atch of edges to insert
            * @param sorted             - sort the edges in the batch
            * @param remove_dups        - removes duplicate edges in the batch
            * @param nn
            * @param apply_walk_updates - decides if walk updates will be executed
            */
            pbbs::sequence<types::WalkID> insert_edges_batch(size_t m, std::tuple<uintV, uintV>* edges, bool sorted = false, bool remove_dups = false, size_t nn = std::numeric_limits<size_t>::max(), bool apply_walk_updates = true, bool run_seq = false)
            {
                auto fl = run_seq ? pbbs::fl_sequential : pbbs::no_flag;

                // 1. Set up
                using Edge = std::tuple<uintV, uintV>;

                auto edges_original = pbbs::make_range(edges, edges + m);
                Edge* edges_deduped = nullptr;

                // 2. Sort the edges in the batch (by source)
                if (!sorted)
                {
                    Malin::sort_edge_batch_by_source(edges, m, nn);
                }

                // 3. Remove duplicate edges
                if (remove_dups)
                {
                    // true when no duplicated edge, false otherwise
                    auto bool_seq = pbbs::delayed_seq<bool>(edges_original.size(), [&] (size_t i)
                    {
                        if(get<0>(edges_original[i]) == get<1>(edges_original[i])) return false;
                        return (i == 0 || edges_original[i] != edges_original[i-1]);
                    });

                    auto E = pbbs::pack(edges_original, bool_seq, fl); // creates a new pbbs::sequence
                    auto m_initial = m;                                // Initial number of generated edges
                    m = E.size();                                      // the size is not the same
                    auto m_final = m;                                  // Final number of edges in batch after duplicate removal
                    edges_deduped = E.to_array();                      // E afterwards is empty and nullptr
                }

                auto E = (edges_deduped) ? pbbs::make_range(edges_deduped, edges_deduped + m) : edges_original;

                // 4. Pack the starts vertices of edges
                auto start_im = pbbs::delayed_seq<size_t>(m, [&] (size_t i)
                {
                    return (i == 0 || (get<0>(E[i]) != get<0>(E[i-1])));
                });

                auto starts = pbbs::pack_index<size_t>(start_im, fl);
                size_t num_starts = starts.size();

                // 5. Build new wharf vertices
                using KV = std::pair<uintV, VertexEntry>;

                // decides to store whatf vertices on stack or heap
                constexpr const size_t stack_size = 20;
                KV kv_stack[stack_size];
                KV* new_verts = kv_stack;
                if (num_starts > stack_size)
                {
                    new_verts = pbbs::new_array<KV>(num_starts);
                }

                // pack the edges in the form: vertex_id - array of new edges
                parallel_for(0, num_starts, [&] (size_t i) {
                    size_t off = starts[i];
                    size_t deg = ((i == (num_starts-1)) ? m : starts[i+1]) - off;
                    uintV v = get<0>(E[starts[i]]);

                    auto S = pbbs::delayed_seq<uintV>(deg, [&] (size_t i) { return get<1>(E[off + i]); });

                    new_verts[i] = make_pair(v, VertexEntry(types::CompressedEdges(S, v, fl), dygrl::CompressedWalks(), new dygrl::SamplerManager(0)));
                });

                types::MapOfChanges rewalk_points = types::MapOfChanges();

                auto replace = [&, run_seq] (const intV& v, const VertexEntry& a, const VertexEntry& b)
                {
                    auto union_edge_tree = tree_plus::uniont(b.compressed_edges, a.compressed_edges, v, run_seq);

                    lists::deallocate(a.compressed_edges.plus);
                    tree_plus::Tree_GC::decrement_recursive(a.compressed_edges.root, run_seq);

                    lists::deallocate(b.compressed_edges.plus);
                    tree_plus::Tree_GC::decrement_recursive(b.compressed_edges.root, run_seq);

                    a.compressed_walks.iter_elms(v, [&](auto value)
                    {
                        auto pair = pairings::Szudzik<types::PairedTriplet>::unpair(value);

                        auto walk_id = pair.first / config::walk_length;
                        auto position = pair.first - (walk_id * config::walk_length);
                        auto next = pair.second;

                        if (!rewalk_points.template contains(walk_id))
                        {
                            rewalk_points.template insert(walk_id, std::make_tuple(position, v, false));
                        }
                        else
                        {

                            types::Position current_min_pos = get<0>(rewalk_points.find(walk_id));

                            if (current_min_pos > position)
                            {
                                rewalk_points.template update(walk_id, std::make_tuple(position, v, false));
                            }
                        }
                    });

                    return VertexEntry(union_edge_tree, a.compressed_walks, b.sampler_manager);
                };

                graph_update_time_on_insert.start();
                this->graph_tree = Graph::Tree::multi_insert_sorted_with_values(this->graph_tree.root, new_verts, num_starts, replace, true, run_seq);
                graph_update_time_on_insert.stop();

                walk_update_time_on_insert.start();
                auto affected_walks = pbbs::sequence<types::WalkID>(rewalk_points.size());
                if (apply_walk_updates) this->batch_walk_update(rewalk_points, affected_walks);
                walk_update_time_on_insert.stop();

                // 6. Deallocate memory
                if (num_starts > stack_size) pbbs::free_array(new_verts);
                if (edges_deduped)           pbbs::free_array(edges_deduped);

                #ifdef MALIN_DEBUG
                    std::cout << "Rewalk points (MapOfChanges): " << std::endl;

                    auto table = rewalk_points.lock_table();

                    for(auto& item : table)
                    {
                        std::cout << "Walk ID: " << item.first
                                  << " Position: "
                                  << (int) std::get<0>(item.second)
                                  << " Vertex: "
                                  << std::get<1>(item.second)
                                  << " Should reset: "
                                  << std::get<2>(item.second)
                                  << std::endl;
                    }

                    table.unlock();
                #endif

                return affected_walks;
            }

            /**
            * @brief Deletes a batch of edges from the graph.
            *
            * @param m - size of the batch
            * @param edges - batch of edges to delete
            * @param sorted - sort the edges in the batch
            * @param remove_dups - removes duplicate edges in the batch
            * @param nn
            * @param run_seq - decides if walk updates will be executed
            */
            pbbs::sequence<types::WalkID> delete_edges_batch(size_t m, tuple<uintV, uintV>* edges, bool sorted = false, bool remove_dups = false, size_t nn = std::numeric_limits<size_t>::max(), bool apply_walk_updates = true, bool run_seq = false)
            {
                auto fl = run_seq ? pbbs::fl_sequential : pbbs::no_flag;

                // 1. Set up
                using Edge = tuple<uintV, uintV>;

                auto edges_original = pbbs::make_range(edges, edges + m);
                Edge* edges_deduped = nullptr;

                // 2. Sort the edges in the batch (by source)
                if (!sorted)
                {
                    Malin::sort_edge_batch_by_source(edges, m, nn);
                }

                // 3. Remove duplicate edges
                if (remove_dups)
                {
                    // true when no duplicated edge, false otherwise
                    auto bool_seq = pbbs::delayed_seq<bool>(edges_original.size(), [&] (size_t i)
                    {
                        if(get<0>(edges_original[i]) == get<1>(edges_original[i])) return false;
                        return (i == 0 || edges_original[i] != edges_original[i-1]);
                    });

                    auto E = pbbs::pack(edges_original, bool_seq, fl); // creates a new pbbs::sequence
                    auto m_initial = m;                                // Initial number of generated edges
                    m = E.size();                                      // the size is not the same
                    auto m_final = m;                                  // Final number of edges in batch after duplicate removal
                    edges_deduped = E.to_array();                      // E afterwards is empty and nullptr
                }

                auto E = (edges_deduped) ? pbbs::make_range(edges_deduped, edges_deduped + m) : edges_original;

                // 4. Pack the starts vertices of edges
                auto start_im = pbbs::delayed_seq<size_t>(m, [&] (size_t i)
                {
                    return (i == 0 || (get<0>(E[i]) != get<0>(E[i-1])));
                });

                auto starts = pbbs::pack_index<size_t>(start_im, fl);
                size_t num_starts = starts.size();

                // 5. Build new wharf vertices
                using KV = std::pair<uintV, VertexEntry>;

                // decides to store whatf vertices on stack or heap
                constexpr const size_t stack_size = 20;
                KV kv_stack[stack_size];
                KV* new_verts = kv_stack;
                if (num_starts > stack_size)
                {
                    new_verts = pbbs::new_array<KV>(num_starts);
                }

                // pack the edges in the form: vertex_id - array of new edges
                parallel_for(0, num_starts, [&] (size_t i) {
                    size_t off = starts[i];
                    size_t deg = ((i == (num_starts-1)) ? m : starts[i+1]) - off;
                    uintV v = get<0>(E[starts[i]]);

                    auto S = pbbs::delayed_seq<uintV>(deg, [&] (size_t i) { return get<1>(E[off + i]); });

                    new_verts[i] = make_pair(v, VertexEntry(types::CompressedEdges(S, v, fl), dygrl::CompressedWalks(), new SamplerManager(0)));
                });

                types::MapOfChanges rewalk_points = types::MapOfChanges();

                auto replace = [&,run_seq] (const intV& v, const VertexEntry& a, const VertexEntry& b)
                {
                    auto difference_edge_tree = tree_plus::difference(b.compressed_edges, a.compressed_edges, v, run_seq);

                    lists::deallocate(a.compressed_edges.plus);
                    tree_plus::Tree_GC::decrement_recursive(a.compressed_edges.root, run_seq);

                    lists::deallocate(b.compressed_edges.plus);
                    tree_plus::Tree_GC::decrement_recursive(b.compressed_edges.root, run_seq);

                    bool should_reset = difference_edge_tree.degree() == 0;

                    a.compressed_walks.iter_elms(v, [&](auto value)
                    {
                        auto pair = pairings::Szudzik<types::PairedTriplet>::unpair(value);

                        auto walk_id = pair.first / config::walk_length;
                        auto position = pair.first - (walk_id * config::walk_length);
                        auto next = pair.second;

                        if (!rewalk_points.template contains(walk_id))
                        {
                            rewalk_points.template insert(walk_id, std::make_tuple(position, v, should_reset));
                        }
                        else
                        {
                            types::Position current_min_pos = get<0>(rewalk_points.find(walk_id));

                            if (current_min_pos > position)
                            {
                                rewalk_points.template update(walk_id, std::make_tuple(position, v, should_reset));
                            }
                        }
                    });

                    return VertexEntry(difference_edge_tree, a.compressed_walks, b.sampler_manager);
                };

                graph_update_time_on_delete.start();
                this->graph_tree = Graph::Tree::multi_insert_sorted_with_values(this->graph_tree.root, new_verts, num_starts, replace, true, run_seq);
                graph_update_time_on_delete.stop();

                walk_update_time_on_delete.start();
                auto affected_walks = pbbs::sequence<types::WalkID>(rewalk_points.size());
                if (apply_walk_updates) this->batch_walk_update(rewalk_points, affected_walks);
                walk_update_time_on_delete.stop();

                // 6. Deallocate memory
                if (num_starts > stack_size) pbbs::free_array(new_verts);
                if (edges_deduped) pbbs::free_array(edges_deduped);

                #ifdef MALIN_DEBUG
                    std::cout << "Rewalk points (MapOfChanges): " << std::endl;

                    auto table = rewalk_points.lock_table();

                    for(auto& item : table)
                    {
                        std::cout << "Walk ID: " << item.first
                                  << " Position: "
                                  << (int) std::get<0>(item.second)
                                  << " Vertex: "
                                  << std::get<1>(item.second)
                                  << " Should reset: "
                                  << std::get<2>(item.second)
                                  << std::endl;
                    }

                    table.unlock();
                #endif

                return affected_walks;
            }

            /**
             * @brief Updates affected walks in batch mode.
             *
             * @param types::MapOfChanges - rewalking points
             */
            void batch_walk_update(types::MapOfChanges& rewalk_points, pbbs::sequence<types::WalkID>& affected_walks)
            {
                types::ChangeAccumulator deletes = types::ChangeAccumulator();
                types::ChangeAccumulator inserts = types::ChangeAccumulator();

                uintV index = 0;

                for(auto& entry : rewalk_points.lock_table())
                {
                    affected_walks[index++] = entry.first;
                }

                auto graph = this->flatten_graph();
                RandomWalkModel* model;

                switch (config::random_walk_model)
                {
                    case types::DEEPWALK:
                        model = new DeepWalk(&graph);
                        break;
                    case types::NODE2VEC:
                        model = new Node2Vec(&graph, config::paramP, config::paramQ);
                        break;
                    default:
                        std::cerr << "Unrecognized random walking model!" << std::endl;
                        std::exit(1);
                }

                parallel_for(0, affected_walks.size(), [&](auto index)
                {
                    auto entry = rewalk_points.template find(affected_walks[index]);

                    auto current_position        = std::get<0>(entry);
                    auto current_vertex_old_walk = std::get<1>(entry);
                    auto should_reset            = std::get<2>(entry);

                    auto current_vertex_new_walk = current_vertex_old_walk;

                    if (should_reset)
                    {
                        current_position = 0;
                        current_vertex_old_walk = current_vertex_old_walk = affected_walks[index] % this->number_of_vertices();
                    }

                    fork_join_scheduler::Job delete_job = [&] ()
                    {
                        types::Position position = current_position;

                        while (current_vertex_old_walk != std::numeric_limits<uint32_t>::max() - 1)
                        {
                            auto vertex = this->graph_tree.find(current_vertex_old_walk);
                            auto next_old_walk = vertex.value.compressed_walks.find_next(affected_walks[index], position, current_vertex_old_walk);

                            types::PairedTriplet hash = pairings::Szudzik<types::Vertex>::pair({affected_walks[index]*config::walk_length + position, next_old_walk});

                            if (!deletes.contains(current_vertex_old_walk)) deletes.insert(current_vertex_old_walk, std::vector<types::PairedTriplet>());

                            deletes.update_fn(current_vertex_old_walk, [&](auto& vector)
                            {
                                vector.push_back(hash);
                            });

                            position++;
                            current_vertex_old_walk = next_old_walk;
                        }
                    }; delete_job();

                    fork_join_scheduler::Job insert_job = [&] ()
                    {
                        if (graph[current_vertex_new_walk].degree == 0)
                        {
                            types::PairedTriplet hash = pairings::Szudzik<types::Vertex>::pair({affected_walks[index]*config::walk_length, std::numeric_limits<uint32_t>::max() - 1});

                            if (!inserts.contains(current_vertex_new_walk)) inserts.insert(current_vertex_new_walk, std::vector<types::PairedTriplet>());
                            inserts.update_fn(current_vertex_new_walk, [&](auto& vector)
                            {
                                vector.push_back(hash);
                            });

                            return;
                        }

                        auto state = model->initial_state(current_vertex_new_walk);

                        if (config::random_walk_model == types::NODE2VEC && current_position > 0)
                        {
                            state.first  = current_vertex_new_walk;
                            state.second = this->vertex_at_walk(affected_walks[index], current_position - 1);
                        }

                        for (types::Position position = current_position; position < config::walk_length; position++)
                        {
                            if (!graph[state.first].samplers->contains(state.second))
                            {
                                graph[state.first].samplers->insert(state.second, MetropolisHastingsSampler(state, model));
                            }

                            state = graph[state.first].samplers->find(state.second).sample(state, model);

                            types::PairedTriplet hash = (position != config::walk_length - 1) ?
                                                        pairings::Szudzik<types::Vertex>::pair({affected_walks[index]*config::walk_length + position, state.first}) :
                                                        pairings::Szudzik<types::Vertex>::pair({affected_walks[index]*config::walk_length + position, std::numeric_limits<uint32_t>::max() - 1});


                            if (!inserts.contains(current_vertex_new_walk)) inserts.insert(current_vertex_new_walk, std::vector<types::PairedTriplet>());

                            inserts.update_fn(current_vertex_new_walk, [&](auto& vector)
                            {
                                vector.push_back(hash);
                            });

                            current_vertex_new_walk = state.first;
                        }

                    }; insert_job();
                });

                using VertexStruct  = std::pair<types::Vertex, VertexEntry>;
                auto insert_walks  = pbbs::sequence<VertexStruct>(inserts.size());
                auto delete_walks  = pbbs::sequence<VertexStruct>(deletes.size());

                fj.pardo([&]()
                {
                    auto ind = 0;

                    for(auto& item : inserts.lock_table())
                    {
                        auto sequence = pbbs::sequence<types::Vertex>(item.second.size());

                        for(auto i = 0; i < item.second.size(); i++)
                            sequence[i] = item.second[i];

                        pbbs::sample_sort_inplace(pbbs::make_range(sequence.begin(), sequence.end()), std::less<>());
                        insert_walks[ind++] = std::make_pair(item.first, VertexEntry(types::CompressedEdges(), dygrl::CompressedWalks(sequence, item.first), new dygrl::SamplerManager(0)));
                    }
                },
                [&]()
                {
                    auto ind = 0;

                    for(auto& item : deletes.lock_table())
                    {
                        auto sequence = pbbs::sequence<types::Vertex>(item.second.size());

                        for(auto i = 0; i < item.second.size(); i++)
                            sequence[i] = item.second[i];

                        pbbs::sample_sort_inplace(pbbs::make_range(sequence.begin(), sequence.end()), std::less<>());
                        delete_walks[ind++] = std::make_pair(item.first, VertexEntry(types::CompressedEdges(), dygrl::CompressedWalks(sequence, item.first), new dygrl::SamplerManager(0)));
                    }
                });

                fj.pardo([&]()
                {
                    pbbs::sample_sort_inplace(pbbs::make_range(delete_walks.begin(), delete_walks.end()), [&](auto& x, auto& y) {
                        return x.first < y.first;
                    });
                }, [&]()
                {
                    pbbs::sample_sort_inplace(pbbs::make_range(insert_walks.begin(), insert_walks.end()), [&](auto& x, auto& y) {
                        return x.first < y.first;
                    });
                });

                auto replaceD = [&] (const uintV& src, const VertexEntry& x, const VertexEntry& y)
                {
                    auto tree_plus = tree_plus::difference(y.compressed_walks, x.compressed_walks, src);

                    // deallocate the memory
                    lists::deallocate(x.compressed_walks.plus);
                    tree_plus::Tree_GC::decrement_recursive(x.compressed_walks.root);
                    lists::deallocate(y.compressed_walks.plus);
                    tree_plus::Tree_GC::decrement_recursive(y.compressed_walks.root);

                    return VertexEntry(x.compressed_edges, dygrl::CompressedWalks(tree_plus.plus, tree_plus.root), x.sampler_manager);
                };

                this->graph_tree = Graph::Tree::multi_insert_sorted_with_values(this->graph_tree.root, delete_walks.begin(), delete_walks.size(), replaceD, true);

                auto replaceI = [&] (const uintV& src, const VertexEntry& x, const VertexEntry& y)
                {
                    auto tree_plus = tree_plus::uniont(x.compressed_walks, y.compressed_walks, src);

                    // deallocate the memory
                    lists::deallocate(x.compressed_walks.plus);
                    tree_plus::Tree_GC::decrement_recursive(x.compressed_walks.root);
                    lists::deallocate(y.compressed_walks.plus);
                    tree_plus::Tree_GC::decrement_recursive(y.compressed_walks.root);

                    return VertexEntry(x.compressed_edges, dygrl::CompressedWalks(tree_plus.plus, tree_plus.root), x.sampler_manager);
                };

                this->graph_tree = Graph::Tree::multi_insert_sorted_with_values(this->graph_tree.root, insert_walks.begin(), insert_walks.size(), replaceI, true);
            }

            /**
             * @brief Prints memory footprint details.
             */
            void memory_footprint() const
            {
                std::cout << std::endl;

                size_t graph_vertices   = this->number_of_vertices();
                size_t graph_edges      = this->number_of_edges();

                size_t vertex_node_size = Graph::node_size();
                size_t c_tree_node_size = types::CompressedTreesLists::node_size();

                size_t edges_heads      = 0;
                size_t walks_heads      = 0;
                size_t edges_bytes      = 0;
                size_t walks_bytes      = 0;
                size_t samplers_bytes   = 0;
                size_t flat_graph_bytes = 0;
                auto flat_graph         = this->flatten_vertex_tree();

                for (auto i = 0; i < flat_graph.size(); i++)
                {
                    flat_graph_bytes += sizeof(flat_graph[i]);

                    edges_heads += flat_graph[i].compressed_edges.edge_tree_nodes();
                    walks_heads += flat_graph[i].compressed_walks.edge_tree_nodes();

                    edges_bytes += flat_graph[i].compressed_edges.size_in_bytes(i);
                    walks_bytes += flat_graph[i].compressed_walks.size_in_bytes(i);

                    for (auto& entry : flat_graph[i].sampler_manager->lock_table())
                    {
                        samplers_bytes += sizeof(entry.first) + sizeof(entry.second);
                    }
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

                std::cout << "Samplers: \n\t"
                          << "Total memory usage: " << utility::MB(samplers_bytes)
                          << " MB = " << utility::GB(samplers_bytes)
                          << " GB" << std::endl;

                std::cout << "Flat graph: \n\t"
                          << "Total memory usage: " << utility::MB(flat_graph_bytes)
                          << " MB = " << utility::GB(flat_graph_bytes)
                          << " GB" << std::endl;

                size_t total_memory = Graph::get_used_bytes()
                        + walks_bytes + walks_heads*c_tree_node_size
                        + edges_bytes + edges_heads*c_tree_node_size + samplers_bytes;

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

            /**
            * Sorts a sequence of batch updates.
            *
            * @param edges       sequence of edges (src, dst) to be sorted by src
            * @param batch_edges total number of edges to be sorted
            * @param nn
            */
            static void sort_edge_batch_by_source(std::tuple<uintV, uintV>* edges, size_t batch_edges, size_t nn = std::numeric_limits<size_t>::max())
            {
                #ifdef MALIN_TIMER
                    timer timer("Malin::SortEdgeBatchBySource");
                #endif

                // 1. Set up
                using Edge = tuple<uintV, uintV>;

                auto edges_original = pbbs::make_range(edges, edges + batch_edges);
                size_t vertex_bits = pbbs::log2_up(nn);     // number of bits to represent a vertex in the graph

                // 2. Induce nn if not given (captures x < number of nodes < y such that x and y are powers of 2)
                if (nn == std::numeric_limits<size_t>::max())
                {
                    auto max_edge_id = pbbs::delayed_seq<size_t>(batch_edges, [&](size_t i)
                    {
                        return std::max(std::get<0>(edges_original[i]), std::get<1>(edges_original[i]));
                    });

                    vertex_bits = pbbs::log2_up(pbbs::reduce(max_edge_id, pbbs::maxm<size_t>()));
                    nn = 1 << vertex_bits;
                }

                // 3. Sort edges by source
                auto edge_to_long = [nn, vertex_bits](Edge e) -> size_t {
                    return (static_cast<size_t>(std::get<0>(e)) << vertex_bits) + static_cast<size_t>(std::get<1>(e));
                };

    //                auto edge_ids_log = pbbs::delayed_seq<size_t>(batch_edges, [&](size_t i) {
    //                    return pbbs::log2_up(edge_to_long(edges_original[i]));
    //                });

                size_t bits = 2 * vertex_bits;

                // Only apply integer sort if it will be work-efficient
                if (nn <= (batch_edges * pbbs::log2_up(batch_edges)))
                {
                    pbbs::integer_sort_inplace(edges_original, edge_to_long, bits);
                }
                else
                {
                    pbbs::sample_sort_inplace(edges_original, std::less<>());
                }

                #ifdef MALIN_TIMER
                    timer.reportTotal("time (seconds)");
                #endif
            }
    };
}

#endif // DYNAMIC_GRAPH_REPRESENTATION_LEARNING_WITH_METROPOLIS_HASTINGS_MALIN_H
