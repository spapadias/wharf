#ifndef WHARF_H
#define WHARF_H

#include <graph/api.h>
#include <cuckoohash_map.hh>
#include <pbbslib/utilities.h>

#include <config.h>
#include <pairings.h>
#include <vertex.h>
#include <snapshot.h>

#include <models/deepwalk.h>
#include <models/node2vec.h>
#include <set>

namespace dynamic_graph_representation_learning_with_metropolis_hastings
{
    /**
     * @brief Wharf represents a structure that stores a graph as an augmented parallel balanced binary tree.
     * Keys in this tree are graph vertices and values are compressed edges, compressed walks, and metropolis hastings samplers.
     */
    class Wharf
    {
        public:
            using Graph = aug_map<dygrl::Vertex>;
            std::atomic<unsigned long> number_of_sampled_vertices;

            /**
             * @brief Wharf constructor.
             *
             * @param graph_vertices - total vertices in a graph
             * @param graph_edges    - total edges in a graph
             * @param offsets        - vertex offsets for its neighbors
             * @param edges          - edges
             * @param free_memory    - free memory excess after graph is loaded
             */
            Wharf(long graph_vertices, long graph_edges, uintE* offsets, uintV* edges, bool free_memory = true)
            {
                #ifdef WHARF_TIMER
                    timer timer("Wharf::Constructor");
                #endif

                // 1. Initialize memory pools
                Wharf::init_memory_pools(graph_vertices, graph_edges);

                // 2. Create an empty vertex sequence
                using VertexStruct = std::pair<types::Vertex, VertexEntry>;
                auto vertices = pbbs::sequence<VertexStruct>(graph_vertices);

                // 3. In parallel construct graph vertices
                parallel_for(0, graph_vertices, [&](long index)
                {
                    size_t off = offsets[index];
                    size_t deg = ((index == (graph_vertices - 1)) ? graph_edges : offsets[index + 1]) - off;
                    auto S = pbbs::delayed_seq<uintV>(deg, [&](size_t j) { return edges[off + j]; });

					vector<dygrl::CompressedWalks> vec_compwalk;
                    if (deg > 0)
                        vertices[index] = std::make_pair(index, VertexEntry(types::CompressedEdges(S, index),
																			vec_compwalk,
                                                                            new dygrl::SamplerManager(0)));
                    else
                        vertices[index] = std::make_pair(index, VertexEntry(types::CompressedEdges(),
                                                                            vec_compwalk,
                                                                            new dygrl::SamplerManager(0)));
                });

                // 4. Construct the graph
                auto replace = [](const VertexEntry& x, const VertexEntry& y) { return y; };
                this->graph_tree = Graph::Tree::multi_insert_sorted(nullptr, vertices.begin(), vertices.size(), replace, true);

				for (auto i = 0; i < 250; i++)
					MAVS2.push_back(types::MapAffectedVertices());

				number_of_sampled_vertices = 0;

                // 5. Memory cleanup
                if (free_memory)
                {
                    pbbs::free_array(offsets);
                    pbbs::free_array(edges);
                }

                vertices.clear();

                #ifdef WHARF_TIMER
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
                #ifdef WHARF_TIMER
                    timer timer("Wharf::FlattenVertexTree");
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

                #ifdef WHARF_TIMER
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
                #ifdef WHARF_TIMER
                    timer timer("Wharf::FlattenGraph");
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

                #ifdef WHARF_TIMER
                    timer.reportTotal("time(seconds)");

                    auto size = flat_graph.size_in_bytes();

                    std::cout << "Wharf::FlattenGraph: Flat graph memory footprint: "
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
             * @brief Destroys inverted index.
             */
            void destroy_index()
            {
                auto map_func = [&] (Graph::E& entry, size_t ind)
                {
					for (auto cw : entry.second.compressed_walks)
					{
						if (cw.plus)
						{
							lists::deallocate(cw.plus);
							cw.plus = nullptr;
						}

						// delete compresed walks - tree part
						if (cw.root)
						{
							auto T = walk_plus::edge_list();

							T.root = cw.root;
							cw.root = nullptr;
						}
					}

                    // delete the sample manager
                    entry.second.sampler_manager->clear();
                };

                this->map_vertices(map_func);
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
                        types::PairedTriplet hash = pairings::Szudzik<types::Vertex>::pair({walk_id * config::walk_length + 0, walk_id % total_vertices});
                        cuckoo.insert(walk_id % total_vertices, std::vector<types::Vertex>());
                        cuckoo.update_fn(walk_id % total_vertices, [&](auto& vector) {
                            vector.push_back(hash);
                        });

                        return;
                    }

                    auto random = config::random;
                    if (config::deterministic_mode)
                        random = utility::Random(walk_id / total_vertices);
                    types::State state  = model->initial_state(walk_id % total_vertices);

                    for(types::Position position = 0; position < config::walk_length; position++)
                    {
                        if (!graph[state.first].samplers->contains(state.second))
                            graph[state.first].samplers->insert(state.second, MetropolisHastingsSampler(state, model));

                        auto new_state = graph[state.first].samplers->find(state.second).sample(state, model);
                        if (config::deterministic_mode)
                            new_state = model->new_state(state, graph[state.first].neighbors[random.irand(graph[state.first].degree)]);

                        if (!cuckoo.contains(state.first))
                            cuckoo.insert(state.first, std::vector<types::Vertex>());

                        types::PairedTriplet hash = (position != config::walk_length - 1) ?
                                pairings::Szudzik<types::Vertex>::pair({walk_id * config::walk_length + position, new_state.first}) :
                                pairings::Szudzik<types::Vertex>::pair({walk_id * config::walk_length + position, state.first});

                        cuckoo.update_fn(state.first, [&](auto& vector) {
                            vector.push_back(hash);
                        });

                        // Assign the new state to the sampler
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
						// assume the initial random walks are created at batch 0
						vector<dygrl::CompressedWalks> vec_compwalks;
						vec_compwalks.push_back(dygrl::CompressedWalks(sequence, vertex, 666, 666, 0));
                        vertices[vertex] = std::make_pair(vertex, VertexEntry(types::CompressedEdges(), vec_compwalks, new dygrl::SamplerManager(0)));
                    }
                    else
                    {
	                    vector<dygrl::CompressedWalks> vec_compwalks; vec_compwalks.push_back(dygrl::CompressedWalks(0)); // at batch 0
	                    vertices[vertex] = std::make_pair(vertex, VertexEntry(types::CompressedEdges(), vec_compwalks,new dygrl::SamplerManager(0)));
                    }
                });

                auto replace = [&] (const uintV& src, const VertexEntry& x, const VertexEntry& y)
                {
					auto x_prime = x;
					x_prime.compressed_walks.push_back(y.compressed_walks.back()); // y has only one walk tree

					return x_prime;
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
            std::string walk_simple_find(types::WalkID walk_id)
            {
                // 1. Grab the first vertex in the walk
                types::Vertex current_vertex = walk_id % this->number_of_vertices();
                std::stringstream string_stream;
                types::Position position = 0;

                // 2. Walk
                types::Vertex previous_vertex;
                while (true)
                {
                    string_stream << current_vertex << " ";

                    auto tree_node = this->graph_tree.find(current_vertex);

                    #ifdef WHARF_DEBUG
                        if (!tree_node.valid)
                        {
                            std::cerr << "Wharf debug error! Wharf::Walk::Vertex="
                                      << current_vertex << " is not found in the vertex tree!"
                                      << std::endl;

                            std::exit(1);
                        }
                    #endif

                    // Cache the previous current vertex before going to the next
                    previous_vertex = current_vertex;

                    // uses only simple find_next
                    current_vertex = tree_node.value.compressed_walks.front().find_next(walk_id, position++, current_vertex); // operate on the front() as only one walk-tree exists after merging

                    if (current_vertex == previous_vertex)
                        break;
                }

                return string_stream.str();
            }

            /**
             * @brief Walks through the walk given walk id.
             *
             * @param walk_id - unique walk ID
             *
             * @return - walk string representation
             */
            void walk_silent(types::WalkID walk_id)
            {
                // 1. Grab the first vertex in the walk
                types::Vertex current_vertex = walk_id % this->number_of_vertices();
                std::stringstream string_stream;
                types::Position position = 0;

                // 2. Walk
                types::Vertex previous_vertex;
                while (true)
                {
                    auto tree_node = this->graph_tree.find(current_vertex);

                    #ifdef WHARF_DEBUG
                        if (!tree_node.valid)
                        {
                            std::cerr << "Wharf debug error! Wharf::Walk::Vertex="
                                      << current_vertex << " is not found in the vertex tree!"
                                      << std::endl;

                            std::exit(1);
                        }
                    #endif

                    // Cache the previous current vertex before going to the next
                    previous_vertex = current_vertex;

                    // uses only simple find_next
					if (config::range_search_mode)
						current_vertex = tree_node.value.compressed_walks.front().find_next_in_range(walk_id, position++, current_vertex);
					else
	                    current_vertex = tree_node.value.compressed_walks.front().find_next(walk_id, position++, current_vertex);

                    if (current_vertex == previous_vertex)
                        break;
                }
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
            pbbs::sequence<types::WalkID> insert_edges_batch(size_t m,
                                                             std::tuple<uintV, uintV>* edges,
															 int batch_num,
                                                             bool sorted = false,
                                                             bool remove_dups = false,
                                                             size_t nn = std::numeric_limits<size_t>::max(),
                                                             bool apply_walk_updates = true,
                                                             bool run_seq = false)
            {
                auto fl = run_seq ? pbbs::fl_sequential : pbbs::no_flag;

                // 1. Set up
                using Edge = std::tuple<uintV, uintV>;

                auto edges_original = pbbs::make_range(edges, edges + m);
                Edge* edges_deduped = nullptr;

                // 2. Sort the edges in the batch (by source)
                if (!sorted)
                {
                    Wharf::sort_edge_batch_by_source(edges, m, nn);
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

                    auto E = pbbs::pack(edges_original, bool_seq, fl); // Creates a new pbbs::sequence
                    auto m_initial = m;                                // Initial number of generated edges
                    m = E.size();                                      // The size is not the same
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

                // Decides to store Wharf vertices on stack or heap
                constexpr const size_t stack_size = 20;
                KV kv_stack[stack_size];
                KV* new_verts = kv_stack;
                if (num_starts > stack_size)
                    new_verts = pbbs::new_array<KV>(num_starts);

                // pack the edges in the form: vertex_id - array of new edges
                parallel_for(0, num_starts, [&] (size_t i) {
                    size_t off = starts[i];
                    size_t deg = ((i == (num_starts-1)) ? m : starts[i+1]) - off;
                    uintV v = get<0>(E[starts[i]]);

                    auto S = pbbs::delayed_seq<uintV>(deg, [&] (size_t i) { return get<1>(E[off + i]); });

					vector<dygrl::CompressedWalks> vec_compwalks; vec_compwalks.push_back(dygrl::CompressedWalks(batch_num));
                    new_verts[i] = make_pair(v, VertexEntry(types::CompressedEdges(S, v, fl), vec_compwalks, new dygrl::SamplerManager(0)));
                });

                types::MapAffectedVertices rewalk_points = types::MapAffectedVertices();

                auto replace = [&, run_seq] (const intV& v, const VertexEntry& a, const VertexEntry& b)
                {
                    auto union_edge_tree = edge_plus::uniont(b.compressed_edges, a.compressed_edges, v, run_seq);

                    lists::deallocate(a.compressed_edges.plus);
                    edge_plus::Tree_GC::decrement_recursive(a.compressed_edges.root, run_seq);

                    lists::deallocate(b.compressed_edges.plus);
                    edge_plus::Tree_GC::decrement_recursive(b.compressed_edges.root, run_seq);

					MAV_time.start();
                    auto triplets_to_delete_pbbs   = pbbs::new_array<std::vector<types::PairedTriplet>>(a.compressed_walks.size());
                    auto triplets_to_delete_vector = std::vector<std::vector<types::PairedTriplet>>();
					parallel_for(0, a.compressed_walks.size(), [&](int index)
					{
					    a.compressed_walks[index].iter_elms(v, [&](auto value)
						{
							auto pair = pairings::Szudzik<types::PairedTriplet>::unpair(value);
							auto walk_id = pair.first / config::walk_length;
							auto position = pair.first - (walk_id * config::walk_length);
							auto next = pair.second;


							auto p_min_global = config::walk_length;
							// find the p_min_global among all existing MAVS for w
							for (auto mav = a.compressed_walks[index].created_at_batch+1; mav < batch_num; mav++)
							{
								if (MAVS2[mav].contains(walk_id))
								{
									auto temp_pos = get<0>((MAVS2[mav]).find(walk_id));
									if (temp_pos < p_min_global)
										p_min_global = temp_pos; // TODO: an accumulated MAV with p_min up to that point might suffice
								}
							} // constructed the p_min_global for this w. preffix of MAVS. preffix tree (trie data structure?)

							// Check the relationship of the triplet with respect to the p_min_global or the w
							if (position < p_min_global)
							{
								// take the triplet under consideration for the MAV and proceed normally
								if (!rewalk_points.contains(walk_id))
								{
									rewalk_points.insert(walk_id, std::make_tuple(position, v, false));
								}
								else
								{
									types::Position current_min_pos = get<0>(rewalk_points.find(walk_id));

									if (current_min_pos > position)
									{
										rewalk_points.update(walk_id, std::make_tuple(position, v, false));
									}
								}
							}
							else
							{
								triplets_to_delete_pbbs[index].push_back(value);
							}

						});
					});

					// Create a new vector of compressed walks
					vector<dygrl::CompressedWalks> vec_compwalks;
					for (auto j = 0; j < a.compressed_walks.size(); j++)
					{
						auto sequence = pbbs::sequence<types::Vertex>(triplets_to_delete_pbbs[j].size());
						parallel_for(0, triplets_to_delete_pbbs[j].size(), [&](auto k) {
					        sequence[k] = triplets_to_delete_pbbs[j][k];
					    });
					    pbbs::sample_sort_inplace(pbbs::make_range(sequence.begin(), sequence.end()), std::less<>());

					    vec_compwalks.push_back(dygrl::CompressedWalks(sequence, v, 666, 666, batch_num-1));
					}

					// Do the differences
					std::vector<dygrl::CompressedWalks> new_compressed_vector;
					for (auto ind = 0; ind < a.compressed_walks.size(); ind++)
					{
				        auto refined_walk_tree = walk_plus::difference(vec_compwalks[ind], a.compressed_walks[ind], v);
						new_compressed_vector.push_back(dygrl::CompressedWalks(refined_walk_tree.plus, refined_walk_tree.root, 666, 666, batch_num-1));

					    // deallocate the memory
                        lists::deallocate(vec_compwalks[ind].plus);
                        walk_plus::Tree_GC::decrement_recursive(vec_compwalks[ind].root);
                        lists::deallocate(a.compressed_walks[ind].plus);
                        walk_plus::Tree_GC::decrement_recursive(a.compressed_walks[ind].root);
					}

					// Merge all the end trees into one
	                std::vector<dygrl::CompressedWalks> final_compressed_vector;
	                final_compressed_vector.push_back(CompressedWalks(batch_num-1));
	                for (auto ind = 0; ind < new_compressed_vector.size(); ind++)
	                {
	                    auto union_all_tree = walk_plus::uniont(new_compressed_vector[ind], final_compressed_vector[0], v);

	                    // deallocate the memory
	                    lists::deallocate(new_compressed_vector[ind].plus);
	                    walk_plus::Tree_GC::decrement_recursive(new_compressed_vector[ind].root); // deallocation is important for performance
	                    lists::deallocate(final_compressed_vector[0].plus);
	                    walk_plus::Tree_GC::decrement_recursive(final_compressed_vector[0].root);

	                    final_compressed_vector[0] = dygrl::CompressedWalks(union_all_tree.plus, union_all_tree.root, 666, 666, batch_num-1);
	                }
                    MAV_time.stop();

	                return VertexEntry(union_edge_tree, final_compressed_vector, b.sampler_manager);
                };

                graph_update_time_on_insert.start();
                this->graph_tree = Graph::Tree::multi_insert_sorted_with_values(this->graph_tree.root, new_verts, num_starts, replace,true, run_seq);
                graph_update_time_on_insert.stop();

				// Calculate the distribution of p_min in the affected walks
				uintV pmin_distribution[config::walk_length];
	            for (auto i = 0; i < config::walk_length; i++)
	                pmin_distribution[i] = 0; // initialize the distribution

				// Store/cache the MAV of each batch
	            for(auto& entry : rewalk_points.lock_table()) // todo: blocking?
	            {
		            MAVS2[batch_num].insert(entry.first, entry.second);
					// populate the distribution
					pmin_distribution[get<0>(entry.second)]++;
	            }
	            assert(rewalk_points.size() == MAVS2[batch_num].size());

				// Print out the distribution
				cout << endl << "pmin distribution" << endl;
				for (auto i = 0; i < config::walk_length; i++)
				{
					cout << pmin_distribution[i] << endl;
				}

                walk_update_time_on_insert.start();
                auto affected_walks = pbbs::sequence<types::WalkID>(rewalk_points.size());
                if (apply_walk_updates)
                        this->batch_walk_update(rewalk_points, affected_walks, batch_num);
                walk_update_time_on_insert.stop();

                // 6. Deallocate memory
                if (num_starts > stack_size) pbbs::free_array(new_verts);
                if (edges_deduped)           pbbs::free_array(edges_deduped);

                #ifdef WHARF_DEBUG
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
            pbbs::sequence<types::WalkID> delete_edges_batch(size_t m,
                                                             tuple<uintV,
                                                             uintV>* edges,
															 int batch_num,
                                                             bool sorted = false,
                                                             bool remove_dups = false,
                                                             size_t nn = std::numeric_limits<size_t>::max(),
                                                             bool apply_walk_updates = true,
                                                             bool run_seq = false)
            {
                auto fl = run_seq ? pbbs::fl_sequential : pbbs::no_flag;

                // 1. Set up
                using Edge = tuple<uintV, uintV>;

                auto edges_original = pbbs::make_range(edges, edges + m);
                Edge* edges_deduped = nullptr;

                // 2. Sort the edges in the batch (by source)
                if (!sorted)
                {
                    Wharf::sort_edge_batch_by_source(edges, m, nn);
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

                    vector<dygrl::CompressedWalks> vec_compwalks; vec_compwalks.push_back(dygrl::CompressedWalks(batch_num));
                    new_verts[i] = make_pair(v, VertexEntry(types::CompressedEdges(S, v, fl), vec_compwalks, new dygrl::SamplerManager(0)));
                });

                types::MapAffectedVertices rewalk_points = types::MapAffectedVertices();

                auto replace = [&,run_seq] (const intV& v, const VertexEntry& a, const VertexEntry& b)
                {
                    auto difference_edge_tree = edge_plus::difference(b.compressed_edges, a.compressed_edges, v, run_seq);

                    lists::deallocate(a.compressed_edges.plus);
                    edge_plus::Tree_GC::decrement_recursive(a.compressed_edges.root, run_seq);

                    lists::deallocate(b.compressed_edges.plus);
                    edge_plus::Tree_GC::decrement_recursive(b.compressed_edges.root, run_seq);

                    bool should_reset = difference_edge_tree.degree() == 0;

                    a.compressed_walks.back().iter_elms(v, [&](auto value)
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

	            for(auto& entry : rewalk_points.lock_table())
	            {
		            MAVS2[batch_num].insert(entry.first, entry.second);
	            }
	            assert(rewalk_points.size() == MAVS2[batch_num].size());

                walk_update_time_on_delete.start();
                auto affected_walks = pbbs::sequence<types::WalkID>(rewalk_points.size());
                if (apply_walk_updates) this->batch_walk_update(rewalk_points, affected_walks, batch_num);
                walk_update_time_on_delete.stop();

                // 6. Deallocate memory
                if (num_starts > stack_size) pbbs::free_array(new_verts);
                if (edges_deduped) pbbs::free_array(edges_deduped);

                #ifdef WHARF_DEBUG
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
             * @param types::MapOfChanges - rewalking points
             */
            void batch_walk_update(types::MapAffectedVertices& rewalk_points, pbbs::sequence<types::WalkID>& affected_walks, int batch_num)
            {
                types::ChangeAccumulator inserts = types::ChangeAccumulator();

                uintV index = 0;

                for(auto& entry : rewalk_points.lock_table()) // todo: blocking?
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

				if (!config::parallel_merge_wu)
				{
					// -----------------------
					// First the walk sampling
					// -----------------------
					Walking_new_sampling_time.start();
					// Parallel Update of Affected Walks
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
						  current_vertex_new_walk = current_vertex_old_walk = affected_walks[index] % this->number_of_vertices();
					  }

					  if (graph[current_vertex_new_walk].degree == 0)
					  {
						  types::PairedTriplet hash = pairings::Szudzik<types::Vertex>::pair({affected_walks[index] * config::walk_length + 0, current_vertex_new_walk});

						  if (!inserts.contains(current_vertex_new_walk)) inserts.insert(current_vertex_new_walk, std::vector<types::PairedTriplet>());
						  inserts.update_fn(current_vertex_new_walk, [&](auto& vector) {
							vector.push_back(hash);
						  });

						  return;
					  }

					  auto random = config::random;
					  if (config::deterministic_mode)
						  random = utility::Random(affected_walks[index] / this->number_of_vertices());
					  auto state = model->initial_state(current_vertex_new_walk);

					  for (types::Position position = current_position; position < config::walk_length; position++)
					  {
						  if (!graph[state.first].samplers->contains(state.second))
							  graph[state.first].samplers->insert(state.second, MetropolisHastingsSampler(state, model));

						  auto temp_state = graph[state.first].samplers->find(state.second).sample(state, model);
						  if (config::deterministic_mode)
						  {
							  state = model->new_state(state, graph[state.first].neighbors[random.irand(graph[state.first].degree)]);
						  }
						  else
							  state = temp_state;

						  number_of_sampled_vertices++;
						  types::PairedTriplet hash = (position != config::walk_length - 1) ?
						                              pairings::Szudzik<types::Vertex>::pair({affected_walks[index] * config::walk_length + position, state.first}) : // new sampled next
						                              pairings::Szudzik<types::Vertex>::pair({affected_walks[index] * config::walk_length + position, current_vertex_new_walk});

						  if (!inserts.contains(current_vertex_new_walk)) inserts.insert(current_vertex_new_walk, std::vector<types::PairedTriplet>());
						  inserts.update_fn(current_vertex_new_walk, [&](auto& vector) {
							vector.push_back(hash);
						  });

						  current_vertex_new_walk = state.first;
					  }

					});
					Walking_new_sampling_time.stop();

					// -----------------------
					// Then the merge process
					// -----------------------
					Merge_time.start();
					if (batch_num % config::merge_frequency == 0)
					{
						cout << "\n*** merge at batch-" << batch_num << endl;
						merge_walk_trees_all_vertices_parallel(batch_num); // ok merge all old nodes
					}
					Merge_time.stop();
				}
				else // parallel resample and merge
				{
					#pragma omp parallel
					{
						#pragma omp task
						{
							Walking_new_sampling_time.start();
							// Parallel Update of Affected Walks
							parallel_for(0, affected_walks.size(), [&](auto index)
							{
							  auto entry = rewalk_points.template find(affected_walks[index]);

							  auto current_position        = std::get<0>(entry);
							  auto current_vertex_old_walk = std::get<1>(entry);
							  auto should_reset            = std::get<2>(entry);

							  auto current_vertex_new_walk = current_vertex_old_walk;

							  if (should_reset) // todo: clear out if this is needed
							  {
								  current_position = 0;
								  current_vertex_new_walk = current_vertex_old_walk = affected_walks[index] % this->number_of_vertices();
							  }

							  if (graph[current_vertex_new_walk].degree == 0)
							  {
								  types::PairedTriplet hash = pairings::Szudzik<types::Vertex>::pair({affected_walks[index] * config::walk_length + 0, current_vertex_new_walk});

								  if (!inserts.contains(current_vertex_new_walk)) inserts.insert(current_vertex_new_walk, std::vector<types::PairedTriplet>());
								  inserts.update_fn(current_vertex_new_walk, [&](auto& vector) {
									vector.push_back(hash);
								  });

								  return;
							  }

							  auto random = config::random;
							  if (config::deterministic_mode)
								  random = utility::Random(affected_walks[index] / this->number_of_vertices());
							  auto state = model->initial_state(current_vertex_new_walk);

							  for (types::Position position = current_position; position < config::walk_length; position++)
							  {
								  if (!graph[state.first].samplers->contains(state.second))
									  graph[state.first].samplers->insert(state.second, MetropolisHastingsSampler(state, model));

								  auto temp_state = graph[state.first].samplers->find(state.second).sample(state, model);
								  if (config::deterministic_mode)
								  {
									  state = model->new_state(state, graph[state.first].neighbors[random.irand(graph[state.first].degree)]);
								  }
								  else
									  state = temp_state;

								  number_of_sampled_vertices++;

								  types::PairedTriplet hash = (position != config::walk_length - 1) ?
								                              pairings::Szudzik<types::Vertex>::pair({affected_walks[index] * config::walk_length + position, state.first}) : // new sampled next
								                              pairings::Szudzik<types::Vertex>::pair({affected_walks[index] * config::walk_length + position, current_vertex_new_walk});

								  if (!inserts.contains(current_vertex_new_walk)) inserts.insert(current_vertex_new_walk, std::vector<types::PairedTriplet>());
								  inserts.update_fn(current_vertex_new_walk, [&](auto& vector) {
									vector.push_back(hash);
								  });

								  // Then, change the current vertex in the new walk
								  current_vertex_new_walk = state.first;
							  }

							});
							Walking_new_sampling_time.stop();
						};

						#pragma omp task
						{
							// ------------
				            // --- Merge non-MAV vertices
				            // ------------
				            Merge_time.start();
						    if (batch_num % config::merge_frequency == 0)
						    {
				                merge_walk_trees_all_vertices_parallel(batch_num);
							}
				            Merge_time.stop();
						};
					}
				}

	            Walking_insert_new_samples.start();
	            using VertexStruct  = std::pair<types::Vertex, VertexEntry>;
                auto insert_walks  = pbbs::sequence<VertexStruct>(inserts.size());

				// Iterate and put the Insert accumulator into a pbbs:sequence
				auto temp_inserts = pbbs::sequence<std::pair<types::Vertex, std::vector<types::PairedTriplet>>>(inserts.size());
				auto skatindex = 0;
				for (auto& item: inserts.lock_table()) // TODO: This cannot be in parallel. Cuckoo-hashmap limitation
				{
					temp_inserts[skatindex++] = std::make_pair(item.first, item.second);
				}

				parallel_for(0, temp_inserts.size(), [&](auto j){
                    auto sequence = pbbs::sequence<types::Vertex>(temp_inserts[j].second.size());

	                parallel_for(0, temp_inserts[j].second.size(), [&](auto i){
						sequence[i] = temp_inserts[j].second[i];
					});

                    pbbs::sample_sort_inplace(pbbs::make_range(sequence.begin(), sequence.end()), std::less<>());
					vector<dygrl::CompressedWalks> vec_compwalks;
					vec_compwalks.push_back(dygrl::CompressedWalks(sequence, temp_inserts[j].first, 666, 666, batch_num));
                    insert_walks[j] = std::make_pair(temp_inserts[j].first, VertexEntry(types::CompressedEdges(), vec_compwalks, new dygrl::SamplerManager(0)));
                });

                pbbs::sample_sort_inplace(pbbs::make_range(insert_walks.begin(), insert_walks.end()), [&](auto& x, auto& y) {
                    return x.first < y.first;
                });

                auto replaceI = [&] (const uintV& src, const VertexEntry& x, const VertexEntry& y)
                {
					// add compressed walks of y to x
					auto x_prime = x.compressed_walks;
					if (y.compressed_walks.back().size() != 0)
						x_prime.push_back(y.compressed_walks.back()); // y has only one walk tree

                    return VertexEntry(x.compressed_edges, x_prime, x.sampler_manager);
                };

				cout << "\n(insert) -- For batch-" << batch_num << " we are touching " << insert_walks.size() << " / " << number_of_vertices() << " vertices" << endl;

                this->graph_tree = Graph::Tree::multi_insert_sorted_with_values(this->graph_tree.root, insert_walks.begin(), insert_walks.size(), replaceI, true);
				Walking_insert_new_samples.stop();
            } // End of batch walk update procedure


			/**
			 * @brief Merges the walk-trees of each vertex in the hybrid-tree such that in the end each vertex has only one walk-tree
			 */
			 void merge_walk_trees_all_vertices_parallel(int num_batches_so_far)
			 {
			    libcuckoo::cuckoohash_map<types::Vertex, std::vector<std::vector<types::PairedTriplet>>> all_to_delete;

	            auto flat_graph = this->flatten_vertex_tree();

				merge_calc_triplets_to_delete.start();
				parallel_for(0, this->number_of_vertices(), [&](size_t i)
				{
					int inc = 0;

					auto triplets_to_delete_pbbs   = pbbs::new_array<std::vector<types::PairedTriplet>>(flat_graph[i].compressed_walks.size());
					auto triplets_to_delete_vector = std::vector<std::vector<types::PairedTriplet>>();

					// traverse each walk-tree and find out the obsolete triplets and create corresponding "deletion" walk-trees
					for (auto wt = flat_graph[i].compressed_walks.begin(); wt != flat_graph[i].compressed_walks.end(); wt++)
					{
						// Define the triplets to delete vector for each walk-tree
						wt->iter_elms(i, [&](auto enc_triplet)
						{
						  auto pair = pairings::Szudzik<types::Vertex>::unpair(enc_triplet);

						  auto walk_id  = pair.first / config::walk_length;
						  auto position = pair.first - (walk_id * config::walk_length);
						  auto next_vertex   = pair.second;

			              auto p_min_global = config::walk_length;
						  for (auto mav = wt->created_at_batch+1; mav < num_batches_so_far+1; mav++)
					      {
   							  if (MAVS2[mav].template contains(walk_id))
							  {
							       auto temp_pos = get<0>((MAVS2[mav]).template find(walk_id)); // it does not always contain this wid
								   if (temp_pos < p_min_global)
									   p_min_global = temp_pos;
							  }
						  } // constructed the p_min_global for this w. preffix of MAVS. preffix tree (trie data structure?)

						  // Check the relationship of the triplet with respect to the p_min_global or the w
						  if (position < p_min_global)
						  {
							; // the triplet is still valid so it stays
						  }
						  else
						  {
							triplets_to_delete_pbbs[inc].push_back(enc_triplet);
						  }

						});

						// pass it to the vector
						triplets_to_delete_vector.push_back(triplets_to_delete_pbbs[inc]);
						inc++;
					}

					// add the triplets to delete for this vertex in a hashmap
                    if (!all_to_delete.contains(i))
						all_to_delete.insert(i, std::vector<std::vector<types::PairedTriplet>>());
                    all_to_delete.update_fn(i, [&](auto& vector) {
                        vector = triplets_to_delete_vector;
                    });
				});
				merge_calc_triplets_to_delete.stop();

				merge_create_delete_walks.start();
				auto temp_deletes = pbbs::sequence<std::pair<types::Vertex, std::vector<std::vector<types::PairedTriplet>>>>(all_to_delete.size());
				auto skatindex = 0;

				for (auto& item: all_to_delete.lock_table())
				{
				  temp_deletes[skatindex++] = std::make_pair(item.first, item.second);
				}

				using VertexStruct = std::pair<types::Vertex, VertexEntry>;
				auto delete_walks  = pbbs::sequence<VertexStruct>(all_to_delete.size());
				parallel_for(0, temp_deletes.size(), [&](auto kkk)
				{
					vector<dygrl::CompressedWalks> vec_compwalks;

					auto vertex_id = temp_deletes[kkk].first;

					for (auto j = 0; j < temp_deletes[kkk].second.size(); j++)
					{
						auto sequence = pbbs::sequence<types::Vertex>(temp_deletes[kkk].second[j].size());
						parallel_for(0, temp_deletes[kkk].second[j].size(), [&](auto k){
					        sequence[k] = temp_deletes[kkk].second[j][k];
						});
						pbbs::sample_sort_inplace(pbbs::make_range(sequence.begin(), sequence.end()), std::less<>());

						vec_compwalks.push_back(dygrl::CompressedWalks(sequence, temp_deletes[kkk].first, 666, 666, num_batches_so_far)); // dummy min,max, batch_num
					}

					delete_walks[kkk] = std::make_pair(temp_deletes[kkk].first, VertexEntry(types::CompressedEdges(), vec_compwalks, new dygrl::SamplerManager(0)));
				});
				merge_create_delete_walks.stop();

				// Sort the delete walks
				pbbs::sample_sort_inplace(pbbs::make_range(delete_walks.begin(), delete_walks.end()), [&](auto& x, auto& y) {
				  return x.first < y.first;
				});

				auto replaceI = [&] (const uintV& src, const VertexEntry& x, const VertexEntry& y)
				{
					assert(x.compressed_walks.size() == y.compressed_walks.size());
				    std::vector<dygrl::CompressedWalks> new_compressed_vector;
				    for (auto ind = 0; ind < x.compressed_walks.size(); ind++)
				    {
				        auto refined_walk_tree = walk_plus::difference(y.compressed_walks[ind], x.compressed_walks[ind], src);
					    new_compressed_vector.push_back(dygrl::CompressedWalks(refined_walk_tree.plus, refined_walk_tree.root, 666, 666, num_batches_so_far));

					    // deallocate the memory
                        lists::deallocate(x.compressed_walks[ind].plus);
		                walk_plus::Tree_GC::decrement_recursive(x.compressed_walks[ind].root);
		                lists::deallocate(y.compressed_walks[ind].plus);
		                walk_plus::Tree_GC::decrement_recursive(y.compressed_walks[ind].root);
				    }

				    // merge the refined walk-trees here
				    std::vector<dygrl::CompressedWalks> final_compressed_vector;

				    final_compressed_vector.push_back(CompressedWalks(num_batches_so_far));
				    for (auto ind = 0; ind < new_compressed_vector.size(); ind++)
				    {
		                auto union_all_tree = walk_plus::uniont(new_compressed_vector[ind], final_compressed_vector[0], src);

					    // deallocate the memory
					    lists::deallocate(new_compressed_vector[ind].plus);
					    walk_plus::Tree_GC::decrement_recursive(new_compressed_vector[ind].root);
					    lists::deallocate(final_compressed_vector[0].plus);
					    walk_plus::Tree_GC::decrement_recursive(final_compressed_vector[0].root);

					    final_compressed_vector[0] = dygrl::CompressedWalks(union_all_tree.plus, union_all_tree.root, 666, 666, num_batches_so_far);
				    }

                    return VertexEntry(x.compressed_edges, final_compressed_vector, x.sampler_manager);
				};

				cout << "\n(merge) -- For batch-" << num_batches_so_far << " we are touching " << delete_walks.size() << " / " << number_of_vertices() << " vertices" << endl;

				merge_multiinsert_ctress.start();
				this->graph_tree = Graph::Tree::multi_insert_sorted_with_values(this->graph_tree.root, delete_walks.begin(), delete_walks.size(), replaceI, true);
				merge_multiinsert_ctress.stop();
			}

			/**
			 * @brief Merges the walk-trees of each vertex in the hybrid-tree such that in the end each vertex has only one walk-tree
			 */
			 void last_merge_all_vertices_parallel_with_minmax(int num_batches_so_far)
			 {
			    libcuckoo::cuckoohash_map<types::Vertex, std::vector<std::vector<types::PairedTriplet>>> all_to_delete;
				auto next_minmax = pbbs::sequence<std::pair<types::Vertex, types::Vertex>>(this->number_of_vertices());

	            auto flat_graph = this->flatten_vertex_tree();

				merge_calc_triplets_to_delete.start();
				parallel_for(0, this->number_of_vertices(), [&](size_t i)
				{
					// Initialization of next-{min, max}
					types::Vertex next_min = UINT32_MAX;
					types::Vertex next_max = 0;

					int inc = 0;

					auto triplets_to_delete_pbbs   = pbbs::new_array<std::vector<types::PairedTriplet>>(flat_graph[i].compressed_walks.size());
					auto triplets_to_delete_vector = std::vector<std::vector<types::PairedTriplet>>();

					// traverse each walk-tree and find out the obsolete triplets and create corresponding "deletion" walk-trees
					for (auto wt = flat_graph[i].compressed_walks.begin(); wt != flat_graph[i].compressed_walks.end(); wt++)
					{
						// Define the triplets to delete vector for each walk-tree
						wt->iter_elms(i, [&](auto enc_triplet)
						{
						  auto pair = pairings::Szudzik<types::Vertex>::unpair(enc_triplet);

						  auto walk_id  = pair.first / config::walk_length;
						  auto position = pair.first - (walk_id * config::walk_length);
						  auto next_vertex   = pair.second;

			              auto p_min_global = config::walk_length;
						  for (auto mav = wt->created_at_batch+1; mav < num_batches_so_far+1; mav++)
					      {
   							  if (MAVS2[mav].template contains(walk_id))
							  {
							       auto temp_pos = get<0>((MAVS2[mav]).template find(walk_id)); // it does not always contain this wid
								   if (temp_pos < p_min_global)
									   p_min_global = temp_pos;
							  }
						  } // constructed the p_min_global for this w. preffix of MAVS.

						  // Check the relationship of the triplet with respect to the p_min_global or the w
						  if (position < p_min_global) // TODO: this accepts all?
						  {
//							; // the triplet is still valid so it stays
							next_min = std::min(next_min, next_vertex);
							next_max = std::max(next_max, next_vertex);
						  }
						  else
						  {
							triplets_to_delete_pbbs[inc].push_back(enc_triplet);
						  }

						});

						// pass it to the vector
						triplets_to_delete_vector.push_back(triplets_to_delete_pbbs[inc]);
						inc++;
					}

					// add the triplets to delete for this vertex in a hashmap
                    if (!all_to_delete.contains(i))
						all_to_delete.insert(i, std::vector<std::vector<types::PairedTriplet>>());
                    all_to_delete.update_fn(i, [&](auto& vector) {
                        vector = triplets_to_delete_vector;
                    });

				    // Cache the minmax
			        next_minmax[i] = std::make_pair(next_min, next_max);

				});
				merge_calc_triplets_to_delete.stop();

				merge_create_delete_walks.start();
				auto temp_deletes = pbbs::sequence<std::pair<types::Vertex, std::vector<std::vector<types::PairedTriplet>>>>(all_to_delete.size());
				auto skatindex = 0;

				for (auto& item: all_to_delete.lock_table())
				{
				  temp_deletes[skatindex++] = std::make_pair(item.first, item.second);
				}

				using VertexStruct = std::pair<types::Vertex, VertexEntry>;
				auto delete_walks  = pbbs::sequence<VertexStruct>(all_to_delete.size());
				cout << "all_to_delete size: " << all_to_delete.size() << endl;
				parallel_for(0, temp_deletes.size(), [&](auto kkk)
				{
					vector<dygrl::CompressedWalks> vec_compwalks;

					auto vertex_id = temp_deletes[kkk].first;

					for (auto j = 0; j < temp_deletes[kkk].second.size(); j++)
					{
						auto sequence = pbbs::sequence<types::Vertex>(temp_deletes[kkk].second[j].size());
						parallel_for(0, temp_deletes[kkk].second[j].size(), [&](auto k){
					        sequence[k] = temp_deletes[kkk].second[j][k];
						});
						pbbs::sample_sort_inplace(pbbs::make_range(sequence.begin(), sequence.end()), std::less<>());

						vec_compwalks.push_back(dygrl::CompressedWalks(sequence, temp_deletes[kkk].first, 666, 666, num_batches_so_far)); // dummy min,max, batch_num
					}

					delete_walks[kkk] = std::make_pair(temp_deletes[kkk].first, VertexEntry(types::CompressedEdges(), vec_compwalks, new dygrl::SamplerManager(0)));
				});
				merge_create_delete_walks.stop();

				// Sort the delete walks
				pbbs::sample_sort_inplace(pbbs::make_range(delete_walks.begin(), delete_walks.end()), [&](auto& x, auto& y) {
				  return x.first < y.first;
				});

				auto replaceI = [&] (const uintV& src, const VertexEntry& x, const VertexEntry& y)
				{
					assert(x.compressed_walks.size() == y.compressed_walks.size());
				    std::vector<dygrl::CompressedWalks> new_compressed_vector;
				    for (auto ind = 0; ind < x.compressed_walks.size(); ind++)
				    {
				        auto refined_walk_tree = walk_plus::difference(y.compressed_walks[ind], x.compressed_walks[ind], src);
					    new_compressed_vector.push_back(dygrl::CompressedWalks(refined_walk_tree.plus, refined_walk_tree.root, 666, 666, num_batches_so_far)); // use dummy min, max, batch_num for now

					    // deallocate the memory
                        lists::deallocate(x.compressed_walks[ind].plus);
		                walk_plus::Tree_GC::decrement_recursive(x.compressed_walks[ind].root);
		                lists::deallocate(y.compressed_walks[ind].plus);
		                walk_plus::Tree_GC::decrement_recursive(y.compressed_walks[ind].root);
				    }

				    // merge the refined walk-trees here
				    std::vector<dygrl::CompressedWalks> final_compressed_vector;

				    final_compressed_vector.push_back(CompressedWalks(num_batches_so_far));
				    for (auto ind = 0; ind < new_compressed_vector.size(); ind++)
				    {
		                auto union_all_tree = walk_plus::uniont(new_compressed_vector[ind], final_compressed_vector[0], src);

					    // deallocate the memory
					    lists::deallocate(new_compressed_vector[ind].plus);
					    walk_plus::Tree_GC::decrement_recursive(new_compressed_vector[ind].root);
					    lists::deallocate(final_compressed_vector[0].plus);
					    walk_plus::Tree_GC::decrement_recursive(final_compressed_vector[0].root);

					    final_compressed_vector[0] = dygrl::CompressedWalks(union_all_tree.plus, union_all_tree.root, 666, 666, num_batches_so_far);
				    }

					// Replace the next min-max correctly
					final_compressed_vector[0].vnext_min = next_minmax[src].first;
					final_compressed_vector[0].vnext_max = next_minmax[src].second;

                    return VertexEntry(x.compressed_edges, final_compressed_vector, x.sampler_manager);
				};

				cout << "\n(Last Merge) -- we are touching " << delete_walks.size() << " / " << number_of_vertices() << " vertices" << endl;

				merge_multiinsert_ctress.start();
				this->graph_tree = Graph::Tree::multi_insert_sorted_with_values(this->graph_tree.root, delete_walks.begin(), delete_walks.size(), replaceI, true);
				merge_multiinsert_ctress.stop();
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

                // measure size of {min,max} bounds for each walk-tree
                size_t vnext_min_bytes  = 0;
                size_t vnext_max_bytes  = 0;

                for (auto i = 0; i < flat_graph.size(); i++)
                {
                    flat_graph_bytes += sizeof(flat_graph[i]);

                    edges_heads += flat_graph[i].compressed_edges.edge_tree_nodes();
					for (auto j = flat_graph[i].compressed_walks.rbegin(); j != flat_graph[i].compressed_walks.rend(); j++)
						walks_heads += j->edge_tree_nodes();

                    edges_bytes += flat_graph[i].compressed_edges.size_in_bytes(i);
	                for (auto j = flat_graph[i].compressed_walks.rbegin(); j != flat_graph[i].compressed_walks.rend(); j++)
                        walks_bytes += j->size_in_bytes(i);

                    for (auto& entry : flat_graph[i].sampler_manager->lock_table())
                    {
                        samplers_bytes += sizeof(entry.first) + sizeof(entry.second);
                    }

	                // measure the size of all {min, max} bounds of range search
	                vnext_min_bytes   += sizeof(flat_graph[i].compressed_walks.front().vnext_min);
	                vnext_max_bytes   += sizeof(flat_graph[i].compressed_walks.front().vnext_max);
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

                std::cout << "Range Search (" << ((config::range_search_mode) ? "ON" : "OFF") << "): \n\t"
                          << "Total {min,max} memory usage: " << utility::MB(((config::range_search_mode) ? (vnext_min_bytes + vnext_max_bytes) : 0))
                          << " MB = " << utility::GB(((config::range_search_mode) ? (vnext_min_bytes + vnext_max_bytes) : 0))
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
                        + edges_bytes + edges_heads*c_tree_node_size + samplers_bytes + ((config::range_search_mode) ? (vnext_min_bytes + vnext_max_bytes) : 0);

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
			std::vector<types::MapAffectedVertices> MAVS2;

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
                #ifdef WHARF_TIMER
                    timer timer("Wharf::SortEdgeBatchBySource");
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

                #ifdef WHARF_TIMER
                    timer.reportTotal("time (seconds)");
                #endif
            }
    };
}

#endif