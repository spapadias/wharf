#include <wharf.h>

void throughput(commandLine& command_line)
{
	string fname             = string(command_line.getOptionValue("-f", default_file_name));
	bool mmap                = command_line.getOption("-m");
	bool is_symmetric        = command_line.getOption("-s");
	bool compressed          = command_line.getOption("-c");
	size_t n_parts           = command_line.getOptionLongValue("-nparts", 1);

	size_t walks_per_vertex = command_line.getOptionLongValue("-w", config::walks_per_vertex);
	size_t length_of_walks  = command_line.getOptionLongValue("-l", config::walk_length);
	string model            = string(command_line.getOptionValue("-model", "deepwalk"));
	double paramP           = command_line.getOptionDoubleValue("-p", config::paramP);
	double paramQ           = command_line.getOptionDoubleValue("-q", config::paramQ);
	string init_strategy    = string(command_line.getOptionValue("-init", "weight"));
	size_t n_trials         = command_line.getOptionLongValue("-trials",  5); //3); todo: changed the number of trials

	string determinism      = string(command_line.getOptionValue("-det", "false"));
	string range_search     = string(command_line.getOptionValue("-rs", "true"));
	size_t num_of_batches   = command_line.getOptionIntValue("-nb", 10);
	size_t half_of_bsize    = command_line.getOptionIntValue("-bs", 5000);
	size_t merge_frequency  = command_line.getOptionIntValue("-mergefreq", 1);
	config::merge_frequency = merge_frequency;
	string mergemode        = string(command_line.getOptionValue("-mergemode", "parallel"));
	if (mergemode == "parallel")
	{
		config::parallel_merge_wu = true;
		std::cout << "Parallel Merge and WU" << std::endl;
	}
	else
	{
		config::parallel_merge_wu = false;
		std::cout << "Serial Merge and WU" << std::endl;
	}

	config::walks_per_vertex = walks_per_vertex;
	config::walk_length      = length_of_walks;

	std::cout << "Walks per vertex: " << (int) config::walks_per_vertex << std::endl;
	std::cout << "Walk length: " <<  (int) config::walk_length << std::endl;

	if (model == "deepwalk")
	{
		config::random_walk_model = types::RandomWalkModelType::DEEPWALK;
		std::cout << "Walking model: DEEPWALK" << std::endl;
	}
	else if (model == "node2vec")
	{
		config::random_walk_model = types::RandomWalkModelType::NODE2VEC;
		config::paramP = paramP;
		config::paramQ = paramQ;
		std::cout << "Walking model: NODE2VEC | Params (p,q) = " << "(" << config::paramP << "," << config::paramQ << ")" << std::endl;
	}
	else
	{
		std::cerr << "Unrecognized walking model! Abort" << std::endl;
		std::exit(1);
	}

	if (init_strategy == "burnin")
	{
		config::sampler_init_strategy = types::SamplerInitStartegy::BURNIN;
		std::cout << "Sampler strategy: BURNIN" << std::endl;
	}
	else if (init_strategy == "weight")
	{
		config::sampler_init_strategy = types::SamplerInitStartegy::WEIGHT;
		std::cout << "Sampler strategy: WEIGHT" << std::endl;
	}
	else if (init_strategy == "random")
	{
		config::sampler_init_strategy = types::SamplerInitStartegy::RANDOM;
		std::cout << "Sampler strategy: RANDOM" << std::endl;
	}
	else
	{
		std::cerr << "Unrecognized sampler init strategy" << std::endl;
		std::exit(1);
	}

	// Set up the range search mode ------
	if (range_search == "true")
		config::range_search_mode = true;
	else
		config::range_search_mode = false;

	// Set up the deterministic mode
	if (determinism == "true")
		config::deterministic_mode = true;
	else
		config::deterministic_mode = false;
	// ------------------------------------

	size_t n;
	size_t m;
	uintE* offsets;
	uintV* edges;
	std::tie(n, m, offsets, edges) = read_unweighted_graph(fname.c_str(), is_symmetric, mmap);

	dygrl::Wharf malin = dygrl::Wharf(n, m, offsets, edges);

	timer initial_walks_from_scratch_timer("GenerateInitialWalksFromScratch", false);
	initial_walks_from_scratch_timer.start();
	malin.generate_initial_random_walks();
	initial_walks_from_scratch_timer.stop();
	cout << "Total time to generate the walk corpus from scratch: " << initial_walks_from_scratch_timer.get_total() << endl;

	int n_batches = num_of_batches; // todo: how many batches per batch size?

	// TODO: Why incorrect numbers when MALIN_DEBUG is off?

	auto batch_sizes = pbbs::sequence<size_t>(1);
	batch_sizes[0] = half_of_bsize; //5;
//	batch_sizes[1] = 50;
//	batch_sizes[2] = 500;
//	batch_sizes[3] = 5000;
//	batch_sizes[4] = 10000;
//	batch_sizes[5] = 15000;
//	batch_sizes[6] = 25000;
//	batch_sizes[7] = 50000;
//  batch_sizes[5] = 500000;


	for (short int i = 0; i < batch_sizes.size(); i++)
	{
		timer insert_timer("InsertTimer");
		timer delete_timer("DeleteTimer");

		graph_update_time_on_insert.reset();
		walk_update_time_on_insert.reset();
		// --- profiling initialization
//		walk_insert_init.reset();
		Walking_new_sampling_time.reset();
		Walking_insert_new_samples.reset();
//		ij.reset();
//		dj.reset();
//		walk_find_in_vertex_tree.reset();
//		walk_find_next_tree.reset();
//		szudzik_hash.reset();
//		fnir_tree_search.reset();
		MAV_time.reset();
		MAV_min = 20000;
		MAV_max = 0;
//		read_access_MAV.reset();
//		bdown_create_vertex_entries.reset();
//		apply_multiinsert_ctrees.reset();
//		linear_cuckoo_acc_scann.reset();
		merge_calc_triplets_to_delete.reset();
		merge_create_delete_walks.reset();
//		ij_sampling.reset();
//		ij_szudzik.reset();
		Merge_time.reset();
		Merge_min = 20000;
		Merge_max = 0;
//		sortAtMergeAll.reset();
//		accumultinsert.reset();
		LastMerge.reset();
//		FindPreviousVertexNode2Vec.reset();
//		FindNextNodeRangeSearch.reset();
//		FindNextNodeSimpleSearch.reset();
		// ---

		std::cout << "Batch size = " << 2 * batch_sizes[i] << " | ";

		double last_insert_time = 0;
		double last_MAV_time    = 0.0;
		double last_Merge_time  = 0.0;
		double last_Walk_sampling_time = 0.0;
		double last_Walk_new_insert_time = 0.0;

		auto latency_insert = pbbs::sequence<double>(n_batches);
		auto latency = pbbs::sequence<double>(n_batches);

		double total_insert_walks_affected = 0;
		double total_delete_walks_affected = 0;

		int batch_seed[n_batches];
		for (auto i = 0; i < n_batches; i++)
			batch_seed[i] = i; // say the seed equals to the #batch

		for (short int b = 0; b < n_batches; b++)
		{
			cout << "batch-" << b << " and batch_seed-" << batch_seed[b] << endl;

			if (config::random_walk_model == types::NODE2VEC)
			{
				// ---------------------------
				cout << "START -- merge before to execute node2vec" << endl;
				malin.last_merge_all_vertices_parallel_with_minmax(b+1);
				cout << "END   -- merge before to execute node2vec" << endl;
				// ---------------------------
			}

			size_t graph_size_pow2 = 1 << (pbbs::log2_up(n) - 1);
			auto edges = utility::generate_batch_of_edges(batch_sizes[i], n, batch_seed[b], false, false);

			std::cout << edges.second << " ";

			insert_timer.start();
			auto x = malin.insert_edges_batch(edges.second, edges.first, b+1, false, true, graph_size_pow2); // pass the batch number as well
			insert_timer.stop();

			total_insert_walks_affected += x.size();

			last_insert_time = walk_update_time_on_insert.get_total() - last_insert_time;
			latency_insert[b] = (double) last_insert_time / x.size();

			latency[b] = latency_insert[b];

			// Delete the batch of edges
			malin.delete_edges_batch(edges.second, edges.first, b+1, false, true, graph_size_pow2, false);

			// Update the MAV min and max
			last_MAV_time = MAV_time.get_total() - last_MAV_time;
			MAV_min = std::min(MAV_min, last_MAV_time);
			MAV_max = std::max(MAV_max, last_MAV_time);
			// Update the Merge min and max
			if ((b+1) % config::merge_frequency == 0) {
				last_Merge_time = Merge_time.get_total() - last_Merge_time;
				Merge_min = std::min(Merge_min, last_Merge_time);
				Merge_max = std::max(Merge_max, last_Merge_time);
			}
			// Update the Walk Update (without Merge) min and max
			last_Walk_sampling_time = Walking_new_sampling_time.get_total() - last_Walk_sampling_time;
			WalkSampling_min = std::min(WalkSampling_min, last_Walk_sampling_time);
			WalkSampling_max = std::max(WalkSampling_max, last_Walk_sampling_time);
			last_Walk_new_insert_time = Walking_insert_new_samples.get_total() - last_Walk_new_insert_time;
			WalkInsert_min = std::min(WalkInsert_min, last_Walk_new_insert_time);
			WalkInsert_max = std::max(WalkInsert_max, last_Walk_new_insert_time);

			pbbs::free_array(edges.first);

			cout << "METRICS AT BATCH-" << b+1 << endl;
			std::cout << "Insert time (avg) = " << insert_timer.get_total() / (b+1) << std::endl;
			std::cout << "GUP (avg) = " << graph_update_time_on_insert.get_total() / (b+1) << std::endl;
			std::cout << "BWUP (avg, includes merge) = " << walk_update_time_on_insert.get_total() / (b+1) << ", average walk affected = " << total_insert_walks_affected / (b+1) << ", sampled vertices = " << malin.number_of_sampled_vertices << std::endl;
			std::cout << "WUP (avg)   = " << (Walking_new_sampling_time.get_total() + Walking_insert_new_samples.get_total()) / (b+1) << "\t(sampling= " << Walking_new_sampling_time.get_total() / (b+1) << ", inserting= " << Walking_insert_new_samples.get_total() / (b+1) << ")" << endl;
			std::cout << "MAV (avg)   = " << MAV_time.get_total() / (b+1) << "\tMAV (min) = " << MAV_min << "\tMAV (max) = " << MAV_max << std::endl;
			std::cout << "Merge (avg," << std::floor(n_batches / merge_frequency) << " times) = " << Merge_time.get_total() / std::floor((b+1) / merge_frequency) << "\tMerge (min) = " << Merge_min << "\tMerge (max) = " << Merge_max << std::endl;
		}
		cout << fixed;
		std::cout << std::endl;

		cout << "METRICS FOR ALL BATCHES" << endl;
		std::cout << "Insert time (avg) = " << insert_timer.get_total() / n_batches << std::endl;
		std::cout << "GUP (avg) = " << graph_update_time_on_insert.get_total() / n_batches << std::endl;
		std::cout << "BWUP (avg, includes merge) = " << walk_update_time_on_insert.get_total() / n_batches << ", average walk affected = " << total_insert_walks_affected / n_batches << ", sampled vertices = " << malin.number_of_sampled_vertices << std::endl;
		std::cout << "WUP (avg)   = " << (Walking_new_sampling_time.get_total() + Walking_insert_new_samples.get_total()) / n_batches << "\t(sampling= " << Walking_new_sampling_time.get_total() / n_batches << ", inserting= " << Walking_insert_new_samples.get_total() / n_batches << ")" << endl;
		std::cout << "MAV (avg)   = " << MAV_time.get_total() / n_batches << "\tMAV (min) = " << MAV_min << "\tMAV (max) = " << MAV_max << std::endl;
		std::cout << "Merge (avg," << std::floor(n_batches / merge_frequency) << " times) = " << Merge_time.get_total() / std::floor(n_batches / merge_frequency) << "\tMerge (min) = " << Merge_min << "\tMerge (max) = " << Merge_max << std::endl;

		// latencies
		double average_latency = 0.0;
		std::cout << "Average walk insert latency = { ";
		for (int i = 0; i < n_batches; i++) {
			average_latency += latency_insert[i];
			std::cout << latency_insert[i] << " ";
		}
		std::cout << "}" << std::endl;

		cout << "(1) throughput: " << fixed << setprecision(8) << total_insert_walks_affected / (walk_update_time_on_insert.get_total() * 1.0) << endl;
		cout << "(2) average latency: " << fixed << setprecision(8) << average_latency / (n_batches * 1.0) << endl;
//		cout << "FindPrev vertex in node2vec: " << FindPreviousVertexNode2Vec.get_total() << endl;

//		cout << "(11) SIMPLE SEARCH throughput: " << fixed << setprecision(8) << total_insert_walks_affected / ((walk_update_time_on_insert.get_total() - FindPreviousVertexNode2Vec.get_total() + FindNextNodeSimpleSearch.get_total()) * 1.0) << endl;
//		cout << "(22) average latency: " << fixed << setprecision(8) << average_latency / (n_batches * 1.0) << endl;
//		cout << "FindPrev SIMPLE in node2vec: " << FindNextNodeSimpleSearch.get_total() << endl;

//		cout << "(111) RANGE SEARCH throughput: " << fixed << setprecision(8) << total_insert_walks_affected / ((walk_update_time_on_insert.get_total() - FindPreviousVertexNode2Vec.get_total() + FindNextNodeRangeSearch.get_total()) * 1.0) << endl;
//		cout << "(222) average latency: " << fixed << setprecision(8) << average_latency / (n_batches * 1.0) << endl;
//		cout << "FindPrev vertex in node2vec: " << FindNextNodeRangeSearch.get_total() << endl;

//		std::cout << "Average walk update latency = { ";
//		for (int i = 0; i < n_batches; i++) {
//			std::cout << latency[i] << " ";
//		}
//		std::cout << "}" << std::endl;

		// Last Merge
		LastMerge.start();
//	malin.merge_walk_trees_all_vertices_parallel(n_batches);
		malin.last_merge_all_vertices_parallel_with_minmax(n_batches);
		LastMerge.stop();
		cout << "Last merge (with MIN-MAX Ranges) time: " << LastMerge.get_total() << endl << endl;
	}



	// Measure time to read-rewalk all walks
/*	ReadWalks.start();
	for (auto i = 0; i < n * config::walks_per_vertex; i++)
//		cout << malin.walk_simple_find(i) << endl;
		malin.walk_silent(i);
	ReadWalks.stop();
	cout << "Read all walks time: " << ReadWalks.get_total() << endl;*/

	// Measure the memory after all batches
	malin.memory_footprint();

//	auto flat_graph = malin.flatten_vertex_tree();
//	for (auto i = 0; i < malin.number_of_vertices(); i++)
//	{
//		cout << "vertex " << i << endl;
////		flat_graph[i].compressed_edges.iter_elms(i, [&](auto edge){
////			cout << edge << " ";
////		});
////		cout << endl;
//
//		cout << "size of walk-tree vector " << flat_graph[i].compressed_walks.size() << endl;
//		int inc = 0;
//		for (auto wt = flat_graph[i].compressed_walks.begin(); wt != flat_graph[i].compressed_walks.end(); wt++) // print the walk-trees in chronological order
//		{
////			inc++;
//			cout << "walk-tree " << inc << endl;
//			inc++;
//			wt->iter_elms(i, [&](auto enc_triplet){
//			  auto pair = pairings::Szudzik<types::Vertex>::unpair(enc_triplet);
//
//			  auto walk_id  = pair.first / config::walk_length;                  // todo: needs floor?
//			  auto position = pair.first - (walk_id * config::walk_length); // todo: position here starts from 0. verify this one!
//			  auto next_vertex   = pair.second;
////				cout << enc_triplet << " ";
////			  cout << "{" << walk_id << ", " << position << ", " << next_vertex << "}" << " " << endl;
//			});
//			cout << endl;
//			cout << "size of walk-tree " << wt->size() << endl;
//		}
//	}
//
//// ----------------------------------------------
//	cout << "(NEW) WALKS" << endl;
//	for (auto i = 0; i < n * config::walks_per_vertex; i++)
//		cout << malin.walk_simple_find(i) << endl;
//// ----------------------------------------------

	std::cout << std::endl;
}

int main(int argc, char** argv)
{
	std::cout << "Running experiment with " << num_workers() << " threads" << std::endl;
	commandLine command_line(argc, argv, "");

	throughput(command_line);
}

