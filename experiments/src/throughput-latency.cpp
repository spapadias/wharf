#include <wharf.h>

void throughput(commandLine& command_line)
{
    string fname             = string(command_line.getOptionValue("-f", default_file_name));
    bool mmap                = command_line.getOption("-m");
    bool is_symmetric        = command_line.getOption("-s");
    bool compressed          = command_line.getOption("-c");
    size_t n_parts           = command_line.getOptionLongValue("-nparts", 1);

    size_t w = command_line.getOptionLongValue("-w", config::walks_per_vertex);
    size_t l  = command_line.getOptionLongValue("-l", config::walk_length);
    string model            = string(command_line.getOptionValue("-model", "deepwalk"));
    double p           = command_line.getOptionDoubleValue("-paramP", config::paramP);
    double q           = command_line.getOptionDoubleValue("-paramQ", config::paramQ);
    string init    = string(command_line.getOptionValue("-init", "weight"));
    size_t n_trials         = command_line.getOptionLongValue("-trials",  5);

    string determinism      = string(command_line.getOptionValue("-d", "false"));
    string range_search     = string(command_line.getOptionValue("-rs", "true"));
	size_t num_of_batches   = command_line.getOptionIntValue("-numbatch", 10);
	size_t half_of_bsize    = command_line.getOptionIntValue("-sizebatch", 5000);
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

	config::walks_per_vertex = w;
    config::walk_length      = l;

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
        config::paramP = p;
        config::paramQ = q;
        std::cout << "Walking model: NODE2VEC | Params (p,q) = " << "(" << config::paramP << "," << config::paramQ << ")" << std::endl;
    }
    else
    {
        std::cerr << "Unrecognized walking model! Abort" << std::endl;
        std::exit(1);
    }

    if (init == "burnin")
    {
        config::sampler_init_strategy = types::SamplerInitStartegy::BURNIN;
        std::cout << "Sampler strategy: BURNIN" << std::endl;
    }
    else if (init == "weight")
    {
        config::sampler_init_strategy = types::SamplerInitStartegy::WEIGHT;
        std::cout << "Sampler strategy: WEIGHT" << std::endl;
    }
    else if (init == "random")
    {
        config::sampler_init_strategy = types::SamplerInitStartegy::RANDOM;
        std::cout << "Sampler strategy: RANDOM" << std::endl;
    }
    else
    {
        std::cerr << "Unrecognized sampler init strategy" << std::endl;
        std::exit(1);
    }

    if (range_search == "true")
        config::range_search_mode = true;
    else
        config::range_search_mode = false;

    if (determinism == "true")
        config::deterministic_mode = true;
    else
        config::deterministic_mode = false;

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
	batch_sizes[0] = half_of_bsize;

	for (short int i = 0; i < batch_sizes.size(); i++)
	{
		timer insert_timer("InsertTimer");
		timer delete_timer("DeleteTimer");

		graph_update_time_on_insert.reset();
		walk_update_time_on_insert.reset();
		Walking_new_sampling_time.reset();
		Walking_insert_new_samples.reset();
		MAV_time.reset();
		MAV_min = 20000;
		MAV_max = 0;
		merge_calc_triplets_to_delete.reset();
		merge_create_delete_walks.reset();
		Merge_time.reset();
		Merge_min = 20000;
		Merge_max = 0;
		LastMerge.reset();

		std::cout << "Batch size = " << 2 * batch_sizes[i] << " | ";

		double last_insert_time = 0.0;
		double last_insert_total = 0.0;
		double last_MAV_time    = 0.0;
		double last_MAV_total   = 0.0;
		double last_Merge_time  = 0.0;
		double last_Merge_total = 0.0;
		double last_Walk_sampling_time = 0.0;
		double last_Walk_sampling_total = 0.0;
		double last_Walk_new_insert_time = 0.0;

		auto latency_insert = pbbs::sequence<double>(n_batches);
		auto latency = pbbs::sequence<double>(n_batches);

		double total_insert_walks_affected = 0;
		double total_delete_walks_affected = 0;

		int batch_seed[n_batches];
		for (auto i = 0; i < n_batches; i++)
			batch_seed[i] = i;

		for (short int b = 0; b < n_batches; b++)
		{
			cout << "batch-" << b << " and batch_seed-" << batch_seed[b] << endl;

			size_t graph_size_pow2 = 1 << (pbbs::log2_up(n) - 1);
			auto edges = utility::generate_batch_of_edges(batch_sizes[i], n, batch_seed[b], false, false);

			std::cout << edges.second << " ";

			insert_timer.start();
			auto x = malin.insert_edges_batch(edges.second, edges.first, b+1, false, true, graph_size_pow2);
			insert_timer.stop();

			total_insert_walks_affected += x.size();

			last_insert_time = walk_update_time_on_insert.get_total() - last_insert_total;
			last_insert_total =  walk_update_time_on_insert.get_total();
			latency_insert[b] = (double) last_insert_time / x.size();

			latency[b] = latency_insert[b];

			// Delete the batch of edges
			malin.delete_edges_batch(edges.second, edges.first, b+1, false, true, graph_size_pow2, false);

			// Update the MAV min and max
			last_MAV_time = MAV_time.get_total() - last_MAV_total;
			last_MAV_total = MAV_time.get_total();
			MAV_min = std::min(MAV_min, last_MAV_time);
			MAV_max = std::max(MAV_max, last_MAV_time);
			// Update the Merge min and max
			if ((b+1) % config::merge_frequency == 0) {
				last_Merge_time = Merge_time.get_total() - last_Merge_total;
				last_Merge_total = Merge_time.get_total();
				Merge_min = std::min(Merge_min, last_Merge_time);
				Merge_max = std::max(Merge_max, last_Merge_time);
			}
			// Update the Walk Update (without Merge) min and max
			last_Walk_sampling_time  = Walking_new_sampling_time.get_total() - last_Walk_sampling_total;
			last_Walk_sampling_total = Walking_new_sampling_time.get_total();
			WalkSampling_min = std::min(WalkSampling_min, last_Walk_sampling_time);
			WalkSampling_max = std::max(WalkSampling_max, last_Walk_sampling_time);
			last_Walk_new_insert_time = Walking_insert_new_samples.get_total() - last_Walk_new_insert_time;
			WalkInsert_min = std::min(WalkInsert_min, last_Walk_new_insert_time);
			WalkInsert_max = std::max(WalkInsert_max, last_Walk_new_insert_time);

			pbbs::free_array(edges.first);

			malin.memory_footprint();
		}
		cout << fixed;
		std::cout << std::endl;

		// Last Merge
		LastMerge.start();
		malin.last_merge_all_vertices_parallel_with_minmax(n_batches);
		LastMerge.stop();
	}

	malin.memory_footprint();

    std::cout << std::endl;
}

int main(int argc, char** argv)
{
    std::cout << "Running experiment with " << num_workers() << " threads" << std::endl;
    commandLine command_line(argc, argv, "");

    throughput(command_line);
}

