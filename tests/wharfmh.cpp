#include <gtest/gtest.h>

#include <wharfmh.h>

class WharfMHTest : public testing::Test
{
    public:
        void SetUp()    override;
        void TearDown() override;

    protected:
        long total_vertices;
        long total_edges;
        uintV* edges;
        uintE* offsets;
        bool mmap = false;
        bool is_symmetric = true;
        std::string default_file_path = "data/cora-graph";
};

void WharfMHTest::SetUp()
{
    std::cout << "-----------------------------------------------------------------------------------------------------" << std::endl;
    std::cout << "WharfMH running with " << num_workers() << " threads" << std::endl;

    // transform an input graph file into an adjacency graph format
    std::string command = "./SNAPtoAdj -s -f " + this->default_file_path + " data/adjacency-graph-format.txt";
    int result = system(command.c_str());

    if (result != 0)
    {
        std::cerr << "WharfMHTest::SetUp::Input file could not be transformed!" << std::endl;
        exit(1);
    }

    std::tie(total_vertices, total_edges, offsets, edges) = read_unweighted_graph("data/adjacency-graph-format.txt", is_symmetric, mmap);
    std::cout << std::endl;
}

void WharfMHTest::TearDown()
{
    // remove adjaceny graph format representation
    int graph = system("rm -rf data/adjacency-graph-format.txt");

    if (graph != 0)
    {
        std::cerr << "WharfMHTest::TearDown::Could not remove static graph input file" << std::endl;
    }

    std::cout << "-----------------------------------------------------------------------------------------------------" << std::endl;
}

TEST_F(WharfMHTest, WharfMHConstructor)
{
    dygrl::Wharf WharfMH = dygrl::Wharf(total_vertices, total_edges, offsets, edges, false);

    // assert the number of vertices and edges in a graph
    ASSERT_EQ(WharfMH.number_of_vertices(), total_vertices);
    ASSERT_EQ(WharfMH.number_of_edges(), total_edges);

    // construct a flat snapshot of a graph
    auto flat_snapshot = WharfMH.flatten_vertex_tree();

    // assert
    parallel_for(0, total_vertices, [&] (long i)
    {
        size_t off = offsets[i];
        size_t degree = ((i == (total_vertices - 1)) ? total_edges : offsets[i+1]) - off;
        auto S = pbbs::delayed_seq<uintV>(degree, [&] (size_t j) { return edges[off + j]; });

        // assert expected degrees
        ASSERT_EQ(flat_snapshot[i].compressed_edges.degree(), degree);

        // assert that compressed_walks_tree is empty
        ASSERT_EQ(flat_snapshot[i].inverted_index.root, nullptr);
        ASSERT_EQ(flat_snapshot[i].inverted_index.size(), 0);

        // assert empty samplers
        ASSERT_EQ(flat_snapshot[i].sampler_manager->size(), 0);

        // assert expected neighbours
        auto edges = flat_snapshot[i].compressed_edges.get_edges(i);

        for(auto j = 0; j < degree; j++)
        {
            bool flag = false;

            for (auto k = 0; k <  S.size(); k++)
            {
                if (S[k] == edges[j]) flag = true;
            }

            ASSERT_EQ(flag, true);
        }
    });
}

TEST_F(WharfMHTest, WharfMHDestructor)
{
    dygrl::Wharf WharfMH = dygrl::Wharf(total_vertices, total_edges, offsets, edges);

    WharfMH.print_memory_pool_stats();
    WharfMH.destroy();
    WharfMH.print_memory_pool_stats();

    // assert vertices and edges
    ASSERT_EQ(WharfMH.number_of_vertices(), 0);
    ASSERT_EQ(WharfMH.number_of_edges(), 0);

    // construct a flat snapshot of a graph
    auto flat_snapshot = WharfMH.flatten_vertex_tree();

    // assert that flat snapshot does not exits
    ASSERT_EQ(flat_snapshot.size(), 0);
}

TEST_F(WharfMHTest, WharfMHDestroyIndex)
{
    dygrl::Wharf WharfMH = dygrl::Wharf(total_vertices, total_edges, offsets, edges);
    WharfMH.generate_initial_random_walks();

    WharfMH.print_memory_pool_stats();
    WharfMH.destroy_index();
    WharfMH.print_memory_pool_stats();

    // assert vertices and edges
    ASSERT_EQ(WharfMH.number_of_vertices(), total_vertices);
    ASSERT_EQ(WharfMH.number_of_edges(), total_edges);

    // construct a flat snapshot of a graph
    auto flat_snapshot = WharfMH.flatten_vertex_tree();

    parallel_for(0, total_vertices, [&] (long i)
    {
        ASSERT_EQ(flat_snapshot[i].inverted_index.size(), 0);
    });
}

TEST_F(WharfMHTest, InsertBatchOfEdges)
{
    // create wharf instance (vertices & edges)
    dygrl::Wharf WharfMH = dygrl::Wharf(total_vertices, total_edges, offsets, edges);
    auto start_edges = WharfMH.number_of_edges();

//    // print random walks
//    for(int i = 0; i < config::walks_per_vertex * WharfMH.number_of_vertices(); i++)
//        std::cout << WharfMH.walk(i) << std::endl;

    // geneate edges
    auto edges = utility::generate_batch_of_edges(100, WharfMH.number_of_vertices(), false, false);

    // insert batch of edges
    WharfMH.insert_edges_batch(edges.second, edges.first, true, false, std::numeric_limits<size_t>::max(), false);

//    // print random walks
//    for(int i = 0; i < config::walks_per_vertex * WharfMH.number_of_vertices(); i++)
//        std::cout << WharfMH.walk(i) << std::endl;

    std::cout << "Edges before batch insert: " << start_edges << std::endl;
    std::cout << "Edges after batch insert: "  << WharfMH.number_of_edges() << std::endl;

    // assert edge insertion
    ASSERT_GE(WharfMH.number_of_edges(), start_edges);
}

TEST_F(WharfMHTest, DeleteBatchOfEdges)
{
    // create wharf instance (vertices & edges)
    dygrl::Wharf WharfMH = dygrl::Wharf(total_vertices, total_edges, offsets, edges);
    auto start_edges = WharfMH.number_of_edges();

    // geneate edges
    auto edges = utility::generate_batch_of_edges(1000000, WharfMH.number_of_vertices(), false, false);

    // insert batch of edges
    WharfMH.delete_edges_batch(edges.second, edges.first, true, false, std::numeric_limits<size_t>::max(), false);

    std::cout << "Edges before batch delete: " << start_edges << std::endl;
    std::cout << "Edges after batch delete: " << WharfMH.number_of_edges() << std::endl;

    // assert edge deletion
    ASSERT_LE(WharfMH.number_of_edges(), start_edges);
}

TEST_F(WharfMHTest, UpdateRandomWalksOnInsertEdges)
{
    // create graph and walks
    dygrl::Wharf WharfMH = dygrl::Wharf(total_vertices, total_edges, offsets, edges);
    WharfMH.generate_initial_random_walks();

    // print random walks
    for(int i = 0; i < config::walks_per_vertex * WharfMH.number_of_vertices(); i++)
    {
        std::cout << WharfMH.walk(i) << std::endl;
    }

    // geneate edges
    auto edges = utility::generate_batch_of_edges(1000, WharfMH.number_of_vertices(), false, false);

    // insert batch of edges
    WharfMH.insert_edges_batch(edges.second, edges.first, true, false);

    // print updated random walks
    for(int i = 0; i < config::walks_per_vertex * WharfMH.number_of_vertices(); i++)
    {
        std::cout << WharfMH.walk(i) << std::endl;
    }
}

TEST_F(WharfMHTest, UpdateRandomWalksOnDeleteEdges)
{
    // create graph and walks
    dygrl::Wharf WharfMH = dygrl::Wharf(total_vertices, total_edges, offsets, edges);
    WharfMH.generate_initial_random_walks();

    // print random walks
    for(int i = 0; i < config::walks_per_vertex * WharfMH.number_of_vertices(); i++)
    {
        std::cout << WharfMH.walk(i) << std::endl;
    }

    // geneate edges
    auto edges = utility::generate_batch_of_edges(1000, WharfMH.number_of_vertices(), false, false);

    // insert batch of edges
    WharfMH.delete_edges_batch(edges.second, edges.first, true, false);

    // print updated random walks
    for(int i = 0; i < config::walks_per_vertex * WharfMH.number_of_vertices(); i++)
    {
        std::cout << WharfMH.walk(i) << std::endl;
    }
}

TEST_F(WharfMHTest, UpdateRandomWalks)
{
    // create graph and walks
    dygrl::Wharf WharfMH = dygrl::Wharf(total_vertices, total_edges, offsets, edges);
    WharfMH.generate_initial_random_walks();

    // print random walks
    for(int i = 0; i < config::walks_per_vertex * WharfMH.number_of_vertices(); i++)
    {
        std::cout << WharfMH.walk(i) << std::endl;
    }

    for(int i = 0; i < 10; i++)
    {
        // geneate edges
        auto edges = utility::generate_batch_of_edges(1000, WharfMH.number_of_vertices(), false, false);

        WharfMH.insert_edges_batch(edges.second, edges.first, true, false);
        WharfMH.delete_edges_batch(edges.second, edges.first, true, false);
    }

    // print random walks
    for(int i = 0; i < config::walks_per_vertex * WharfMH.number_of_vertices(); i++)
    {
        std::cout << WharfMH.walk(i) << std::endl;
    }
}

TEST_F(WharfMHTest, WharfMHMemoryFootprint)
{
    dygrl::Wharf WharfMH = dygrl::Wharf(total_vertices, total_edges, offsets, edges);
    WharfMH.generate_initial_random_walks();
    WharfMH.memory_footprint();
}

TEST_F(WharfMHTest, WharfMHThroughputLatency)
{
    dygrl::Wharf WharfMH = dygrl::Wharf(total_vertices, total_edges, offsets, edges);
    WharfMH.generate_initial_random_walks();
    int n_trials = 3;

    auto batch_sizes = pbbs::sequence<size_t>(6);
    batch_sizes[0] = 5;
    batch_sizes[1] = 50;
    batch_sizes[2] = 500;
    batch_sizes[3] = 5000;
    batch_sizes[4] = 50000;
    batch_sizes[5] = 500000;

    for (short int i = 0; i < batch_sizes.size(); i++)
    {
        timer insert_timer("InsertTimer");
        timer delete_timer("DeleteTimer");

        graph_update_time_on_insert.reset();
        walk_update_time_on_insert.reset();
        graph_update_time_on_delete.reset();
        walk_update_time_on_delete.reset();

        std::cout << "Batch size = " << 2*batch_sizes[i] << " | ";

        double last_insert_time = 0;
        double last_delete_time = 0;

        auto latency_insert = pbbs::sequence<double>(n_trials);
        auto latency_delete = pbbs::sequence<double>(n_trials);
        auto latency        = pbbs::sequence<double>(n_trials);

        double total_insert_walks_affected = 0;
        double total_delete_walks_affected = 0;

        for (short int trial = 0; trial < n_trials; trial++)
        {
            size_t graph_size_pow2 = 1 << (pbbs::log2_up(total_vertices) - 1);
            auto edges = utility::generate_batch_of_edges(batch_sizes[i], total_vertices, false, false);

            std::cout << edges.second << " ";

            insert_timer.start();
            auto x = WharfMH.insert_edges_batch(edges.second, edges.first, false, true, graph_size_pow2);
            insert_timer.stop();

            total_insert_walks_affected += x.size();

            last_insert_time = walk_update_time_on_insert.get_total() - last_insert_time;
            latency_insert[trial] = last_insert_time / x.size();

            delete_timer.start();
            auto y = WharfMH.delete_edges_batch(edges.second, edges.first, false, true, graph_size_pow2);
            delete_timer.stop();

            total_delete_walks_affected += y.size();

            last_delete_time = walk_update_time_on_delete.get_total() - last_delete_time;
            latency_delete[trial] = last_delete_time / y.size();

            latency[trial] = (last_insert_time + last_delete_time) / (x.size() + y.size());

            // free edges
            pbbs::free_array(edges.first);
        }

        std::cout << std::endl;

        std::cout << "Average insert time = " << insert_timer.get_total() / n_trials << std::endl;
        std::cout << "Average graph update insert time = " << graph_update_time_on_insert.get_total() / n_trials << std::endl;
        std::cout << "Average walk update insert time = "
                  << walk_update_time_on_insert.get_total() / n_trials
                  << ", average walk affected = " << total_insert_walks_affected / n_trials << std::endl;

        std::cout << "Average delete time = " << delete_timer.get_total() / n_trials << std::endl;
        std::cout << "Average graph update delete time = " << graph_update_time_on_delete.get_total() / n_trials << std::endl;
        std::cout << "Average walk update delete time = "
                  << walk_update_time_on_delete.get_total() / n_trials
                  << ", average walk affected = " << total_delete_walks_affected / n_trials << std::endl;

        std::cout << "Average walk insert latency = { ";
        for(int i = 0; i < n_trials; i++)
        {
            std::cout << latency_insert[i] << " ";
        }
        std::cout << "}" << std::endl;

        std::cout << "Average walk delete latency = { ";
        for(int i = 0; i < n_trials; i++)
        {
            std::cout << latency_delete[i] << " ";
        }
        std::cout << "}" << std::endl;

        std::cout << "Average walk update latency = { ";
        for(int i = 0; i < n_trials; i++)
        {
            std::cout << latency[i] << " ";
        }
        std::cout << "}" << std::endl;
    }

    WharfMH.destroy_index();
    timer generate_initial_walks("Generate Initial Random Walks", false);

    for (int i = 0; i < n_trials; i++)
    {
        generate_initial_walks.start();
        WharfMH.generate_initial_random_walks();
        generate_initial_walks.stop();

        WharfMH.destroy_index();
    }


    std::cout << std::endl
              << "Average time to generate random walks from scratch = "
              << generate_initial_walks.get_total() / n_trials << std::endl;

    std::cout << std::endl;
}

TEST_F(WharfMHTest, WharfMHThroughputLatency2)
{
	dygrl::Wharf WharfMH = dygrl::Wharf(total_vertices, total_edges, offsets, edges);
	WharfMH.generate_initial_random_walks();
	int n_trials = 1; //3;

	cout << "WALK CORPUS" << endl;
	for (auto i = 0; i < total_vertices * config::walks_per_vertex; i++)
		cout << WharfMH.walk(i) << endl;

	auto batch_sizes = pbbs::sequence<size_t>(1);
	batch_sizes[0] = 10;
//	batch_sizes[1] = 50;
//	batch_sizes[2] = 500;
//	batch_sizes[3] = 5000;
//	batch_sizes[4] = 50000;
//	batch_sizes[5] = 500000;

	for (short int i = 0; i < batch_sizes.size(); i++)
	{
		timer insert_timer("InsertTimer");
		timer delete_timer("DeleteTimer");

		graph_update_time_on_insert.reset();
		walk_update_time_on_insert.reset();
		graph_update_time_on_delete.reset();
		walk_update_time_on_delete.reset();

		std::cout << "Batch size = " << 2*batch_sizes[i] << " | ";

		double last_insert_time = 0;
		double last_delete_time = 0;

		auto latency_insert = pbbs::sequence<double>(n_trials);
		auto latency_delete = pbbs::sequence<double>(n_trials);
		auto latency        = pbbs::sequence<double>(n_trials);

		double total_insert_walks_affected = 0;
		double total_delete_walks_affected = 0;

		for (short int trial = 0; trial < n_trials; trial++)
		{
			size_t graph_size_pow2 = 1 << (pbbs::log2_up(total_vertices) - 1);
			cout << "graph size pow: " << graph_size_pow2 << endl;
			cout << "batch size: " << batch_sizes[i] << endl;
			cout << "total V: " << total_vertices << endl;
			auto edges = utility::generate_batch_of_edges(batch_sizes[i], total_vertices, false, false);
			// ---
			for (auto i = 0; i < edges.second; i++)
				cout << "edge-" << i + 1 << " is [" << get<0>(edges.first[i]) << ", " << get<1>(edges.first[i]) << "]" << endl;

			std::cout << edges.second << " ";

			insert_timer.start();
			auto x = WharfMH.insert_edges_batch(edges.second, edges.first, false, true, graph_size_pow2);
			insert_timer.stop();

			total_insert_walks_affected += x.size();

			last_insert_time = walk_update_time_on_insert.get_total() - last_insert_time;
			latency_insert[trial] = last_insert_time / x.size();

//			delete_timer.start();
//			auto y = WharfMH.delete_edges_batch(edges.second, edges.first, false, true, graph_size_pow2);
//			delete_timer.stop();

//			total_delete_walks_affected += y.size();

//			last_delete_time = walk_update_time_on_delete.get_total() - last_delete_time;
//			latency_delete[trial] = last_delete_time / y.size();

//			latency[trial] = (last_insert_time + last_delete_time) / (x.size() + y.size());
			latency[trial] = latency_insert[trial];

			// free edges
			pbbs::free_array(edges.first);
		}

		std::cout << std::endl;

		std::cout << "Average insert time = " << insert_timer.get_total() / n_trials << std::endl;
		std::cout << "Average graph update insert time = " << graph_update_time_on_insert.get_total() / n_trials << std::endl;
		std::cout << "Average walk update insert time = "
		          << walk_update_time_on_insert.get_total() / n_trials
		          << ", average walk affected = " << total_insert_walks_affected / n_trials << std::endl;

//		std::cout << "Average delete time = " << delete_timer.get_total() / n_trials << std::endl;
//		std::cout << "Average graph update delete time = " << graph_update_time_on_delete.get_total() / n_trials << std::endl;
//		std::cout << "Average walk update delete time = "
//		          << walk_update_time_on_delete.get_total() / n_trials
//		          << ", average walk affected = " << total_delete_walks_affected / n_trials << std::endl;

		std::cout << "Average walk insert latency = { ";
		for(int i = 0; i < n_trials; i++)
		{
			std::cout << latency_insert[i] << " ";
		}
		std::cout << "}" << std::endl;

//		std::cout << "Average walk delete latency = { ";
//		for(int i = 0; i < n_trials; i++)
//		{
//			std::cout << latency_delete[i] << " ";
//		}
//		std::cout << "}" << std::endl;

		std::cout << "Average walk update latency = { ";
		for(int i = 0; i < n_trials; i++)
		{
			std::cout << latency[i] << " ";
		}
		std::cout << "}" << std::endl;
	}

	cout << "(NEW) WALK CORPUS" << endl;
	for (auto i = 0; i < total_vertices * config::walks_per_vertex; i++)
		cout << WharfMH.walk(i) << endl;

	WharfMH.destroy_index();
	timer generate_initial_walks("Generate Initial Random Walks", false);

	for (int i = 0; i < n_trials; i++)
	{
		generate_initial_walks.start();
		WharfMH.generate_initial_random_walks();
		generate_initial_walks.stop();

		WharfMH.destroy_index();
	}


	std::cout << std::endl
	          << "Average time to generate random walks from scratch = "
	          << generate_initial_walks.get_total() / n_trials << std::endl;

	std::cout << std::endl;
}

// ----------------------------
// --- INSERT ONLY WORKLOAD ---
// ----------------------------

TEST_F(WharfMHTest, WharfInsertOnlyWorkload) {
	dygrl::Wharf malin = dygrl::Wharf(total_vertices, total_edges, offsets, edges);
	malin.generate_initial_random_walks();
	int n_batches = 1; // todo: how many batches per batch size?

	auto batch_sizes = pbbs::sequence<size_t>(1);
	batch_sizes[0] = 5; //5;
//	batch_sizes[1] = 50;
//	batch_sizes[2] = 500;
//	batch_sizes[3] = 5000;
//	batch_sizes[4] = 50000;
//  batch_sizes[5] = 500000;

	for (short int i = 0; i < batch_sizes.size(); i++)
	{
		timer insert_timer("InsertTimer");
		timer delete_timer("DeleteTimer");

		graph_update_time_on_insert.reset();
		walk_update_time_on_insert.reset();
		graph_update_time_on_delete.reset();
		walk_update_time_on_delete.reset();
//		// --- profiling initialization
//		walk_insert_init.reset();
//		walk_insert_2jobs.reset();
//		walk_insert_2accs.reset();
//		ij.reset();
//		dj.reset();
//		walk_find_in_vertex_tree.reset();
//		walk_find_next_tree.reset();
//		szudzik_hash.reset();
//		fnir_tree_search.reset();
//		MAV_time.reset();
//		// ---

		std::cout << "Batch size = " << 2 * batch_sizes[i] << " | ";

		double last_insert_time = 0;

		auto latency_insert = pbbs::sequence<double>(n_batches);
		auto latency = pbbs::sequence<double>(n_batches);

		double total_insert_walks_affected = 0;
		double total_delete_walks_affected = 0;

		int batch_seed[n_batches];
		for (auto i = 0; i < n_batches; i++)
			batch_seed[i] = i; // say the seed equals to the #batch todo: produce a different batch each time

		for (short int b = 0; b < n_batches; b++)
		{
			cout << "batch-" << b << " and batch_seed-" << batch_seed[b] << endl;

			size_t graph_size_pow2 = 1 << (pbbs::log2_up(total_vertices) - 1);
			auto edges = utility::generate_batch_of_edges(batch_sizes[i], total_vertices, batch_seed[b], false, false);

			std::cout << edges.second << " ";
			for (auto i = 0; i < edges.second; i++)
				cout << "edge-" << i + 1 << " is [" << get<0>(edges.first[i]) << ", " << get<1>(edges.first[i]) << "]" << endl;

			insert_timer.start();
			auto x = malin.insert_edges_batch(edges.second, edges.first, b+1, false, true, graph_size_pow2); // pass the batch number as well
			insert_timer.stop();

			total_insert_walks_affected += x.size();

			last_insert_time = walk_update_time_on_insert.get_total() - last_insert_time;
			latency_insert[b] = (double) last_insert_time / x.size();

			latency[b] = latency_insert[b];

			// free edges
			pbbs::free_array(edges.first);
		}

		std::cout << std::endl;

		std::cout << "Average insert time = "
		          << insert_timer.get_total() / n_batches << std::endl;
		std::cout << "Average graph update insert time = "
		          << graph_update_time_on_insert.get_total() / n_batches
		          << std::endl;
		std::cout << "Average walk update insert time = "
		          << walk_update_time_on_insert.get_total() / n_batches
		          << ", average walk affected = "
		          << total_insert_walks_affected / n_batches << std::endl;

		std::cout << "Average delete time = "
		          << delete_timer.get_total() / n_batches << std::endl;
		std::cout << "Average graph update delete time = "
		          << graph_update_time_on_delete.get_total() / n_batches
		          << std::endl;
		std::cout << "Average walk update delete time = "
		          << walk_update_time_on_delete.get_total() / n_batches
		          << ", average walk affected = "
		          << total_delete_walks_affected / n_batches << std::endl;

		// MAV time
//		std::cout << "Average MAV (we are not deleting obsolete parts) = "
//		          << MAV_time.get_total() / n_batches
//		          << std::endl;
//
//		// --- profiling ---
//		std::cout << "{ total profiling for insert and delete" << std::endl;
//		std::cout << "Initialization: "
//		          << walk_insert_init.get_total() / n_batches << " ("
//		          << (walk_insert_init.get_total() * 100) /
//		             (walk_insert_init.get_total() +
//		              walk_insert_2jobs.get_total() +
//		              walk_insert_2accs.get_total()) << "%)" << std::endl;
//		std::cout << "Insert/Delete Jobs: "
//		          << walk_insert_2jobs.get_total() / n_batches << " ("
//		          << (walk_insert_2jobs.get_total() * 100) /
//		             (walk_insert_init.get_total() +
//		              walk_insert_2jobs.get_total() +
//		              walk_insert_2accs.get_total()) << "%)" << std::endl;
//		std::cout << "InsertJob: " << ij.get_total() / n_batches
//		          << " | DeleteJob: " << dj.get_total() / n_batches << std::endl;
//		std::cout << "FindInVertexTree in DeleteJob total: "
//		          << walk_find_in_vertex_tree.get_total() / n_batches
//		          << std::endl;
//		std::cout << "FindNext in DeleteJob total: "
//		          << walk_find_next_tree.get_total() / n_batches << std::endl;
//		std::cout << "FindNext (search of the tree): "
//		          << fnir_tree_search.get_total() / n_batches << std::endl;
//		std::cout << "Sudzik total: " << szudzik_hash.get_total() / n_batches
//		          << std::endl;
//
//		std::cout << "Accumulators: "
//		          << walk_insert_2accs.get_total() / n_batches << " ("
//		          << (walk_insert_2accs.get_total() * 100) /
//		             (walk_insert_init.get_total() +
//		              walk_insert_2jobs.get_total() +
//		              walk_insert_2accs.get_total()) << "%)" << std::endl;
//		std::cout << "}" << std::endl;
//		// --- profiling ---

		// latencies
		std::cout << "Average walk insert latency = { ";
		for (int i = 0; i < n_batches; i++) {
			std::cout << latency_insert[i] << " ";
		}
		std::cout << "}" << std::endl;

		std::cout << "Average walk update latency = { ";
		for (int i = 0; i < n_batches; i++) {
			std::cout << latency[i] << " ";
		}
		std::cout << "}" << std::endl;
	}

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
//			inc++;
//			cout << "walk-tree " << inc << endl;
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
//		}
//	}
}