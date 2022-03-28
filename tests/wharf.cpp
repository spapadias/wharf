#include <gtest/gtest.h>

#include <wharf.h>

class WharfTest : public testing::Test
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
        std::string default_file_path = "data/aspen-paper-graph";
};

void WharfTest::SetUp()
{
    std::cout << "-----------------------------------------------------------------------------------------------------" << std::endl;
    std::cout << "Wharf running with " << num_workers() << " threads" << std::endl;

    // transform an input graph file into an adjacency graph format
    std::string command = "./SNAPtoAdj -s -f " + this->default_file_path + " data/adjacency-graph-format.txt";
    int result = system(command.c_str());

    if (result != 0)
    {
        std::cerr << "WharfTest::SetUp::Input file could not be transformed!" << std::endl;
        exit(1);
    }

    std::tie(total_vertices, total_edges, offsets, edges) = read_unweighted_graph("data/adjacency-graph-format.txt", is_symmetric, mmap);
    std::cout << std::endl;
}

void WharfTest::TearDown()
{
    // remove adjaceny graph format representation
    int graph = system("rm -rf data/adjacency-graph-format.txt");

    if (graph != 0)
    {
        std::cerr << "WharfTest::TearDown::Could not remove static graph input file" << std::endl;
    }

    std::cout << "-----------------------------------------------------------------------------------------------------" << std::endl;
}

TEST_F(WharfTest, WharfConstructor)
{
    dygrl::Wharf wharf = dygrl::Wharf(total_vertices, total_edges, offsets, edges, false);

    // assert the number of vertices and edges in a graph
    ASSERT_EQ(wharf.number_of_vertices(), total_vertices);
    ASSERT_EQ(wharf.number_of_edges(), total_edges);

    // construct a flat snapshot of a graph
    auto flat_snapshot = wharf.flatten_vertex_tree();

    // assert
    parallel_for(0, total_vertices, [&] (long i)
    {
        size_t off = offsets[i];
        size_t degree = ((i == (total_vertices - 1)) ? total_edges : offsets[i+1]) - off;
        auto S = pbbs::delayed_seq<uintV>(degree, [&] (size_t j) { return edges[off + j]; });

        // assert expected degrees
        ASSERT_EQ(flat_snapshot[i].compressed_edges.degree(), degree);

        // assert that compressed_walks_tree is empty
        ASSERT_EQ(flat_snapshot[i].compressed_walks.size(), 0);

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

TEST_F(WharfTest, DockDestructor)
{
    dygrl::Wharf wharf = dygrl::Wharf(total_vertices, total_edges, offsets, edges);

    wharf.print_memory_pool_stats();
    wharf.destroy();
    wharf.print_memory_pool_stats();

    // assert vertices and edges
    ASSERT_EQ(wharf.number_of_vertices(), 0);
    ASSERT_EQ(wharf.number_of_edges(), 0);

    // construct a flat snapshot of a graph
    auto flat_snapshot = wharf.flatten_vertex_tree();

    // assert that flat snapshot does not exits
    ASSERT_EQ(flat_snapshot.size(), 0);
}

TEST_F(WharfTest, WharfDestroyIndex)
{
    dygrl::Wharf wharf = dygrl::Wharf(total_vertices, total_edges, offsets, edges);
    wharf.generate_initial_random_walks();

    wharf.print_memory_pool_stats();
    wharf.destroy_index();
    wharf.print_memory_pool_stats();

    // assert vertices and edges
    ASSERT_EQ(wharf.number_of_vertices(), total_vertices);
    ASSERT_EQ(wharf.number_of_edges(), total_edges);

    // construct a flat snapshot of a graph
    auto flat_snapshot = wharf.flatten_vertex_tree();

    parallel_for(0, total_vertices, [&] (long i)
    {
//        ASSERT_EQ(flat_snapshot[i].compressed_walks.size(), 0); // todo: removed it
        ASSERT_EQ(flat_snapshot[i].sampler_manager->size(), 0);
    });
}

TEST_F(WharfTest, InsertBatchOfEdges)
{
    // create wharf instance (vertices & edges)
    dygrl::Wharf wharf = dygrl::Wharf(total_vertices, total_edges, offsets, edges);
    auto start_edges = wharf.number_of_edges();

    // geneate edges
    auto edges = utility::generate_batch_of_edges(1000000, wharf.number_of_vertices(), false, false);

    // insert batch of edges
    wharf.insert_edges_batch(edges.second, edges.first, true, false, std::numeric_limits<size_t>::max(), false);

    std::cout << "Edges before batch insert: " << start_edges << std::endl;
    std::cout << "Edges after batch insert: "  << wharf.number_of_edges() << std::endl;

    // assert edge insertion
    ASSERT_GE(wharf.number_of_edges(), start_edges);
}

// Pump up the test to debug diff updates approach
TEST_F(WharfTest, InsertBatchOfEdgesPlayground)
{
	// create wharf instance (vertices & edges)
	dygrl::Wharf wharf = dygrl::Wharf(total_vertices, total_edges, offsets, edges);
	auto start_edges = wharf.number_of_edges();

	wharf.generate_initial_random_walks();

	// geneate edges
	auto edges = utility::generate_batch_of_edges(10, wharf.number_of_vertices(), 10, false, false);

	// insert batch of edges
	wharf.insert_edges_batch(edges.second, edges.first, true, false/*, std::numeric_limits<size_t>::max(), false*/);

	auto flat_graph = wharf.flatten_vertex_tree();
	for (auto i = 0; i < wharf.number_of_vertices(); i++)
	{
		cout << "vertex " << i << endl;

		int inc = 0;
		for (auto wt = flat_graph[i].compressed_walks.begin(); wt != flat_graph[i].compressed_walks.end(); wt++)
		{
			inc++;
			cout << "walk-tree " << inc << endl;
			wt->iter_elms(i, [&](auto enc_triplet){
				auto pair = pairings::Szudzik<types::Vertex>::unpair(enc_triplet);

				auto walk_id  = pair.first / config::walk_length;
				auto position = pair.first - (walk_id * config::walk_length);
				auto next_vertex   = pair.second;
				cout << "{" << walk_id << ", " << position << ", " << next_vertex << "}" << " " << endl;
			});
			cout << endl;
		}
	}

	std::cout << "Edges before batch insert: " << start_edges << std::endl;
	std::cout << "Edges after batch insert: "  << wharf.number_of_edges() << std::endl;

	// assert edge insertion
	ASSERT_GE(wharf.number_of_edges(), start_edges);
}

TEST_F(WharfTest, DeleteBatchOfEdges)
{
    // create wharf instance (vertices & edges)
    dygrl::Wharf wharf = dygrl::Wharf(total_vertices, total_edges, offsets, edges);
    auto start_edges = wharf.number_of_edges();

    // geneate edges
    auto edges = utility::generate_batch_of_edges(1000000, wharf.number_of_vertices(), false, false);

    // insert batch of edges
    wharf.delete_edges_batch(edges.second, edges.first, true, false, std::numeric_limits<size_t>::max(), false);

    std::cout << "Edges before batch delete: " << start_edges << std::endl;
    std::cout << "Edges after batch delete: "  << wharf.number_of_edges() << std::endl;

    // assert edge deletion
    ASSERT_LE(wharf.number_of_edges(), start_edges);
}

TEST_F(WharfTest, UpdateRandomWalksOnInsertEdges)
{
    // create graph and walks
    dygrl::Wharf wharf = dygrl::Wharf(total_vertices, total_edges, offsets, edges);
    wharf.generate_initial_random_walks();

    // geneate edges
    auto edges = utility::generate_batch_of_edges(1000, wharf.number_of_vertices(), false, false);

    // insert batch of edges
    wharf.insert_edges_batch(edges.second, edges.first, true, false);
}

TEST_F(WharfTest, UpdateRandomWalksOnDeleteEdges)
{
    // create graph and walks
    dygrl::Wharf wharf = dygrl::Wharf(total_vertices, total_edges, offsets, edges);
    wharf.generate_initial_random_walks();

    // geneate edges
    auto edges = utility::generate_batch_of_edges(1000, wharf.number_of_vertices(), false, false);

    // insert batch of edges
    wharf.delete_edges_batch(edges.second, edges.first, true, false);
}

TEST_F(WharfTest, UpdateRandomWalks)
{
    // create graph and walks
    dygrl::Wharf wharf = dygrl::Wharf(total_vertices, total_edges, offsets, edges);
    wharf.generate_initial_random_walks();

    for(int i = 0; i < 10; i++)
    {
        auto edges = utility::generate_batch_of_edges(1000, wharf.number_of_vertices(), false, false);

        wharf.insert_edges_batch(edges.second, edges.first, true, false);
        wharf.delete_edges_batch(edges.second, edges.first, true, false);
    }
}

// -------------------------------//
// ---- tests for range search ---//
// -------------------------------//

TEST_F(WharfTest, GenerateAndPrintInitialRW)
{
    // create graph and walks
    dygrl::Wharf wharf = dygrl::Wharf(total_vertices, total_edges, offsets, edges);
    wharf.generate_initial_random_walks();

    auto walk_printing_timer = timer("TotalTimeToRewalk", false);
    double time_simple, time_range;

    walk_printing_timer.start();
    time_simple = walk_printing_timer.get_total();

    cout << "-- now print out the walks with the find_in_range operation" << endl;

    walk_printing_timer.reset(); walk_printing_timer.start();
    time_range = walk_printing_timer.get_total();

    cout << "Time to print all walk corpus with simple find: " << time_simple << endl
         << " and to print all walk corpus with range  find: " << time_range  << endl;
}

// -----------------------------------
// -----------------------------------

TEST_F(WharfTest, UpdateRandomWalksWithRangeSearch)
{
    // create graph and walks
    dygrl::Wharf wharf = dygrl::Wharf(total_vertices, total_edges, offsets, edges);
    wharf.generate_initial_random_walks();

    // geneate edges
    auto edges = utility::generate_batch_of_edges(1000, wharf.number_of_vertices(), false, false);
    wharf.insert_edges_batch(edges.second, edges.first, true, false);
}

TEST_F(WharfTest, BatchGenerator)
{
	cout << "first batch" << endl;
	auto edges = utility::generate_batch_of_edges(10, 6, false, false);
	for (auto i = 0; i < edges.second; i++)
		cout << "edge-" << i + 1 << " is [" << get<0>(edges.first[i]) << ", " << get<1>(edges.first[i]) << "]" << endl;

	cout << "second batch" << endl;
	edges = utility::generate_batch_of_edges(10, 6, false, false);
	for (auto i = 0; i < edges.second; i++)
		cout << "edge-" << i + 1 << " is [" << get<0>(edges.first[i]) << ", " << get<1>(edges.first[i]) << "]" << endl;

	cout << "third batch" << endl;
	edges = utility::generate_batch_of_edges(10, 6, false, false);
	for (auto i = 0; i < edges.second; i++)
		cout << "edge-" << i + 1 << " is [" << get<0>(edges.first[i]) << ", " << get<1>(edges.first[i]) << "]" << endl;
}