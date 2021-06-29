#include <gtest/gtest.h>

#include <malin.h>

class MalinTest : public testing::Test
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
        std::string default_file_path = "data/youtube-graph";
};

void MalinTest::SetUp()
{
    std::cout << "-----------------------------------------------------------------------------------------------------" << std::endl;
    std::cout << "Malin running with " << num_workers() << " threads" << std::endl;

    // transform an input graph file into an adjacency graph format
    std::string command = "./SNAPtoAdj -s -f " + this->default_file_path + " data/adjacency-graph-format.txt";
    int result = system(command.c_str());

    if (result != 0)
    {
        std::cerr << "MalinTest::SetUp::Input file could not be transformed!" << std::endl;
        exit(1);
    }

    std::tie(total_vertices, total_edges, offsets, edges) = read_unweighted_graph("data/adjacency-graph-format.txt", is_symmetric, mmap);
    std::cout << std::endl;
}

void MalinTest::TearDown()
{
    // remove adjaceny graph format representation
    int graph = system("rm -rf data/adjacency-graph-format.txt");

    if (graph != 0)
    {
        std::cerr << "MalinTest::TearDown::Could not remove static graph input file" << std::endl;
    }

    std::cout << "-----------------------------------------------------------------------------------------------------" << std::endl;
}

TEST_F(MalinTest, MalinConstructor)
{
    dygrl::Malin malin = dygrl::Malin(total_vertices, total_edges, offsets, edges, false);

    // assert the number of vertices and edges in a graph
    ASSERT_EQ(malin.number_of_vertices(), total_vertices);
    ASSERT_EQ(malin.number_of_edges(), total_edges);

    // construct a flat snapshot of a graph
    auto flat_snapshot = malin.flatten_vertex_tree();

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

TEST_F(MalinTest, MalinDestructor)
{
    dygrl::Malin malin = dygrl::Malin(total_vertices, total_edges, offsets, edges);

    malin.print_memory_pool_stats();
    malin.destroy();
    malin.print_memory_pool_stats();

    // assert vertices and edges
    ASSERT_EQ(malin.number_of_vertices(), 0);
    ASSERT_EQ(malin.number_of_edges(), 0);

    // construct a flat snapshot of a graph
    auto flat_snapshot = malin.flatten_vertex_tree();

    // assert that flat snapshot does not exits
    ASSERT_EQ(flat_snapshot.size(), 0);
}

TEST_F(MalinTest, InsertBatchOfEdges)
{
    // create wharf instance (vertices & edges)
    dygrl::Malin malin = dygrl::Malin(total_vertices, total_edges, offsets, edges);
    auto start_edges = malin.number_of_edges();

    // geneate edges
    auto edges = utility::generate_batch_of_edges(1000000, malin.number_of_vertices(), false, false);

    // insert batch of edges
    malin.insert_edges_batch(edges.second, edges.first, true, false, std::numeric_limits<size_t>::max(), false);

    std::cout << "Edges before batch insert: " << start_edges << std::endl;
    std::cout << "Edges after batch insert: "  << malin.number_of_edges() << std::endl;

    // assert edge insertion
    ASSERT_GE(malin.number_of_edges(), start_edges);
}

TEST_F(MalinTest, DeleteBatchOfEdges)
{
    // create wharf instance (vertices & edges)
    dygrl::Malin malin = dygrl::Malin(total_vertices, total_edges, offsets, edges);
    auto start_edges = malin.number_of_edges();

    // geneate edges
    auto edges = utility::generate_batch_of_edges(1000000, malin.number_of_vertices(), false, false);

    // insert batch of edges
    malin.delete_edges_batch(edges.second, edges.first, true, false, std::numeric_limits<size_t>::max(), false);

    std::cout << "Edges before batch delete: " << start_edges << std::endl;
    std::cout << "Edges after batch delete: " << malin.number_of_edges() << std::endl;

    // assert edge deletion
    ASSERT_LE(malin.number_of_edges(), start_edges);
}

TEST_F(MalinTest, UpdateRandomWalksOnInsertEdges)
{
    // create graph and walks
    dygrl::Malin malin = dygrl::Malin(total_vertices, total_edges, offsets, edges);
    malin.generate_initial_random_walks();

    // print random walks
    for(int i = 0; i < config::walks_per_vertex * malin.number_of_vertices(); i++)
    {
        std::cout << malin.walk(i) << std::endl;
    }

    // geneate edges
    auto edges = utility::generate_batch_of_edges(1000000, malin.number_of_vertices(), false, false);

    // insert batch of edges
    malin.insert_edges_batch(edges.second, edges.first, true, false);

    // print updated random walks
    for(int i = 0; i < config::walks_per_vertex * malin.number_of_vertices(); i++)
    {
        std::cout << malin.walk(i) << std::endl;
    }
}

TEST_F(MalinTest, UpdateRandomWalksOnDeleteEdges)
{
    // create graph and walks
    dygrl::Malin malin = dygrl::Malin(total_vertices, total_edges, offsets, edges);
    malin.generate_initial_random_walks();

    // print random walks
    for(int i = 0; i < config::walks_per_vertex * malin.number_of_vertices(); i++)
    {
        std::cout << malin.walk(i) << std::endl;
    }

    // geneate edges
    auto edges = utility::generate_batch_of_edges(1000000, malin.number_of_vertices(), false, false);

    // insert batch of edges
    malin.delete_edges_batch(edges.second, edges.first, true, false);

    // print updated random walks
    for(int i = 0; i < config::walks_per_vertex * malin.number_of_vertices(); i++)
    {
        std::cout << malin.walk(i) << std::endl;
    }
}

TEST_F(MalinTest, UpdateRandomWalks)
{
    // create graph and walks
    dygrl::Malin malin = dygrl::Malin(total_vertices, total_edges, offsets, edges);
    malin.generate_initial_random_walks();

    // print random walks
    for(int i = 0; i < config::walks_per_vertex * malin.number_of_vertices(); i++)
    {
        std::cout << malin.walk(i) << std::endl;
    }

    for(int i = 0; i < 10; i++)
    {
        // geneate edges
        auto edges = utility::generate_batch_of_edges(1000000, malin.number_of_vertices(), false, false);

        malin.insert_edges_batch(edges.second, edges.first, true, false);
        malin.delete_edges_batch(edges.second, edges.first, true, false);
    }

    // print random walks
    for(int i = 0; i < config::walks_per_vertex * malin.number_of_vertices(); i++)
    {
        std::cout << malin.walk(i) << std::endl;
    }
}

TEST_F(MalinTest, MalinMemoryFootprint)
{
    dygrl::Malin malin = dygrl::Malin(total_vertices, total_edges, offsets, edges);
    malin.generate_initial_random_walks();
    malin.memory_footprint();
}

TEST_F(MalinTest, MalinThroughputLatency)
{
    dygrl::Malin malin = dygrl::Malin(total_vertices, total_edges, offsets, edges);
    malin.generate_initial_random_walks();
    int n_trials = 5;

    auto batch_sizes = pbbs::sequence<size_t>(7);
    batch_sizes[0] = std::pow(10, 1);           // 10
    batch_sizes[1] = std::pow(10, 2);           // 100
    batch_sizes[2] = std::pow(10, 3);           // 1.000
    batch_sizes[3] = std::pow(10, 4);           // 10.000
    batch_sizes[4] = std::pow(10, 5);           // 100.000
    batch_sizes[5] = std::pow(10, 6);           // 1.000.000
    batch_sizes[6] = std::pow(10, 7);           // 10.000.000

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

        for (short int trial = 0; trial < n_trials; trial++)
        {
            size_t graph_size_pow2 = 1 << (pbbs::log2_up(total_vertices) - 1);
            auto edges = utility::generate_batch_of_edges(batch_sizes[i], total_vertices, false, false);

            std::cout << edges.second << " ";

            insert_timer.start();
            auto x = malin.insert_edges_batch(edges.second, edges.first, false, true, graph_size_pow2);
            insert_timer.stop();

            last_insert_time = walk_update_time_on_insert.get_total() - last_insert_time;
            latency_insert[trial] = last_insert_time / x.size();

            delete_timer.start();
            auto y = malin.delete_edges_batch(edges.second, edges.first, false, true, graph_size_pow2);
            delete_timer.stop();

            last_delete_time = walk_update_time_on_delete.get_total() - last_delete_time;
            latency_delete[trial] = last_delete_time / y.size();

            latency[trial] = (last_insert_time + last_delete_time) / (x.size() + y.size());

            // free edges
            pbbs::free_array(edges.first);
        }

        std::cout << std::endl;

        std::cout << "Average insert time = " << insert_timer.get_total() / n_trials << std::endl;
        std::cout << "Average graph update insert time = " << graph_update_time_on_insert.get_total() / n_trials << std::endl;
        std::cout << "Average walk update insert time = " << walk_update_time_on_insert.get_total() / n_trials << std::endl;

        std::cout << "Average delete time = " << delete_timer.get_total() / n_trials << std::endl;
        std::cout << "Average graph update delete time = " << graph_update_time_on_delete.get_total() / n_trials << std::endl;
        std::cout << "Average walk update delete time = " << walk_update_time_on_delete.get_total() / n_trials << std::endl;

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

        std::cout << std::endl;
    }
}