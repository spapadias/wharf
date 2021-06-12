#include <gtest/gtest.h>

#include <dock.h>

class DockTest : public testing::Test
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
        std::string default_file_path = "data/email-graph";
};

void DockTest::SetUp()
{
    std::cout << "-----------------------------------------------------------------------------------------------------" << std::endl;
    std::cout << "Dock running with " << num_workers() << " threads" << std::endl;

    // transform an input graph file into an adjacency graph format
    std::string command = "./SNAPtoAdj -s -f " + this->default_file_path + " data/adjacency-graph-format.txt";
    int result = system(command.c_str());

    if (result != 0)
    {
        std::cerr << "DockTest::SetUp::Input file could not be transformed!" << std::endl;
        exit(1);
    }

    std::tie(total_vertices, total_edges, offsets, edges) = read_unweighted_graph("data/adjacency-graph-format.txt", is_symmetric, mmap);
    std::cout << std::endl;
}

void DockTest::TearDown()
{
    // remove adjaceny graph format representation
    int graph = system("rm -rf data/adjacency-graph-format.txt");

    if (graph != 0)
    {
        std::cerr << "DockTest::TearDown::Could not remove static graph input file" << std::endl;
    }

    std::cout << "-----------------------------------------------------------------------------------------------------" << std::endl;
}

TEST_F(DockTest, DockConstructor)
{
    dygrl::Dock dock = dygrl::Dock(total_vertices, total_edges, offsets, edges, false);

    // assert the number of vertices and edges in a graph
    ASSERT_EQ(dock.number_of_vertices(), total_vertices);
    ASSERT_EQ(dock.number_of_edges(), total_edges);

    // construct a flat snapshot of a graph
    auto flat_snapshot = dock.flatten_vertex_tree();

    // assert
    parallel_for(0, total_vertices, [&] (long i)
    {
        size_t off = offsets[i];
        size_t degree = ((i == (total_vertices - 1)) ? total_edges : offsets[i+1]) - off;
        auto S = pbbs::delayed_seq<uintV>(degree, [&] (size_t j) { return edges[off + j]; });

        // assert expected degrees
        ASSERT_EQ(flat_snapshot[i].compressed_edges.degree(), degree);

        // assert that compressed_walks_tree is empty
        ASSERT_EQ(flat_snapshot[i].compressed_walks.root, nullptr);
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

TEST_F(DockTest, DockDestructor)
{
    dygrl::Dock dock = dygrl::Dock(total_vertices, total_edges, offsets, edges);

    dock.print_memory_pool_stats();
    dock.destroy();
    dock.print_memory_pool_stats();

    // assert vertices and edges
    ASSERT_EQ(dock.number_of_vertices(), 0);
    ASSERT_EQ(dock.number_of_edges(), 0);

    // construct a flat snapshot of a graph
    auto flat_snapshot = dock.flatten_vertex_tree();

    // assert that flat snapshot does not exits
    ASSERT_EQ(flat_snapshot.size(), 0);
}

TEST_F(DockTest, InsertBatchOfEdges)
{
    // create wharf instance (vertices & edges)
    dygrl::Dock dock = dygrl::Dock(total_vertices, total_edges, offsets, edges);
    auto start_edges = dock.number_of_edges();

    // geneate edges
    auto edges = utility::generate_batch_of_edges(1000000, dock.number_of_vertices(), false, false);
    auto generated_edges = edges.first;
    auto edges_generated = edges.second;

    // insert batch of edges
    dock.insert_edges_batch(edges_generated, generated_edges, true, false, std::numeric_limits<size_t>::max(), false);

    std::cout << "Edges before batch insert: " << start_edges << std::endl;
    std::cout << "Edges after batch insert: "  << dock.number_of_edges() << std::endl;

    // assert edge insertion
    ASSERT_GE(dock.number_of_edges(), start_edges);
}

TEST_F(DockTest, DeleteBatchOfEdges)
{
    // create wharf instance (vertices & edges)
    dygrl::Dock dock = dygrl::Dock(total_vertices, total_edges, offsets, edges);
    auto start_edges = dock.number_of_edges();

    // geneate edges
    auto edges = utility::generate_batch_of_edges(1000000, dock.number_of_vertices(), false, false);
    auto generated_edges = edges.first;
    auto edges_generated = edges.second;

    // insert batch of edges
    dock.delete_edges_batch(edges_generated, generated_edges, true, false, std::numeric_limits<size_t>::max(), false);

    std::cout << "Edges before batch delete: " << start_edges << std::endl;
    std::cout << "Edges after batch delete: " << dock.number_of_edges() << std::endl;

    // assert edge deletion
    ASSERT_LE(dock.number_of_edges(), start_edges);
}

TEST_F(DockTest, UpdateRandomWalksOnInsertEdges)
{
    // create graph and walks
    dygrl::Dock dock = dygrl::Dock(total_vertices, total_edges, offsets, edges);
    dock.create_random_walks();

    // print random walks
    for(int i = 0; i < config::walks_per_vertex * dock.number_of_vertices(); i++)
    {
        std::cout << dock.rewalk(i) << std::endl;
    }

    // geneate edges
    auto edges = utility::generate_batch_of_edges(1000000, dock.number_of_vertices(), false, false);
    auto generated_edges = edges.first;
    auto edges_generated = edges.second;

    // insert batch of edges
    dock.insert_edges_batch(edges_generated, generated_edges, true, false);

    // print updated random walks
    for(int i = 0; i < config::walks_per_vertex * dock.number_of_vertices(); i++)
    {
        std::cout << dock.rewalk(i) << std::endl;
    }
}

TEST_F(DockTest, UpdateRandomWalksOnDeleteEdges)
{
    // create graph and walks
    dygrl::Dock dock = dygrl::Dock(total_vertices, total_edges, offsets, edges);
    dock.create_random_walks();

    // print random walks
    for(int i = 0; i < config::walks_per_vertex * dock.number_of_vertices(); i++)
    {
        std::cout << dock.rewalk(i) << std::endl;
    }

    // geneate edges
    auto edges = utility::generate_batch_of_edges(1000000, dock.number_of_vertices(), false, false);
    auto generated_edges = edges.first;
    auto edges_generated = edges.second;

    // insert batch of edges
    dock.delete_edges_batch(edges_generated, generated_edges, true, false);

    // print updated random walks
    for(int i = 0; i < config::walks_per_vertex * dock.number_of_vertices(); i++)
    {
        std::cout << dock.rewalk(i) << std::endl;
    }
}

TEST_F(DockTest, UpdateRandomWalks)
{
    // create graph and walks
    dygrl::Dock dock = dygrl::Dock(total_vertices, total_edges, offsets, edges);
    dock.create_random_walks();

    // print random walks
    for(int i = 0; i < config::walks_per_vertex * dock.number_of_vertices(); i++)
    {
        std::cout << dock.rewalk(i) << std::endl;
    }

    for(int i = 0; i < 10; i++)
    {
        // geneate edges
        auto edges = utility::generate_batch_of_edges(1000000, dock.number_of_vertices(), false, false);
        auto generated_edges = edges.first;
        auto edges_generated = edges.second;

        dock.insert_edges_batch(edges_generated, generated_edges, true, false);
        dock.delete_edges_batch(edges_generated, generated_edges, true, false);
    }

    // print random walks
    for(int i = 0; i < config::walks_per_vertex * dock.number_of_vertices(); i++)
    {
        std::cout << dock.rewalk(i) << std::endl;
    }
}