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
        bool mmap = false;          // TODO @Djordjije: do we need this?
        bool is_symmetric = true;   // TODO @Djordjije: do we need this?
        std::string default_file_path = "data/email-graph";
};

void DockTest::SetUp()
{
    std::cout << "-----------------------------------------------------------------------------------------------------"
              << std::endl;
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

    std::cout << "-----------------------------------------------------------------------------------------------------"
              << std::endl;
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

        auto edges = flat_snapshot[i].compressed_edges.get_edges(i);

        // assert expected neighbours
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

TEST_F(DockTest, DockCreateRandomWalks)
{
    dygrl::Dock dock = dygrl::Dock(total_vertices, total_edges, offsets, edges);
    dock.create_random_walks();

    for(auto i = 0; i < dock.number_of_vertices()*config::walks_per_vertex; i++)
    {
        std::cout << dock.rewalk(i) << std::endl;
    }
}

