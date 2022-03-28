#include <gtest/gtest.h>

#include <wharf.h>
#include <math.h>

class SamplerTest : public testing::Test
{
    public:
        void SetUp()    final;
        void TearDown() final;

    protected:
        long  total_vertices{};
        long  total_edges{};
        uintE* offsets{};
        uintV* edges{};
};

void SamplerTest::SetUp()
{
    std::cout << "---------------------------------------------------------------------------------------" << std::endl;
    std::cout << "Running using " << num_workers() << " threads" << std::endl;

    total_vertices = 6;
    total_edges = 18;

    offsets = new uintE[total_vertices];
    offsets[0] = 0; offsets[3] = 10;
    offsets[1] = 2; offsets[4] = 13;
    offsets[2] = 5; offsets[5] = 15;

    edges = new uintV[total_edges];
    edges[0] = 1; edges[3] = 2; edges[6] = 1; edges[9] = 5;  edges[12] = 5; edges[15] = 2;
    edges[1] = 2; edges[4] = 3; edges[7] = 3; edges[10] = 1; edges[13] = 2; edges[16] = 3;
    edges[2] = 0; edges[5] = 0; edges[8] = 4; edges[11] = 2; edges[14] = 5; edges[17] = 4;
}

void SamplerTest::TearDown()
{
    std::cout << "---------------------------------------------------------------------------------------" << std::endl;
}

TEST_F(SamplerTest, DeepWalk)
{
    dygrl::Wharf malin = dygrl::Wharf(total_vertices, total_edges, offsets, edges);
    dygrl::FlatGraph graph = malin.flatten_graph();
    dygrl::RandomWalkModel* model = new dygrl::DeepWalk(&graph);
    std::map<types::Vertex, std::pair<float, int>> hash_table;

    types::State state = model->initial_state(config::random.irand(6));
    dygrl::MetropolisHastingsSampler sampler = dygrl::MetropolisHastingsSampler(state, model);

    std::cout << "Walker init state: (" << state.first << "," << state.second << ")" << std::endl;
    std::cout << "Previously sampled vertex: " << sampler.previously_sampled_vertex() << std::endl;

    auto neighbors = graph.neighbors(state.first);
    float total = 0;

    for (auto i = 0; i < std::get<1>(neighbors); i++)
    {
        auto neighbor = std::get<0>(neighbors)[i];
        auto weight = model->weight(state, neighbor);

        total += weight;
        hash_table[neighbor].first = weight;
        hash_table[neighbor].second = 0;
    }

    if (std::get<2>(neighbors)) pbbs::free_array(std::get<0>(neighbors));

    float iters = 10000;
    for (int i = 0; i < iters; i++)
    {
        auto neighbor = sampler.sample(state, model).first;
        hash_table[neighbor].second += 1;
    }

    for(const auto& item : hash_table)
    {
        std::cout << "Vertex: " << item.first
                  << " Ideal sample: " << (item.second.first / total) * iters
                  << " Sample: " << item.second.second << std::endl;
    }
}

TEST_F(SamplerTest, Node2Vec)
{
    dygrl::Wharf malin = dygrl::Wharf(total_vertices, total_edges, offsets, edges);
    dygrl::FlatGraph graph = malin.flatten_graph();
    dygrl::RandomWalkModel* model = new dygrl::Node2Vec(&graph, 0.7, 0.2);
    std::map<types::Vertex, std::pair<float, int>> hash_table;

    types::State state = model->initial_state(config::random.irand(6));
    dygrl::MetropolisHastingsSampler sampler = dygrl::MetropolisHastingsSampler(state, model);

    std::cout << "Walker init state: (" << state.first << "," << state.second << ")" << std::endl;
    std::cout << "Previously sampled vertex: " << sampler.previously_sampled_vertex() << std::endl;

    auto neighbors = graph.neighbors(state.first);
    float total = 0;

    for (auto i = 0; i < std::get<1>(neighbors); i++)
    {
        auto neighbor = std::get<0>(neighbors)[i];
        auto weight = model->weight(state, neighbor);

        total += weight;
        hash_table[neighbor].first = weight;
        hash_table[neighbor].second = 0;
    }

    if (std::get<2>(neighbors)) pbbs::free_array(std::get<0>(neighbors));

    int iters = 10000;
    for (int i = 0; i < iters; i++)
    {
        auto neighbor = sampler.sample(state, model).first;
        hash_table[neighbor].second += 1;
    }

    for(const auto& item : hash_table)
    {
        std::cout << "Vertex: " << item.first
                  << " Ideal sample: " << (item.second.first / total) * iters
                  << " Sample: " << item.second.second << std::endl;
    }
}
