#include <dock.h>

void vertex_classification(commandLine& command_line)
{
    string fname             = string(command_line.getOptionValue("-f", "email-graph"));
    bool mmap                = command_line.getOption("-m");
    bool is_symmetric        = command_line.getOption("-s");

    size_t walks_per_vertex = command_line.getOptionLongValue("-w", config::walks_per_vertex);
    size_t length_of_walks  = command_line.getOptionLongValue("-l", config::walk_length);
    string model            = string(command_line.getOptionValue("-model", "deepwalk"));
    double paramP           = command_line.getOptionDoubleValue("-paramP", config::paramP);
    double paramQ           = command_line.getOptionDoubleValue("-paramQ", config::paramQ);
    string init_strategy    = string(command_line.getOptionValue("-init", "weight"));

    config::walks_per_vertex = walks_per_vertex;
    config::walk_length      = length_of_walks;

    std::cout << "Walks per vertex: " << (int) config::walks_per_vertex << std::endl;
    std::cout << "Walk length: " << (int) config::walk_length << std::endl;

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

    size_t n;
    size_t m;
    uintE* offsets;
    uintV* edges;
    std::tie(n, m, offsets, edges) = read_unweighted_graph(fname.c_str(), is_symmetric, mmap);

    dygrl::Dock dock = dygrl::Dock(n, m, offsets, edges);
    dock.create_random_walks();

    std::cout << "Learning embeddings..." << std::endl;

    stringstream stream;
    stream << "printf '";

    for(int i = 0; i < n * walks_per_vertex; i++)
    {
        stream << dock.rewalk(i) << std::endl;
    }

    stream << "'";
    stream << " | yskip - model";
    system(stream.str().c_str());
}

int main(int argc, char** argv)
{
    std::cout << std::endl
              << "Running experiment with: "
              << num_workers()
              << " threads."
              << std::endl
              << std::endl;

    commandLine command_line(argc, argv, "");
    vertex_classification(command_line);
}



