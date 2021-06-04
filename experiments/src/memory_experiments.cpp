#include <dock.h>

void memory_footprint(commandLine& command_line)
{
    string fname             = string(command_line.getOptionValue("-f", default_file_name));
    bool mmap                = command_line.getOption("-m");
    bool is_symmetric        = command_line.getOption("-s");
    bool compressed          = command_line.getOption("-c");
    size_t n_parts           = command_line.getOptionLongValue("-nparts", 1);

    size_t walks_per_vertex = command_line.getOptionLongValue("-w", config::walks_per_vertex);
    size_t length_of_walks  = command_line.getOptionLongValue("-l", config::walk_length);
    string model            = string(command_line.getOptionValue("-model", "deepwalk"));
    double paramP           = command_line.getOptionDoubleValue("-paramP", config::paramP);
    double paramQ           = command_line.getOptionDoubleValue("-paramQ", config::paramQ);
    string init_strategy    = string(command_line.getOptionValue("-init", "weight"));


    config::walks_per_vertex = walks_per_vertex;
    config::walk_length      = length_of_walks;

    if (model == "deepwalk")
    {
        config::random_walk_model = types::RandomWalkModelType::DEEPWALK;
    }
    else if (model == "node2vec")
    {
        config::random_walk_model = types::RandomWalkModelType::NODE2VEC;
        config::paramP = paramP;
        config::paramQ = paramQ;
    }
    else
    {
        std::cout << "Unrecognized walking model" << std::endl;
        std::exit(1);
    }

    if (init_strategy == "burnin")
    {
        config::sampler_init_strategy = types::SamplerInitStartegy::BURNIN;
    }
    else if (init_strategy == "weight")
    {
        config::sampler_init_strategy = types::SamplerInitStartegy::WEIGHT;
    }
    else if (init_strategy == "random")
    {
        config::sampler_init_strategy = types::SamplerInitStartegy::RANDOM;
    }
    else
    {
        std::cout << "Unrecognized sampler init strategy" << std::endl;
        std::exit(1);
    }

    size_t n;
    size_t m;
    uintE* offsets;
    uintV* edges;
    std::cout << "Reading Unweighted Graph" << std::endl;
    std::tie(n, m, offsets, edges) = read_unweighted_graph(fname.c_str(), is_symmetric, mmap);
    std::cout << "Read Unweighted Graph" << std::endl;

    dygrl::Dock dock = dygrl::Dock(n, m, offsets, edges);
    dock.create_random_walks();
    dock.memory_footprint();

    ofstream walks;
    walks.open ("data/walks.txt");

    for(int i = 0; i < config::walks_per_vertex * dock.number_of_vertices(); i++)
    {
        walks << dock.rewalk(i) << std::endl;
    }

    walks.close();

    // transform an input graph file into an adjacency graph format
    std::string command = "./inverted_index_pam -f data/walks.txt ";
    system(command.c_str());
}

int main(int argc, char** argv)
{
    std::cout << "Running experiment with: " << num_workers() << " threads." << std::endl;
    commandLine command_line(argc, argv, "");
    memory_footprint(command_line);
}

