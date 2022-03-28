#include <wharf.h>

using Edge = std::tuple<uintV, uintV>;

void preprocess(commandLine& command_line, std::vector<Edge*>& stream)
{
    string fname      = string(command_line.getOptionValue("-f", "cora-graph"));
    int parition_size = (int) command_line.getOptionLongValue("-p", 1000);
    std::vector<std::vector<Edge>> temp_stream;

    uintV first, second;
    std::ifstream file(fname);

    int i = 0;
    while (file >> first >> second)
    {
        auto index = i % parition_size;

        if (index == 0)
        {
            temp_stream.emplace_back(std::vector<Edge>());
        }

        temp_stream.back().push_back(std::make_tuple(first, second));
        i += 1;
    }

    for(int i = 0; i < temp_stream.size(); i++)
    {
        stream.push_back(new Edge[temp_stream[i].size()]);

        for(int j = 0; j < temp_stream[i].size(); j++)
        {
            stream.back()[j] = temp_stream[i][j];
        }
    }

    std::cout << stream.size() << std::endl;
}

void vertex_classification(commandLine& command_line, const std::vector<Edge*>& stream)
{
    // 1. collect command line parameters
    size_t walks_per_vertex  = command_line.getOptionLongValue("-w", config::walks_per_vertex);
    size_t length_of_walks   = command_line.getOptionLongValue("-l", config::walk_length);
    string model             = string(command_line.getOptionValue("-model", "deepwalk"));
    double paramP            = command_line.getOptionDoubleValue("-paramP", config::paramP);
    double paramQ            = command_line.getOptionDoubleValue("-paramQ", config::paramQ);
    string init_strategy     = string(command_line.getOptionValue("-init", "weight"));
    size_t vector_dimension  = command_line.getOptionLongValue("-d", 128);
    size_t learning_strategy = command_line.getOptionLongValue("-le", 2);

    // 2. set up command line parameters
    config::walks_per_vertex = walks_per_vertex;
    config::walk_length      = length_of_walks;

    std::cout << "Walks per vertex: " << (int) config::walks_per_vertex << std::endl;
    std::cout << "Walk length: " << (int) config::walk_length << std::endl;
    std::cout << "Vector dimension: " << vector_dimension << std::endl;
    std::cout << "Learning strategy: ";

    if(learning_strategy == 1)
    {
        std::cout << "ONLINE" << std::endl;
    }
    else if(learning_strategy == 2)
    {
        std::cout << "MINI-BATCH" << std::endl;
    }
    else
    {
        std::cerr << "Unrecognized learning strategy! Abort" << std::endl;
        std::exit(1);
    }

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
}

int main(int argc, char** argv)
{
    std::cout << std::endl
              << "Running experiment with: "
              << num_workers()
              << " threads."
              << std::endl
              << std::endl;

    std::vector<Edge*> stream;

    commandLine command_line(argc, argv, "");
    preprocess(command_line, stream);
    vertex_classification(command_line, stream);
}



