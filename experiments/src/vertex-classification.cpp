#include <wharfmh.h>

using Edge = std::tuple<types::Vertex, types::Vertex>;

void create_edge_stream(commandLine& command_line, std::vector<std::pair<Edge*, uintV>>& stream)
{
//    string fname = string(command_line.getOptionValue("-f", "stream"));
    string fname = "data/stream";
    size_t edge_parition_size = command_line.getOptionLongValue("-eps", 5000);

    types::Vertex firstV, secondV;
    std::vector<std::vector<Edge>> temp;

    std::ifstream file(fname); int iter = 0;

    while (file >> firstV >> secondV)
    {
        auto index = iter % edge_parition_size;

        if (index == 0)
        {
            temp.emplace_back(std::vector<Edge>());
        }

        temp.back().push_back(std::make_tuple(firstV, secondV));
        temp.back().push_back(std::make_tuple(secondV, firstV));
        iter += 1;
    }

    for(auto& partition : temp)
    {
        stream.emplace_back(new Edge[partition.size()], partition.size());

        for(int j = 0; j < partition.size(); j++)
        {
            stream.back().first[j] = partition[j];
        }
    }
}

void vertex_classification_incremental(commandLine& command_line, const std::vector<std::pair<Edge*, uintV>>& stream)
{
    size_t walks_per_vertex  = command_line.getOptionLongValue("-w", config::walks_per_vertex);
    size_t length_of_walks   = command_line.getOptionLongValue("-l", config::walk_length);
    string model             = string(command_line.getOptionValue("-model", "deepwalk"));
    double paramP            = command_line.getOptionDoubleValue("-paramP", config::paramP);
    double paramQ            = command_line.getOptionDoubleValue("-paramQ", config::paramQ);
    string init_strategy     = string(command_line.getOptionValue("-init", "random"));
    size_t vector_dimension  = command_line.getOptionLongValue("-d", 128);
    size_t learning_strategy = command_line.getOptionLongValue("-le", 2);

    config::walks_per_vertex = walks_per_vertex;
    config::walk_length      = length_of_walks;

    std::cout << "Incremental Learning" << std::endl;
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

    string fname             = string(command_line.getOptionValue("-f", default_file_name));
    bool mmap                = command_line.getOption("-m");
    bool is_symmetric        = command_line.getOption("-s");

    size_t n; size_t m;
    uintE* offsets; uintV* edges;
    stringstream ss; ss << fname << ".adj";

    std::tie(n, m, offsets, edges) = read_unweighted_graph(ss.str().c_str(), is_symmetric, mmap);
//    pbbs::free_array(offsets);
//    pbbs::free_array(edges);

    timer incremental_timer("IncrementalTimer", false);

    /******************************************************************************************************************/

    // 1. load all vertices in isolation (every vertex has a degree of 0)
    incremental_timer.start();
    dygrl::WharfMH WharfMH = dygrl::WharfMH(n, m, offsets, edges);
    WharfMH.generate_initial_random_walks();
    incremental_timer.stop();

    // 2. train initial embeddings
    std::cout << "Learning initial embeddings" << std::endl;
    stringstream command; std::ofstream file("walks.txt");

    incremental_timer.start();
    for (int i = 0; i < n * config::walks_per_vertex; i++)
    {
        command << WharfMH.walk(i) << std::endl;
    }

    file << command.str(); file.close(); command.str(std::string());
    command << "yskip --thread-num="
            << num_workers()
            << " -d " << vector_dimension
            << " -l " << learning_strategy
            << " walks.txt model;";

    system(command.str().c_str());
    incremental_timer.stop();

    command.str(std::string());
    command << "perl to_word2vec.pl < model > model.w2v; python3 vertex-classification.py";
    system(command.str().c_str());

    // 3. train embeddings incrementally
    command.str(std::string());
    for (auto& edge_batch : stream)
    {
        incremental_timer.start();
        file.open("walks.txt");
        auto walks = WharfMH.insert_edges_batch(edge_batch.second, edge_batch.first, false, true);

        for (auto& walk_id : walks)
        {
            command << WharfMH.walk(walk_id) << std::endl;
        }

        file << command.str(); file.close(); command.str(std::string());

        command << "yskip --thread-num="
                << num_workers()
                << " -d " << vector_dimension
                << " -l " << learning_strategy
                << " walks.txt model;";
        system(command.str().c_str());
        incremental_timer.stop();

        command.str(std::string());
        command << "perl to_word2vec.pl < model > model.w2v; python3 vertex-classification.py";
        system(command.str().c_str());
        command.str(std::string());
    }

    incremental_timer.reportTotal("Total");
    WharfMH.destroy();
}

void vertex_classification_static(commandLine& command_line, const std::vector<std::pair<Edge*, uintV>>& stream)
{
    size_t walks_per_vertex  = command_line.getOptionLongValue("-w", config::walks_per_vertex);
    size_t length_of_walks   = command_line.getOptionLongValue("-l", config::walk_length);
    string model             = string(command_line.getOptionValue("-model", "deepwalk"));
    double paramP            = command_line.getOptionDoubleValue("-paramP", config::paramP);
    double paramQ            = command_line.getOptionDoubleValue("-paramQ", config::paramQ);
    string init_strategy     = string(command_line.getOptionValue("-init", "random"));
    size_t vector_dimension  = command_line.getOptionLongValue("-d", 128);
    size_t learning_strategy = 0;

    config::walks_per_vertex = walks_per_vertex;
    config::walk_length      = length_of_walks;

    std::cout << "Static Learning" << std::endl;
    std::cout << "Walks per vertex: " << (int) config::walks_per_vertex << std::endl;
    std::cout << "Walk length: " << (int) config::walk_length << std::endl;
    std::cout << "Vector dimension: " << vector_dimension << std::endl;
    std::cout << "Learning strategy: ";

    if(learning_strategy == 0)
    {
        std::cout << "BATCH" << std::endl;
    }
    else if(learning_strategy == 1)
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

    string fname             = string(command_line.getOptionValue("-f", default_file_name));
    bool mmap                = command_line.getOption("-m");
    bool is_symmetric        = command_line.getOption("-s");

    size_t n; size_t m;
    uintE* offsets; uintV* edges;
    stringstream ss; ss << fname << ".adj";

    std::tie(n, m, offsets, edges) = read_unweighted_graph(ss.str().c_str(), is_symmetric, mmap);
//    pbbs::free_array(offsets);
//    pbbs::free_array(edges);

    timer static_timer("StaticTimer", false);

    /******************************************************************************************************************/

    // 1. load all vertices in isolation (every vertex has a degree of 0)
    static_timer.start();
    dygrl::WharfMH WharfMH = dygrl::WharfMH(n, m, offsets, edges);
    WharfMH.generate_initial_random_walks();
    static_timer.stop();

    // 2. train initial embeddings
    std::cout << "Learning initial embeddings" << std::endl;
    stringstream command; std::ofstream file("walks.txt");

    static_timer.start();
    for (int i = 0; i < n * config::walks_per_vertex; i++)
    {
        command << WharfMH.walk(i) << std::endl;
    }

    file << command.str(); file.close(); command.str(std::string());
    command << "yskip --thread-num="
    << num_workers()
    << " -d " << vector_dimension
    << " -l " << learning_strategy
    << " walks.txt model;";

    system(command.str().c_str());
    static_timer.stop();

    command.str(std::string());
    command << "perl to_word2vec.pl < model > model.w2v; python3 vertex-classification.py";
    system(command.str().c_str());

    // 3. train embeddings incrementally
    command.str(std::string());
    for (auto& edge_batch : stream)
    {
        WharfMH.destroy_index();

        static_timer.start();
        file.open("walks.txt");
        WharfMH.insert_edges_batch(edge_batch.second, edge_batch.first, false, true, std::numeric_limits<size_t>::max(), false);
        WharfMH.generate_initial_random_walks();

        for (int i = 0; i < n * config::walks_per_vertex; i++)
        {
            command << WharfMH.walk(i) << std::endl;
        }

        file << command.str(); file.close(); command.str(std::string());

        command << "yskip --thread-num="
                << num_workers()
                << " --initial-model=model"
                << " -d " << vector_dimension
                << " -l " << learning_strategy
                << " walks.txt model;";
                system(command.str().c_str());

        static_timer.stop();

        command.str(std::string());
        command << "perl to_word2vec.pl < model > model.w2v; python3 vertex-classification.py";
        system(command.str().c_str());
        command.str(std::string());
    }

    static_timer.reportTotal("Total");
}

int main(int argc, char** argv)
{
    std::cout << " - Running vertex classification with "
              << num_workers() << " threads." << std::endl << std::endl;

    std::vector<std::pair<Edge*, uintV>> stream;
    commandLine command_line(argc, argv, "");

    create_edge_stream(command_line, stream);
    vertex_classification_incremental(command_line, stream);
    vertex_classification_static(command_line, stream);
}



