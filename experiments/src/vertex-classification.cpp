#include <dock.h>

void vertex_classification(commandLine& command_line)
{
    // 1. collect command line parameters
    string fname             = string(command_line.getOptionValue("-f", "email-graph"));
    bool mmap                = command_line.getOption("-m");
    bool is_symmetric        = command_line.getOption("-s");

    size_t walks_per_vertex  = command_line.getOptionLongValue("-w", config::walks_per_vertex);
    size_t length_of_walks   = command_line.getOptionLongValue("-l", config::walk_length);
    string model             = string(command_line.getOptionValue("-model", "deepwalk"));
    double paramP            = command_line.getOptionDoubleValue("-paramP", config::paramP);
    double paramQ            = command_line.getOptionDoubleValue("-paramQ", config::paramQ);
    string init_strategy     = string(command_line.getOptionValue("-init", "weight"));
    size_t vector_dimension  = command_line.getOptionLongValue("-d", 100);
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

    // 3. load graph and generate walks
    size_t n; size_t m;
    uintE* offsets; uintV* edges;
    std::tie(n, m, offsets, edges) = read_unweighted_graph(fname.c_str(), is_symmetric, mmap);

    dygrl::Dock dock = dygrl::Dock(n, m, offsets, edges);
    dock.create_random_walks();

    // 4. learn emebeddings
    std::cout << "Learning embeddings..." << std::endl;
    stringstream command;
//    command << "printf '";

    std::ofstream file("file.txt");

    for(int i = 0; i < n * walks_per_vertex; i++)
    {
        file << dock.rewalk(i) << std::endl;
//        command << dock.rewalk(i) << std::endl;
    }

    file.close();

//    command << "' | yskip --thread-num="
//            << num_workers()
//            << " -d " << vector_dimension
//            << " -l " << learning_strategy
//            << " - model;"
//            << "perl to_word2vec.pl < model > model.w2v";

    command << "yskip file.txt model; perl to_word2vec.pl < model > model.w2v";
    auto res = system(command.str().c_str());

//    // 5.
//    for(int i = 0; i < 1; i++)
//    {
//        // geneate edges
//        auto edges = utility::generate_batch_of_edges(10, dock.number_of_vertices(), false, false);
//        auto map = dock.insert_edges_batch(edges.second, edges.first, true, false);
//        command << "printf '";
//
//        for(auto& entry : map.lock_table())
//        {
//            command << dock.rewalk(entry.first) << std::endl;
//        }
//
//        command << "' | yskip --thread-num="
//                << num_workers()
//                << " --initial-model=model"
//                << " -d " << vector_dimension
//                << " -l " << learning_strategy
//                << " - model;"
//                << "perl to_word2vec.pl < model > model.w2v";
//
//        system(command.str().c_str());
//        command.str("");
//
//        map = dock.delete_edges_batch(edges.second, edges.first, true, false);
//        command << "printf '";
//
//        for(auto& entry : map.lock_table())
//        {
//            command << dock.rewalk(entry.first) << std::endl;
//        }
//
//        command << "' | yskip --thread-num="
//                << num_workers()
//                << " --initial-model=model"
//                << " -d " << vector_dimension
//                << " -l " << learning_strategy
//                << " - model;"
//                << "perl to_word2vec.pl < model > model.w2v";
//
//        system(command.str().c_str());
//        command.str("");
//    }
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



