#include <malin.h>

using Edge = std::tuple<types::Vertex, types::Vertex>;

void create_edge_stream(commandLine& command_line, std::vector<Edge*>& stream)
{
    string fname = string(command_line.getOptionValue("-f", "cora-graph"));
    size_t edge_parition_size = command_line.getOptionLongValue("-eps", 1000);

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
        iter += 1;
    }

    for(auto& partition : temp)
    {
        stream.push_back(new Edge[partition.size()]);

        for(int j = 0; j < partition.size(); j++)
        {
            stream.back()[j] = partition[j];
        }
    }
}

void vertex_classification(commandLine& command_line, const std::vector<Edge*>& stream)
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
    std::tie(n, m, offsets, edges) = read_unweighted_graph(fname.c_str(), is_symmetric, mmap);

    pbbs::free_array(offsets);
    pbbs::free_array(edges);

    
    dygrl::Malin malin = dygrl::Malin(n, m);

//    malin.generate_initial_random_walks();
//
//    // 4. learn emebeddings
//    std::cout << "Learning embeddings..." << std::endl;
//    stringstream command;
////    command << "printf '";
//
//    std::ofstream file("file.txt");
//
//    for(int i = 0; i < n * walks_per_vertex; i++)
//    {
//        file << dock.rewalk(i) << std::endl;
////        command << dock.rewalk(i) << std::endl;
//    }
//
//    file.close();
//
////    command << "' | yskip --thread-num="
////            << num_workers()
////            << " -d " << vector_dimension
////            << " -l " << learning_strategy
////            << " - model;"
////            << "perl to_word2vec.pl < model > model.w2v";
//
//    command << "yskip file.txt model; perl to_word2vec.pl < model > model.w2v";
//    auto res = system(command.str().c_str());
//
////    // 5.
////    for(int i = 0; i < 1; i++)
////    {
////        // geneate edges
////        auto edges = utility::generate_batch_of_edges(10, dock.number_of_vertices(), false, false);
////        auto map = dock.insert_edges_batch(edges.second, edges.first, true, false);
////        command << "printf '";
////
////        for(auto& entry : map.lock_table())
////        {
////            command << dock.rewalk(entry.first) << std::endl;
////        }
////
////        command << "' | yskip --thread-num="
////                << num_workers()
////                << " --initial-model=model"
////                << " -d " << vector_dimension
////                << " -l " << learning_strategy
////                << " - model;"
////                << "perl to_word2vec.pl < model > model.w2v";
////
////        system(command.str().c_str());
////        command.str("");
////
////        map = dock.delete_edges_batch(edges.second, edges.first, true, false);
////        command << "printf '";
////
////        for(auto& entry : map.lock_table())
////        {
////            command << dock.rewalk(entry.first) << std::endl;
////        }
////
////        command << "' | yskip --thread-num="
////                << num_workers()
////                << " --initial-model=model"
////                << " -d " << vector_dimension
////                << " -l " << learning_strategy
////                << " - model;"
////                << "perl to_word2vec.pl < model > model.w2v";
////
////        system(command.str().c_str());
////        command.str("");
////    }
}

int main(int argc, char** argv)
{
    std::cout << " - Running vertex classification with "
              << num_workers() << " threads." << std::endl << std::endl;

    std::vector<Edge*> stream;
    commandLine command_line(argc, argv, "");

    create_edge_stream(command_line, stream);
    vertex_classification(command_line, stream);
}



