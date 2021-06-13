#include <dock.h>

void throughput(commandLine& command_line)
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
    size_t n_trials         = command_line.getOptionLongValue("-trials", 5);

    config::walks_per_vertex = walks_per_vertex;
    config::walk_length      = length_of_walks;

    std::cout << "Walks per vertex: " << (int) config::walks_per_vertex << std::endl;
    std::cout << "Walk length: " <<  (int) config::walk_length << std::endl;

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

    auto batch_sizes = pbbs::sequence<size_t>(7);
    batch_sizes[0] = std::pow(10, 1);
    batch_sizes[1] = std::pow(10, 2);
    batch_sizes[2] = std::pow(10, 3);
    batch_sizes[3] = std::pow(10, 4);
    batch_sizes[4] = std::pow(10, 5);
    batch_sizes[5] = std::pow(10, 6);
    batch_sizes[6] = std::pow(10, 7);
    batch_sizes[7] = std::pow(10, 8);
    batch_sizes[8] = std::pow(10, 9);

    for (short int i = 0; i < batch_sizes.size(); i++)
    {
        timer insert_timer("InsertTimer");
        timer delete_timer("DeleteTimer");

        graph_update_time_on_insert.reset();
        walk_update_time_on_insert.reset();
        graph_update_time_on_delete.reset();
        walk_update_time_on_delete.reset();

        std::cout << "Batch size = " << 2*batch_sizes[i] << " | ";

        for (short int trial = 0; trial < n_trials; trial++)
        {
            size_t graph_size_pow2 = 1 << (pbbs::log2_up(n) - 1);
            auto edges = utility::generate_batch_of_edges(batch_sizes[i], n, false, false);

            std::cout << edges.second << " ";

            insert_timer.start();
            dock.insert_edges_batch(edges.second, edges.first, false, true, graph_size_pow2);
            insert_timer.stop();

            delete_timer.start();
            dock.delete_edges_batch(edges.second, edges.first, false, true, graph_size_pow2);
            delete_timer.stop();
        }

        std::cout << std::endl;

        std::cout << "Average insert time = " << insert_timer.get_total() / n_trials << std::endl;
        std::cout << "Average graph update insert time = " << graph_update_time_on_insert.get_total() / n_trials << std::endl;
        std::cout << "Average walk update insert time = " << walk_update_time_on_insert.get_total() / n_trials << std::endl;

        std::cout << "Average delete time = " << delete_timer.get_total() / n_trials << std::endl;
        std::cout << "Average graph update delete time = " << graph_update_time_on_delete.get_total() / n_trials << std::endl;
        std::cout << "Average walk update delete time = " << walk_update_time_on_delete.get_total() / n_trials << std::endl;

        std::cout << std::endl;
    }
}

int main(int argc, char** argv)
{
    std::cout << "Running experiment with " << num_workers() << " threads" << std::endl;
    commandLine command_line(argc, argv, "");

    throughput(command_line);
}

