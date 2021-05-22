#include <dock.h>

extern "C"
{
    extern void train_main(int argc, char **argv);
}

void learn_embeddings(commandLine& command_line)
{
    string fname             = string(command_line.getOptionValue("-f", default_file_name));
    bool mmap                = command_line.getOption("-m");
    bool is_symmetric        = command_line.getOption("-s");
    bool compressed          = command_line.getOption("-c");
    size_t n_parts           = command_line.getOptionLongValue("-nparts", 1);
    size_t walks_per_vertex  = command_line.getOptionLongValue("-w", config::walks_per_vertex);
    size_t length_of_walks   = command_line.getOptionLongValue("-l", config::walk_length);

    config::walks_per_vertex = walks_per_vertex;
    config::walk_length = length_of_walks;

    size_t n;
    size_t m;
    uintE* offsets;
    uintV* edges;

    std::tie(n, m, offsets, edges) = read_unweighted_graph(fname.c_str(), is_symmetric, mmap);

    dygrl::Dock dock = dygrl::Dock(n, m, offsets, edges);
    dock.create_random_walks();


    ofstream myfile;
    myfile.open("txt/all");

    for(int i = 0; i < config::walks_per_vertex * dock.number_of_vertices(); i++)
    {
        myfile << dock.rewalk(i) << "\n";
//        std::cout << dock.rewalk(i) << std::endl;
    }

    myfile.close();

    train_main(command_line.argc, command_line.argv);
}

int main(int argc, char** argv)
{
    std::cout << "Running experiment with: " << num_workers() << " threads." << std::endl;
    commandLine command_line(argc, argv, "??");
    learn_embeddings(command_line);
}

