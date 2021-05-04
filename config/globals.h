#ifndef DYNAMIC_GRAPH_REPRESENTATION_LEARNING_WITH_METROPOLIS_HASTINGS_GLOBALS_H
#define DYNAMIC_GRAPH_REPRESENTATION_LEARNING_WITH_METROPOLIS_HASTINGS_GLOBALS_H

namespace config
{
    // determines the number of walks per node to be generated by the walkers
    const unsigned char walks_per_vertex = 200;

    // determines the length of one random walk
    const unsigned char walk_length      = 100;
}

#endif // DYNAMIC_GRAPH_REPRESENTATION_LEARNING_WITH_METROPOLIS_HASTINGS_GLOBALS_H
