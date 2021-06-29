#ifndef DYNAMIC_GRAPH_REPRESENTATION_LEARNING_WITH_METROPOLIS_HASTINGS_GLOBALS_H
#define DYNAMIC_GRAPH_REPRESENTATION_LEARNING_WITH_METROPOLIS_HASTINGS_GLOBALS_H

namespace config
{
    // determines the number of walks per vertex to be generated by the walkers
    uint8_t walks_per_vertex   = 1;

    // determines the length of one random walk
    uint8_t walk_length        = 7;

    // determines the type of random walk model to use
    auto random_walk_model     = types::RandomWalkModelType::DEEPWALK;

    // determines parameter P for node2vec model
    float paramP               = 0.7;

    // determines parameter Q for node2vec model
    float paramQ               = 0.3;

    // determines the initialization strategy for metropolis hastings samplers
    auto sampler_init_strategy = types::SamplerInitStartegy::BURNIN;

    // random number generator
    auto random                = utility::Random(std::time(nullptr));
}

#endif // DYNAMIC_GRAPH_REPRESENTATION_LEARNING_WITH_METROPOLIS_HASTINGS_GLOBALS_H
