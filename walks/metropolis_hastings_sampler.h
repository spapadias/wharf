#ifndef DYNAMIC_GRAPH_REPRESENTATION_LEARNING_WITH_METROPOLIS_HASTINGS_METROPOLIS_HASTINGS_SAMPLERS_H
#define DYNAMIC_GRAPH_REPRESENTATION_LEARNING_WITH_METROPOLIS_HASTINGS_METROPOLIS_HASTINGS_SAMPLERS_H

#include <random_walk_model.h>

namespace dynamic_graph_representation_learning_with_metropolis_hastings
{
    class MetropolisHastingsSampler
    {
        /**
         * @brief MetropolisHastingsSampler constructor.
         */
        MetropolisHastingsSampler()
        {
            this->initialized = false;
            this->last_sampled_vertex = 0;
        }

        void sample(types::State& state, const dygrl::RandomWalkModel* model)
        {
            // 1. Initialize the sampler according to configured initialization strategy
            if(!this->initialized)
            {
//                this->Initialize(state, neighbours, model, config::SAMPLER_INIT_STRATEGY);
//                this->initialized = true;
            }
        }


        private:
            bool initialized;
            types::Vertex last_sampled_vertex;
    };
}

#endif // DYNAMIC_GRAPH_REPRESENTATION_LEARNING_WITH_METROPOLIS_HASTINGS_METROPOLIS_HASTINGS_SAMPLERS_H
