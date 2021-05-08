#ifndef DYNAMIC_GRAPH_REPRESENTATION_LEARNING_WITH_METROPOLIS_HASTINGS_METROPOLIS_HASTINGS_SAMPLERS_H
#define DYNAMIC_GRAPH_REPRESENTATION_LEARNING_WITH_METROPOLIS_HASTINGS_METROPOLIS_HASTINGS_SAMPLERS_H

#include <random_walk_model.h>

namespace dynamic_graph_representation_learning_with_metropolis_hastings
{
    class MetropolisHastingsSampler
    {
        public:
             /**
              * @brief MetropolisHastingsSampler constructor.
              *
              * @param current_state - current state of the walker
              * @param model         - model of random walk
              */
            MetropolisHastingsSampler(types::State& current_state, dygrl::RandomWalkModel* model)
            {
                this->last_sampled_vertex = 0;
                this->init(current_state, model, config::sampler_init_strategy);
            }

            /**
             * @brief Sample new vertex.
             *
             * @param state - current walker state
             * @param model - random walk model
             *
             * @return - new state
             */
            types::State sample(types::State& state, dygrl::RandomWalkModel* model)
            {
                // 1. Propose new candidate and calculate weights
                auto candidate_sample = model->propose_vertex(state);
                float new_weight      = model->weight(state, candidate_sample);
                float previous_weight = model->weight(state, this->last_sampled_vertex);

                // 2. Try to accept the candidate
                if (this->accept(previous_weight, new_weight))
                {
                    this->last_sampled_vertex = candidate_sample;
                }

                return model->new_state(state, this->last_sampled_vertex);
            }


        private:
            types::Vertex last_sampled_vertex;

            /**
             * @brief Metropolis Hastings sampler initializer.
             *
             * @param current_state - current state of the walker
             * @param model         - model of random walk
             * @param init_startegy - initialization strategy
             */
            void init(types::State& current_state, dygrl::RandomWalkModel* model, types::SamplerInitStartegy init_startegy)
            {
                /* Random initialization of MH sampler */
                if (init_startegy == types::RANDOM)
                {
                    this->last_sampled_vertex = model->propose_vertex(current_state);
                }
                /* Burn-in initialization of MH sampler */
                else if (init_startegy == types::BURNIN)
                {
                    this->init(current_state, model, types::RANDOM);

                    for (int iter = 0; iter < 100; iter++)
                    {
                        this->last_sampled_vertex = this->sample(current_state, model).first;
                    }
                }
                /* High weight initialization of MH sampler */
                else if (init_startegy == types::WEIGHT)
                {
                    this->init(current_state, model, types::RANDOM);

                    types::Vertex best_sample = this->last_sampled_vertex;
                    float best_sample_weight  = model->weight(current_state, best_sample);

                    for (int iter = 0; iter < 20; iter++)
                    {
                       types::Vertex candidate_sample = model->propose_vertex(current_state);
                       float weight = model->weight(current_state, candidate_sample);

                        if (weight > best_sample_weight)
                        {
                            best_sample_weight = weight;
                            best_sample = candidate_sample;
                        }
                    }

                    this->last_sampled_vertex = best_sample;
                }
            }

            bool accept(float previous_weight, float new_weight)
            {
//                auto random = pbbs::random(time(nullptr));

                if (previous_weight < new_weight) return true;
                return ((double) rand() / (RAND_MAX)) <= (double)(new_weight) / (double)(previous_weight);
            }
    };
}

#endif // DYNAMIC_GRAPH_REPRESENTATION_LEARNING_WITH_METROPOLIS_HASTINGS_METROPOLIS_HASTINGS_SAMPLERS_H
