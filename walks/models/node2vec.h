#ifndef NODE2VEC_H
#define NODE2VEC_H

#include <random_walk_model.h>

namespace dynamic_graph_representation_learning_with_metropolis_hastings
{
    /**
     * @brief Node2Vec random walk model implementation.
     * @details https://arxiv.org/abs/1607.00653
     */
    class Node2Vec : public RandomWalkModel
    {
        public:
            /**
             * @brief Node2Vec constructor.
             *
             * @param snapshot - graph snapshot
             */
            explicit Node2Vec(dygrl::Snapshot* snapshot, float paramP, float paramQ)
            {
                this->paramP = paramP;
                this->paramQ = paramQ;
                this->snapshot = snapshot;
            }

            /**
            * @brief Node2Vec destructor.
            */
            ~Node2Vec()
            {
                delete this->snapshot;
            }

            /**
            * @brief Determines an initial state of the walker.
            *
            * @param vertex - graph vertex
            *
            * @return - an initial state of the walker
            */
            types::State initial_state(types::Vertex vertex) final
            {
                auto neighbors   = this->snapshot->neighbors(vertex);
                auto prev_vertex = std::get<0>(neighbors)[config::random.irand(std::get<1>(neighbors))];

                if (std::get<2>(neighbors)) pbbs::free_array(std::get<0>(neighbors));

                return types::State(vertex, prev_vertex);
            }

            /**
            * @brief The transition of states is crucial for the walking.
            * Based on the next vertex and the current state we update the state.
            *
            * @param state  - current state of the walker
            * @param vertex - next vertex to go
            *
            * @return - a new state of the walker
            */
            types::State new_state(const types::State& state, types::Vertex vertex) final
            {
                return types::State(vertex, state.first);
            }

            /**
            * @brief Calculates the edge weight based on the current state and the potentially proposed vertex.
            *
            * @param state  - current state of the walker
            * @param vertex - potentially proposed vertex
            *
            * @return - dynamically calculated weight
            */
            float weight(const types::State& state, types::Vertex vertex) final
            {
                if (vertex == state.second)
                {
                    return 1 / this->paramP;
                }
                else if (this->has_edge(state.second, vertex))
                {
                    return 1;
                }
                else
                {
                    return 1 / this->paramQ;
                }
            }

            /**
            * @brief Propose next vertex given current state.
            *
            * @param vertex - current walker state
            *
            * @return - proposed vertex
            */
            types::Vertex propose_vertex(const types::State& state) final
            {
                auto neighbors = this->snapshot->neighbors(state.first);
                auto vertex    = std::get<0>(neighbors)[config::random.irand(std::get<1>(neighbors))];

                if (std::get<2>(neighbors)) pbbs::free_array(std::get<0>(neighbors));

                return vertex;
            }

        private:
            Snapshot* snapshot;
            float paramP;
            float paramQ;

            bool has_edge(types::Vertex prev, types::Vertex next)
            {
                auto neighbors = this->snapshot->neighbors(prev);
                bool res = binary_search(std::get<0>(neighbors), std::get<0>(neighbors) + std::get<1>(neighbors), next);
                if (std::get<2>(neighbors)) pbbs::free_array(std::get<0>(neighbors));

                return res;
            }
    };
}

#endif
