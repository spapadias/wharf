#ifndef DEEPWALK_H
#define DEEPWALK_H

#include <random_walk_model.h>

namespace dynamic_graph_representation_learning_with_metropolis_hastings
{
    /**
     * @brief DeepWalk random walk model implementation.
     * @details https://dl.acm.org/doi/abs/10.1145/2623330.2623732
     */
    class DeepWalk : public RandomWalkModel
    {
        public:
            /**
             * @brief DeepWalk constructor.
             *
             * @param snapshot - graph snapshot
             */
            explicit DeepWalk(dygrl::Snapshot* snapshot)
            {
                this->snapshot = snapshot;
            }

            /**
            * @brief DeepWalk destructor.
            */
            ~DeepWalk()
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
                return types::State(vertex, vertex);
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
                return types::State(vertex, vertex);
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
                return 1.0;
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
    };
}

#endif
