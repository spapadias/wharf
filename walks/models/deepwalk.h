#ifndef DYNAMIC_GRAPH_REPRESENTATION_LEARNING_WITH_METROPOLIS_HASTINGS_DEEPWALK_H
#define DYNAMIC_GRAPH_REPRESENTATION_LEARNING_WITH_METROPOLIS_HASTINGS_DEEPWALK_H

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
            explicit DeepWalk(const dygrl::Snapshot* snapshot)
            {
                this->snapshot = snapshot;
            };

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
            types::State update_state(const types::State& state, types::Vertex vertex) final
            {
                return types::State(vertex, vertex);
            }

            /**
            * @brief Explains how to calculate the edge weight based on the current state and the potentially proposed vertex.
            *
            * @param state  - current state of the walker
            * @param vertex - potentially proposed vertex
            *
            * @return float - dynamically calculated weight
            */
            float weight(const types::State& state, types::Vertex vertex) final
            {
                return 1.0;
            }

        private:
            const Snapshot* snapshot;
    };
}

#endif // DYNAMIC_GRAPH_REPRESENTATION_LEARNING_WITH_METROPOLIS_HASTINGS_DEEPWALK_H
