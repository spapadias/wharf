#ifndef DYNAMIC_GRAPH_REPRESENTATION_LEARNING_WITH_METROPOLIS_HASTINGS_RANDOM_WALK_MODEL_H
#define DYNAMIC_GRAPH_REPRESENTATION_LEARNING_WITH_METROPOLIS_HASTINGS_RANDOM_WALK_MODEL_H

namespace dynamic_graph_representation_learning_with_metropolis_hastings
{
    /**
    * RandomWalkModel is the interface that all random walk models (e.g DeepWalk, Node2Vec, etc.) need to implement.
    */
    class RandomWalkModel
    {
        public:
            /**
            * @brief Determines an initial state of the walker.
            *
            * @param vertex - graph vertex
            *
            * @return - an initial state of the walker
            */
            virtual types::State initial_state(types::Vertex vertex) = 0;

            /**
            * @brief The transition of states is crucial for the walking.
            * Based on the next vertex and the current state we update the state.
            *
            * @param state  - current state of the walker
            * @param vertex - next vertex to go
            *
            * @return - a new state of the walker
            */
            virtual types::State update_state(const types::State& state, types::Vertex vertex) = 0;

            /**
            * @brief Explains how to calculate the edge weight based on the current state and the potentially proposed vertex.
            *
            * @param state  - current state of the walker
            * @param vertex - potentially proposed vertex
            *
            * @return float - dynamically calculated weight
            */
            virtual float weight(const types::State& state, types::Vertex vertex) = 0;
    };
}

#endif // DYNAMIC_GRAPH_REPRESENTATION_LEARNING_WITH_METROPOLIS_HASTINGS_RANDOM_WALK_MODEL_H
