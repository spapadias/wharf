#ifndef DYNAMIC_GRAPH_REPRESENTATION_LEARNING_WITH_METROPOLIS_HASTINGS_RANDOM_WALK_MODEL_H
#define DYNAMIC_GRAPH_REPRESENTATION_LEARNING_WITH_METROPOLIS_HASTINGS_RANDOM_WALK_MODEL_H

namespace dynamic_graph_representation_learning_with_metropolis_hastings
{
    /**
     * @brief RandomWalkModel is the interface that all random walk models (e.g DeepWalk, Node2Vec, etc.) need to implement.
     */
    class RandomWalkModel
    {
        public:
            /**
             * RandomWalkModel constructor.
             *
             * @param snapshot
             */
            explicit RandomWalkModel(const Snapshot* snapshot)
            {
                this->snapshot = snapshot;
            }

            virtual types::State initial_state(dygrl::Vertex vertex) = 0;

            virtual types::State new_state(types::State state, dygrl::Vertex vertex) = 0;

            virtual float weight(types::State state, dygrl::Vertex vertex) = 0;

        private:
            const Snapshot* snapshot = nullptr;
    };
}

#endif // DYNAMIC_GRAPH_REPRESENTATION_LEARNING_WITH_METROPOLIS_HASTINGS_RANDOM_WALK_MODEL_H
