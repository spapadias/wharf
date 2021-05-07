#ifndef DYNAMIC_GRAPH_REPRESENTATION_LEARNING_WITH_METROPOLIS_HASTINGS_SNAPSHOT_H
#define DYNAMIC_GRAPH_REPRESENTATION_LEARNING_WITH_METROPOLIS_HASTINGS_SNAPSHOT_H

namespace dynamic_graph_representation_learning_with_metropolis_hastings
{
    /**
     * Snapshot is an interface that should be implemented by all types of snapshots
     * (e.g FlatVertexTree, FlatGraph, etc.).
     */
    class Snapshot
    {
        public:
            /**
             * Snapshot size.
             *
             * @return - number of entries in a snapshot
             */
            virtual types::Vertex size() = 0;

            /**
             * Snapshot size in bytes. The memory footprint of a snapshot.
             *
             * @return - size of a snapshot in bytes
             */
            virtual size_t size_in_bytes() = 0;
    };

    /**
     * FlatVertexTree stores vertex entries in a flatten array.
     * This allows for O(1) access to vertex content and promises improved cache locality.
     */
    class FlatVertexTree : public Snapshot
    {
        public:
            /**
             * FlatVertexTree constructor.
             *
             * @param vertices - number of vertices in a graph
             */
            explicit FlatVertexTree(types::Vertex vertices)
            {
                this->snapshot = pbbs::sequence<VertexEntry>(vertices);
            }

            /**
            * FlatVertexTree subscript operator overloading.
            *
            * @param vertex - graph vertex
            *
            * @return - snapshot entry for a given vertex
            */
            VertexEntry& operator[](types::Vertex vertex)
            {
                return this->snapshot[vertex];
            }

            /**
            * FlatVertexTree size.
            *
            * @return - number of entries in a snapshot
            */
            types::Vertex size() final
            {
                return this->snapshot.size();
            }

            /**
            * FlatVertexTree size in bytes. The memory footprint of a snapshot.
            *
            * @return - size of a snapshot in bytes
            */
            size_t size_in_bytes() final
            {
                return this->snapshot.size() * sizeof(this->snapshot[0]);
            }

        private:
            pbbs::sequence<VertexEntry> snapshot;
    };

    /**
     * FlatGraphEntry represents one entry in the flat graph snapshot.
     */
    struct FlatGraphEntry
    {
        types::Vertex* neighbors;
        types::Degree degrees;
        const SamplerManager* samplers;

        /**
         * FlatGraphEntry destructor.
         */
        ~FlatGraphEntry()
        {
            pbbs::free_array(neighbors);
        }
    };

    /**
     * FlatGraph stores for each vertex its neighbors, degree, and reference to its sampler manager.
     * This allows for O(1) access to the vertex, its edges and sampler manager and promises improved cache locality
     * at the expense of memory.
     */
    class FlatGraph : public Snapshot
    {
        public:
            /**
             * FlatGraph constructor.
             *
             * @param vertices - number of vertices in a graph
             */
            explicit FlatGraph(types::Vertex vertices)
            {
                this->snapshot = pbbs::sequence<FlatGraphEntry>(vertices);
            }

            /**
            * FlatGraph subscript operator overloading.
            *
            * @param vertex - graph vertex
            *
            * @return - snapshot entry for a given vertex
            */
            FlatGraphEntry& operator[](types::Vertex vertex)
            {
                return this->snapshot[vertex];
            }

            /**
            * FlatGraph size.
            *
            * @return - number of entries in a snapshot
            */
            types::Vertex size() final
            {
                return this->snapshot.size();
            }

            /**
            * FlatGraph size in bytes. The memory footprint of a snapshot.
            *
            * @return - size of a snapshot in bytes
            */
            size_t size_in_bytes() final
            {
                std::atomic<size_t> total  = 0;

                parallel_for(0, this->snapshot.size(), [&](auto index)
                {
                    total += sizeof(this->snapshot[index]) + (sizeof(this->snapshot[index].neighbors) * this->snapshot[index].degrees);
                });

                return total;
            }

        private:
            pbbs::sequence<FlatGraphEntry> snapshot;
    };
}

#endif // DYNAMIC_GRAPH_REPRESENTATION_LEARNING_WITH_METROPOLIS_HASTINGS_SNAPSHOT_H
