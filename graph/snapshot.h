#ifndef SNAPSHOT_H
#define SNAPSHOT_H

namespace dynamic_graph_representation_learning_with_metropolis_hastings
{
    /**
     * @brief Snapshot is an interface that should be implemented by all types of snapshots
     * (e.g FlatVertexTree, FlatGraph, etc.).
     */
    class Snapshot
    {
        public:
            /**
            * @brief Given a graph vertex it returns its neighbors and degree.
            *
            * @param vertex - graph vertex
            *
            * @return - neighbors and degree of a given vertex
            */
            virtual types::Neighbors neighbors(types::Vertex vertex) = 0;

            /**
             * @brief Snapshot size.
             *
             * @return - number of entries in a snapshot
             */
            virtual types::Vertex size() = 0;

            /**
             * @brief Snapshot size in bytes. The memory footprint of a snapshot.
             *
             * @return - size of a snapshot in bytes
             */
            virtual size_t size_in_bytes() = 0;
    };

    /**
     * @brief FlatVertexTree stores vertex entries in an array.
     * This allows for O(1) access to vertex content and promises improved cache locality.
     */
    class FlatVertexTree : public Snapshot
    {
        public:
            /**
             * @brief FlatVertexTree constructor.
             *
             * @param vertices - number of vertices in a graph
             */
            explicit FlatVertexTree(types::Vertex vertices)
            {
                this->snapshot = pbbs::sequence<VertexEntry>(vertices);
            }

            /**
            * @brief FlatVertexTree subscript operator overloading.
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
            * @brief FlatVertexTree size.
            *
            * @return - number of entries in a snapshot
            */
            types::Vertex size() final
            {
                return this->snapshot.size();
            }

            /**
            * @brief FlatVertexTree size in bytes. The memory footprint of a snapshot.
            *
            * @return - size of a snapshot in bytes
            */
            size_t size_in_bytes() final
            {
                return this->snapshot.size() * sizeof(this->snapshot[0]);
            }

            /**
            * @brief Given a graph vertex it returns its neighbors and degree.
            *
            * @param vertex - graph vertex
            *
            * @return - neighbors and degree of a given vertex
            */
            types::Neighbors neighbors(types::Vertex vertex) final
            {
                auto neighbors = snapshot[vertex].compressed_edges.get_edges(vertex);
                auto degrees   = snapshot[vertex].compressed_edges.degree();

                return std::make_tuple(neighbors, degrees, true);
            }

        private:
            pbbs::sequence<VertexEntry> snapshot;
    };

    /**
     * @brief FlatGraphEntry represents one entry in the flat graph snapshot.
     */
    struct FlatGraphEntry
    {
        types::Vertex* neighbors;
        types::Degree degree;
        SamplerManager* samplers;

        /**
         * @brief FlatGraphEntry destructor.
         */
        ~FlatGraphEntry()
        {
            pbbs::free_array(neighbors);
        }
    };

    /**
     * @brief FlatGraph stores for each vertex its neighbors, degree, and reference to its sampler manager.
     * This allows for O(1) access to the vertex, its edges and sampler manager and promises improved cache locality
     * at the expense of memory.
     */
    class FlatGraph : public Snapshot
    {
        public:
            /**
             * @brief FlatGraph constructor.
             *
             * @param vertices - number of vertices in a graph
             */
            explicit FlatGraph(types::Vertex vertices)
            {
                this->snapshot = pbbs::sequence<FlatGraphEntry>(vertices);
            }

            /**
            * @brief FlatGraph subscript operator overloading.
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
            * @brief FlatGraph size.
            *
            * @return - number of entries in a snapshot
            */
            types::Vertex size() final
            {
                return this->snapshot.size();
            }

            /**
            * @brief FlatGraph size in bytes. The memory footprint of a snapshot.
            *
            * @return - size of a snapshot in bytes
            */
            size_t size_in_bytes() final
            {
                std::atomic<size_t> total  = 0;

                parallel_for(0, this->snapshot.size(), [&](auto index)
                {
                    total += sizeof(this->snapshot[index]) + (sizeof(this->snapshot[index].neighbors) * this->snapshot[index].degree);
                });

                return total;
            }

            /**
            * @brief Returns neighbors of a given vertex.
            *
            * @param vertex - graph vertex
            *
            * @return - neighbors of a given vertex
            */
            types::Neighbors neighbors(types::Vertex vertex) final
            {
                return std::make_tuple(snapshot[vertex].neighbors, snapshot[vertex].degree, false);
            }

        private:
            pbbs::sequence<FlatGraphEntry> snapshot;
    };
}

#endif
