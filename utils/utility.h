#ifndef DYNAMIC_GRAPH_REPRESENTATION_LEARNING_WITH_METROPOLIS_HASTINGS_UTILITY_H
#define DYNAMIC_GRAPH_REPRESENTATION_LEARNING_WITH_METROPOLIS_HASTINGS_UTILITY_H

#include <rmat_util.h>

namespace utility
{
    /**
     * @brief Converts B to GB.
     *
     * @param bytes - size in B
     *
     * @return size in GB
     */
    double GB(size_t bytes)
    {
        double gb = bytes;
        gb /= 1024.0;
        gb /= 1024.0;
        gb /= 1024.0;

        return gb;
    }

    /**
     * @brief Converts B to MB.
     *
     * @param bytes - size in B
     *
     * @return size in MB
     */
    double MB(size_t bytes)
    {
        double mb = bytes;
        mb /= 1024.0;
        mb /= 1024.0;

        return mb;
    }

    /**
    * @brief Generates batch of edges for insertion or deletion.
    *
    * @param edges_number - number of edges to generate
    * @param vertices_number - number of vertices in the graph
    * @param self_loops - determines whether self loops are allowed
    * @param directed - determines whether it needs to generate directed or undirected edges
    * @param a
    * @param b
    * @param c
    * @param run_seq - determines whether to run part of the code sequantially
    *
    * @return - batch of generated edges
    */
    auto generate_batch_of_edges
    (
        size_t edges_number,
        size_t vertices_number,
        bool self_loops = false,
        bool directed = true,
        double a = 0.5,
        double b = 0.2,
        double c = 0.1,
        bool run_seq = false
    )
    {
        #ifdef MALIN_TIMER
            timer timer("Utility::GenerateBatchOfEdges");
        #endif

        using Edge = std::tuple<uintV, uintV>;  // (vertex1, vertex2)

        // 1. Set up for the edge generation
        auto rand              = pbbs::random(0);
//        auto rand              = pbbs::random(std::time(nullptr));
        size_t graph_size_pow2 = 1 << (pbbs::log2_up(vertices_number) - 1);
        auto rmat              = rMat<uintV>(graph_size_pow2, rand.ith_rand(0), a, b, c);
        auto edges             = (directed) ? pbbs::new_array<Edge>(edges_number) : pbbs::new_array<Edge>(2 * edges_number);

        // 2. Generate edges in parallel
        parallel_for(0, edges_number, [&] (size_t i)
        {
            edges[i] = rmat(i);
        });

        if (!directed)
        {
            parallel_for(0, edges_number, [&] (size_t i)
            {
                edges[i + edges_number] = std::make_tuple(std::get<1>(edges[i]), std::get<0>(edges[i]));;
            });
        }

        // 3. Sort edges by source
        auto edges_sorted = (directed) ? pbbs::make_range(edges, edges + edges_number) : pbbs::make_range(edges, edges + 2 * edges_number);
        auto node_bits    = pbbs::log2_up(graph_size_pow2);

        auto edge_to_long = [graph_size_pow2, node_bits](Edge e) -> size_t {
            return (static_cast<size_t>(std::get<0>(e)) << node_bits) + static_cast<size_t>(std::get<1>(e));
        };

        size_t bits = 2 * node_bits;

        if (graph_size_pow2 <= (edges_sorted.size() * pbbs::log2_up(edges_sorted.size())))
        {
            pbbs::integer_sort_inplace(edges_sorted, edge_to_long, bits);
        }
        else
        {
            pbbs::sample_sort_inplace(edges_sorted, std::less<>());
        }

        // 4. Remove duplicate edges
        Edge* generated_edges = nullptr;

        // Remove duplicate edges
        auto bool_seq = pbbs::delayed_seq<bool>(edges_sorted.size(), [&] (size_t i)
        {
            if (!self_loops && std::get<0>(edges_sorted[i]) == std::get<1>(edges_sorted[i])) return false;
            return (i == 0 || edges_sorted[i] != edges_sorted[i - 1]);
        });

        auto pack = pbbs::pack(edges_sorted, bool_seq, run_seq ? pbbs::fl_sequential : pbbs::no_flag);
        auto edges_generated = pack.size();
        generated_edges = pack.to_array();

        pbbs::free_array(edges);

        #ifdef MALIN_DEBUG
            std::cout << edges_generated << " / "
                      << ((directed) ? edges_number : 2 * edges_number)
                      << " distinct edges (direction, duplicate removal and self-loop processing)"
                      << std::endl;
        #endif

        #ifdef MALIN_TIMER
            timer.reportTotal("time(seconds)");
        #endif

        return std::make_pair(generated_edges, edges_generated);
    }

    /**
     * @brief Random number generator.
     * @details http://xoroshiro.di.unimi.it/#shootout
     */
    class Random
    {
        public:
            uint64_t rng_seed0, rng_seed1;

            explicit Random(uint64_t seed)
            {
                for (int i = 0; i < 2; i++)
                {
                    long long z = seed += UINT64_C(0x9E3779B97F4A7C15);

                    z = (z ^ z >> 30) * UINT64_C(0xBF58476D1CE4E5B9);
                    z = (z ^ z >> 27) * UINT64_C(0x94D049BB133111EB);

                    if (i == 0)
                        rng_seed0 = z ^ (z >> 31);
                    else
                        rng_seed1 = z ^ (z >> 31);
                }
            }

            void reinit(uint64_t seed)
            {
                for (int i = 0; i < 2; i++)
                {
                    long long z = seed += UINT64_C(0x9E3779B97F4A7C15);

                    z = (z ^ z >> 30) * UINT64_C(0xBF58476D1CE4E5B9);
                    z = (z ^ z >> 27) * UINT64_C(0x94D049BB133111EB);

                    if (i == 0)
                        rng_seed0 = z ^ (z >> 31);
                    else
                        rng_seed1 = z ^ (z >> 31);
                }
            }

            static inline uint64_t rotl(const uint64_t x, int k)
            {
                return (x << k) | (x >> (64 - k));
            }

            uint64_t lrand()
            {
                const uint64_t s0 = rng_seed0;
                uint64_t s1 = rng_seed1;

                const uint64_t result = s0 + s1;
                s1 ^= s0;

                rng_seed0 = rotl(s0, 55) ^ s1 ^ (s1 << 14);
                rng_seed1 = rotl(s1, 36);

                return result;
            }

            double drand()
            {
                const union un
                {
                    uint64_t i;
                    double d;
                }
                a = {UINT64_C(0x3FF) << 52 | lrand() >> 12};

                return a.d - 1.0;
            }

            int irand(int max) { return lrand() % max; }

            int irand(int min, int max) { return lrand() % (max - min) + min; }
    };
}

#endif // DYNAMIC_GRAPH_REPRESENTATION_LEARNING_WITH_METROPOLIS_HASTINGS_UTILITY_H
