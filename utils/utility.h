#ifndef DYNAMIC_GRAPH_REPRESENTATION_LEARNING_WITH_METROPOLIS_HASTINGS_UTILITY_H
#define DYNAMIC_GRAPH_REPRESENTATION_LEARNING_WITH_METROPOLIS_HASTINGS_UTILITY_H

namespace utility
{
    /**
     * Converts B to GB
     *
     * @param bytes
     *
     * @return object size in GB
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
     * Converts B to MB
     *
     * @param bytes
     *
     * @return object size in MB
     */
    double MB(size_t bytes)
    {
        double mb = bytes;
        mb /= 1024.0;
        mb /= 1024.0;

        return mb;
    }

    /**
     * Random number generator.
     */
    class Random
    {
        public:
            // http://xoroshiro.di.unimi.it/#shootout
            // We use xoroshiro128+, the fastest generator available
            uint64_t rng_seed0, rng_seed1;

            Random(uint64_t seed)
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
