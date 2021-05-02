#ifndef DYNAMIC_GRAPH_REPRESENTATION_LEARNING_WITH_METROPOLIS_HASTINGS_UTILITY_H
#define DYNAMIC_GRAPH_REPRESENTATION_LEARNING_WITH_METROPOLIS_HASTINGS_UTILITY_H

namespace utility
{
    /**
     * Converts B to GB
     *
     * @param bytes
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
     * @return object size in MB
     */
    double MB(size_t bytes)
    {
        double mb = bytes;
        mb /= 1024.0;
        mb /= 1024.0;

        return mb;
    }
}

#endif // DYNAMIC_GRAPH_REPRESENTATION_LEARNING_WITH_METROPOLIS_HASTINGS_UTILITY_H
