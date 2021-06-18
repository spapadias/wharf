#ifndef DYNAMIC_GRAPH_REPRESENTATION_LEARNING_WITH_METROPOLIS_HASTINGS_INVERTED_INDEX_H
#define DYNAMIC_GRAPH_REPRESENTATION_LEARNING_WITH_METROPOLIS_HASTINGS_INVERTED_INDEX_H

#include <../compressed_trees/trees/augmented_map.h>

namespace dynamic_graph_representation_learning_with_metropolis_hastings
{
    /**
    * @brief Structure of a walk index entry.
    * Walk entry is a pair ((walk id, vertex position in a walk), next vertex).
    */
    struct WalkIndexEntry
    {
        using key_t = std::pair<types::WalkID, types::Position>; // key: pair(walk id, position in the walk)
        using val_t = types::Vertex;                             // value: the next vertex in the walk
        using aug_t = char;                                      // aug_t: not used but needed for augmented maps

        using entry_t = std::pair<key_t, val_t>;                 // entry: ((walk id, position in the walk), next vertex)

        // key x key -> key
        static bool comp(key_t first_key, key_t second_key)
        {
            return first_key.first == second_key.first ?
                first_key.second < second_key.second : first_key.first < second_key.first;
        }

        // key x value -> augmentation
        static aug_t from_entry(const key_t& key, const val_t& value) { return '\0'; }

        // augmentation x augmentation -> augmentation
        static aug_t combine(const aug_t& augX, const aug_t& augY) { return '\0'; }

        // empty -> augmentation (default augmentation)
        static aug_t get_empty() { return 0; }

        // copy existing entry
        static entry_t copy_entry(const entry_t& entry)
        {
            return std::make_pair(std::make_pair(entry.first.first, entry.first.second), entry.second);
        }

        // delete an entry
        static void del(entry_t& entry) {}
    };

    /**
     * @brief TODO:comment
     */
    class InvertedIndex : public aug_map<WalkIndexEntry>
    {
        public:

        private:
    };

}

#endif //DYNAMIC_GRAPH_REPRESENTATION_LEARNING_WITH_METROPOLIS_HASTINGS_INVERTED_INDEX_H
