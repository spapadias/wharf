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
        using key_t = uintV;                       // key: pair(walk id, position in the walk)
        using val_t = types::Vertex;               // value: the next vertex in the walk
        using aug_t = char;                        // aug_t: not used but needed for augmented maps

        using entry_t = std::pair<key_t, val_t>;   // entry: ((walk id, position in the walk), next vertex)

        // key x key -> key
        static bool comp(const key_t& keyX, const key_t& keyY) { return keyX < keyY; }

        // key x value -> augmentation
        static aug_t from_entry(const key_t& key, const val_t& value) { return '\0'; }

        // augmentation x augmentation -> augmentation
        static aug_t combine(const aug_t& augX, const aug_t& augY) { return '\0'; }

        // empty -> augmentation (default augmentation)
        static aug_t get_empty() { return 0; }

        // copy existing entry
        static entry_t copy_entry(const entry_t& entry) { return std::make_pair(entry.first, entry.second); }

        // delete an entry
        static void del(entry_t& entry) {}
    };

    /**
     * @brief Inverted Index as an augmented parallel balanced binary tree. It stores parts of the walk in the form
     * of a pair ((walk id, vertex position in a walk), next vertex).
     */
    class InvertedIndex : public aug_map<WalkIndexEntry>
    {
        public:

            /**
             * @brief InvertedIndex default constructor.
             */
            InvertedIndex() : aug_map<WalkIndexEntry>() {};

            /**
             * @brief InvertedIndex constructor.
             *
             * @param walk_index_entries - inverted index entries
             */
            explicit InvertedIndex(const pbbs::sequence<dygrl::WalkIndexEntry::entry_t>& walk_index_entries) : aug_map<WalkIndexEntry>(walk_index_entries) {};

            /**
             * @brief InvertedIndex copy constructor.
             *
             * @param other - inverted index
             */
            explicit InvertedIndex(const aug_map<WalkIndexEntry>& other) : aug_map<WalkIndexEntry>(other) {};

            /**
            * @brief Finds the next vertex in the walk given walk id and position
            *
            * @param walk_id  - unique walk id
            * @param position - position in the walk
            *
            * @return - next vertex in the walk
            */
            types::Vertex find_next(types::WalkID walk_id, types::Position position)
            {
                auto result = this->find(walk_id*config::walk_length + position);

                #ifdef MALIN_DEBUG
                    if (!result.valid)
                    {
                        std::cerr << "Malin debug error! InvertedIndex::FindNext::walk_id = "
                                  << walk_id
                                  << ", position = "
                                  << (int) position
                                  << std::endl;

                        std::exit(1);
                    }
                #endif

                return result.value;
            }
    };

}

#endif // DYNAMIC_GRAPH_REPRESENTATION_LEARNING_WITH_METROPOLIS_HASTINGS_INVERTED_INDEX_H
