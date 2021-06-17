#ifndef DYNAMIC_GRAPH_REPRESENTATION_LEARNING_WITH_METROPOLIS_HASTINGS_INVERTED_INDEX_H
#define DYNAMIC_GRAPH_REPRESENTATION_LEARNING_WITH_METROPOLIS_HASTINGS_INVERTED_INDEX_H

#include <../compressed_trees/trees/augmented_map.h>

namespace dynamic_graph_representation_learning_with_metropolis_hastings
{
    class InvertedIndex
    {
        public:

            /**
            * @brief Structure of a walk entry.
            * Walk entry is a pair <<walk id, current vertex position in a walk>, next vertex>.
            */
            struct WalkEntry
            {
                using key_t   = std::pair<types::WalkID, types::Position>;  // key: pair(walk id, position in the walk)
                using val_t   = uintV;                                      // value: the next vertex in the walk
                using aug_t   = char;                                       // aug_t: dummy

                using entry_t = std::pair<key_t, val_t>;                    // entry: <pair(walk id, position in the walk), next vertex>

                // key x key -> key
                static bool comp(key_t first_key, key_t second_key)
                {
                    return first_key.first == second_key.first ?
                        first_key.second < second_key.second : first_key.first < second_key.first;
                }

                // key x value -> augmentation
                static aug_t from_entry(const key_t& key, const val_t& value)
                {
                    return 'a';
                }

                // augmentation x augmentation -> augmentation
                static aug_t combine(const aug_t& augX, const aug_t& augY)
                {
                    return 'a';
                }

                // empty -> augmentation (default augmentation)
                static aug_t get_empty()
                {
                    return 0;
                }

                // copy existing entry
                static entry_t copy_entry(const entry_t& entry)
                {
                    return std::make_pair(std::make_pair(entry.first.first, entry.first.second), entry.second);
                }

                // delete an entry
                static void del(entry_t& entry) {}
            };

            using Index = aug_map<WalkEntry>;  // all random walk pairs for one vertex

            /**
            * Inverted index object
            */
            Index proxy;

            /**
             * @brief Inverted index default constructor.
             */
            InvertedIndex()
            {
                this->proxy.root = NULL;
                Index::GC::init();
            }

            /**
            * @brief Inverted index constructor.
            *
            * @param S sequence of std::pair<unsigned int, walk_triplet> to insert into a index
            */
            explicit InvertedIndex(const pbbs::sequence<WalkEntry::entry_t> &walk_entries)
            {
                this->proxy = Index::multi_insert(this->proxy.root, walk_entries);
            }

            /**
             * @brief Inverted index copy constructor.
             *
             * @param inverted_index inverted index to copy
             */
            InvertedIndex(const InvertedIndex& inverted_index)
            {
                this->proxy.root = inverted_index.proxy.root;
                Index::GC::increment(this->proxy.root);
            }

            /**
            * @brief Inverted index copy constructor.
            *
            * @param inverted_index inverted index to copy
            */
            InvertedIndex(const Index& index)
            {
                this->proxy.root = index.root;
                Index::GC::increment(this->proxy.root);
            }
    };

}

#endif //DYNAMIC_GRAPH_REPRESENTATION_LEARNING_WITH_METROPOLIS_HASTINGS_INVERTED_INDEX_H
