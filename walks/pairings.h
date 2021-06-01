#ifndef DYNAMIC_GRAPH_REPRESENTATION_LEARNING_WITH_METROPOLIS_HASTINGS_PAIRINGS_H
#define DYNAMIC_GRAPH_REPRESENTATION_LEARNING_WITH_METROPOLIS_HASTINGS_PAIRINGS_H

#include <cmath>

namespace pairings
{
    /**
      * @brief Szudzik's pairing function.
      * This function takes two n bit numbers and maps them to a unique 2*n bit number.
      *
      * @details https://www.vertexfragment.com/ramblings/cantor-szudzik-pairing-functions/
      *
      * @tparam Type
      */
    template <class Type>
    class Szudzik
    {
        public:
            /**
             * @brief Maps pair of numbers to a unique number.
             *
             * @param pair - pair of n bit numbers
             *
             * @return - uniquely paired number
             */
            static Type pair(const std::pair<Type, Type>& pair)
            {
                if ( pair.second >= pair.first )
                {
                    #ifdef DOCK_DEBUG
                        Type encoded_value = pair.second * ( pair.second + 1 ) + pair.first;
                        auto decoded_value = unpair(encoded_value);

                        if (decoded_value.first != pair.first || decoded_value.second != pair.second)
                        {
                            std::cerr << "(1) Overflow detected in Szudzik's pairing function. Paired (x="
                                      << pair.first
                                      << ",y="
                                      << pair.second
                                      << ") decoded as (x="
                                      << decoded_value.first
                                      << ",y="
                                      << decoded_value.second
                                      << ")"
                                      << std::endl;

                            std::exit(1);
                        }
                    #endif

                    return pair.second * ( pair.second + 1 ) + pair.first;
                }
                else
                {
                    #ifdef DOCK_DEBUG
                        Type encoded_value = pair.first * pair.first + pair.second;
                        auto decoded_value = unpair(encoded_value);

                        if (decoded_value.first != pair.first || decoded_value.second != pair.second)
                        {
                            std::cerr << "(2) Overflow detected in Szudzik's pairing function. Paired (x="
                                      << pair.first
                                      << ",y="
                                      << pair.second
                                      << ") decoded as (x="
                                      << decoded_value.first
                                      << ",y="
                                      << decoded_value.second
                                      << ")"
                                      << std::endl;

                            std::exit(1);
                        }
                    #endif

                    return pair.first * pair.first + pair.second;
                }
            }

            /**
            * @brief Maps tuple of 3 numbers to a unique number.
            *
            * @param tuple - tuple of 3 n bit numbers
             *
            * @return - uniquely paired number
            */
            static Type pair_triplet(const std::tuple<Type, Type, Type>& tuple)
            {
                return pair({pair({std::get<0>(tuple), std::get<1>(tuple)}), std::get<2>(tuple)});
            };

            /**
             * @brief Converts a unique number to a unique pair of numbers.
             *
             * @param encoded_value - encoded number
             *
             * @return - unique pair of numbers
             */
            static std::pair<Type, Type> unpair(const Type& encoded_value)
            {
                Type sq_z = Type( std::floor(std::sqrt( encoded_value )) );

                if ( sq_z * sq_z > encoded_value )
                    sq_z--;

                Type t = encoded_value - sq_z * sq_z;

                if (t < sq_z)
                    return { sq_z, t };
                else
                    return { t - sq_z, sq_z };
            }

            /**
            * @brief Converts a unique number to a unique tuple of numbers.
            *
            * @param encoded_value - encoded number
            *
            * @return - unique tuple of numbers
            */
            static std::tuple<Type, Type, Type> unpair_triplet(const Type& encoded_value)
            {
                auto first_pair  = unpair(encoded_value);
                auto second_pair = unpair(first_pair.first);

                return std::make_tuple(second_pair.first, second_pair.second, first_pair.second);
            }
    };
}

#endif // DYNAMIC_GRAPH_REPRESENTATION_LEARNING_WITH_METROPOLIS_HASTINGS_PAIRINGS_H
