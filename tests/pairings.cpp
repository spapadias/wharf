
#include <gtest/gtest.h>

#include <pairings.h>
#include <pbbslib/parallel.h>
#include <pbbslib/get_time.h>

class PairingsUnitTest : public testing::Test
{
    public:
        void SetUp()    override;
        void TearDown() override;
};

void PairingsUnitTest::SetUp()
{
    std::cout << "-----------------------------------------------------------------------------------------------------"
              << std::endl;
    std::cout << "Dock running with " << num_workers() << " threads" << std::endl;
}

void PairingsUnitTest::TearDown()
{
    std::cout << "-----------------------------------------------------------------------------------------------------"
              << std::endl;
}

TEST_F(PairingsUnitTest, SzudzikPairingFunctionPairs)
{
    auto pair = std::pair<unsigned int, unsigned int>(65535, 65535);

    auto encoded_value = pairings::Szudzik<unsigned int>::pair(pair);
    std::cout << "Encoded value (Szudzik pairing): " << encoded_value << std::endl;
    auto decoded_value = pairings::Szudzik<unsigned int>::unpair(encoded_value);
    std::cout << "Decoded value (Szudzik pairing): " << decoded_value.first << " " << decoded_value.second << std::endl;

    ASSERT_EQ(decoded_value.first, 65535);
    ASSERT_EQ(decoded_value.second, 65535);

    std::cout << std::numeric_limits<uint32_t>::max() << std::endl;
}

TEST_F(PairingsUnitTest, SzudzikPairingFunction10M)
{
    auto array_size = 10000000;
    auto pairs = new std::pair<unsigned int, unsigned int>[array_size];

    parallel_for(0, array_size, [&](size_t i) {
        pairs[i] = std::make_pair(rand() % 65536, rand() % 65536);
    });

    timer timer("Szudzik pairing two elements", false);

    parallel_for(0, array_size, [&](size_t i) {
        timer.start();
        auto encoded = pairings::Szudzik<unsigned int>::pair(pairs[i]);
        auto decoded = pairings::Szudzik<unsigned int>::unpair(encoded);
        timer.stop();

        if (decoded.first != pairs[i].first || decoded.second != pairs[i].second)
        {
            std::cerr << "Unpairing failed: " <<  pairs[i].first << " " << pairs[2].second << std::endl;
        }
    });

    std::cout << "Total time needed to pair and unpair "
              << array_size
              << " pairs is: "
              << timer.get_total()
              << " seconds"
              << std::endl;
}

TEST_F(PairingsUnitTest, SzudzikPairingFunctionTriples)
{
    auto tuple = std::tuple<unsigned int, unsigned int, unsigned int>(123, 25, 200);

    auto encoded_value = pairings::Szudzik<unsigned int>::pair_triplet(tuple);
    std::cout << "Encoded value (Szudzik pairing for triplets): " << encoded_value << std::endl;
    auto decoded_value = pairings::Szudzik<unsigned int>::unpair_triplet(encoded_value);
    std::cout << "Decoded value (Szudzik pairing for triplets): "
              << std::get<0>(decoded_value)
              << " "
              << std::get<1>(decoded_value)
              << " "
              << std::get<2>(decoded_value)
              << std::endl;

    ASSERT_EQ(std::get<0>(decoded_value), 123);
    ASSERT_EQ(std::get<1>(decoded_value), 25);
    ASSERT_EQ(std::get<2>(decoded_value), 200);
}

TEST_F(PairingsUnitTest, SzudzikPairingFunctionTriples10M)
{
    std::cout << "PairingsUnitTest::SzudzikPairingFunctionForTriplets10M" << std::endl;
    auto array_size = 10000000;
    auto pairs = new std::tuple<size_t, size_t, size_t>[array_size];

    parallel_for(0, array_size, [&](size_t i)
    {
        pairs[i] = std::make_tuple(rand() % 65536, rand() % 65536, rand() % 4294967295);
    });

    timer timer("Szudzik pairing three elements", false);

    parallel_for(0, array_size, [&](size_t i)
    {
        timer.start();
        auto encoded = pairings::Szudzik<size_t>::pair_triplet(pairs[i]);
        auto decoded = pairings::Szudzik<size_t>::unpair_triplet(encoded);
        timer.stop();

        if (std::get<0>(decoded) != std::get<0>(pairs[i]) || std::get<1>(decoded) != std::get<1>(pairs[i]) || std::get<2>(decoded) != std::get<2>(pairs[i]))
        {
            std::cerr << "Unpairing triplet failed: " <<  std::get<0>(pairs[i]) << " "
                      << std::get<1>(pairs[i]) << std::get<2>(pairs[i]) << std::endl;
        }
    });

    std::cout << "Total time needed to pair and unpair "
              << array_size
              << " triplets is: "
              << timer.get_total()  << std::endl;
}

