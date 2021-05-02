#include <gtest/gtest.h>
#include <dock.h>

class DockTest : public testing::Test
{
    public:
        void SetUp()    override;
        void TearDown() override;

    protected:
};

void DockTest::SetUp()
{
    std::cout << "SETUP" << std::endl;
}

void DockTest::TearDown()
{
    std::cout << "TEARDOWN" << std::endl;
}

