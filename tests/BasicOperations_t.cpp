#include <gtest/gtest.h>
#include <Matriz.hpp>

TEST(BasicOperations, Sum)
{
    Matriz a{1, 3, {1, 2, 3}};
    Matriz b{1, 3, {1, 2, 3}};
    Matriz c{1, 3, {2, 4, 6}};

    EXPECT_EQ(a + b, c);
}

TEST(BasicOperations, Subtract)
{
    Matriz a{1, 3, {1, 2, 3}};
    Matriz b{1, 3, {1, 2, 3}};
    Matriz c{1, 3, {0, 0, 0}};

    EXPECT_EQ(a - b, c);
}

TEST(BasicOperations, Multiply)
{
    Matriz a{1, 3, {1, 2, 3}};
    Matriz b{3, 3, {1, 2, 3, 4, 5, 6, 7, 8, 9}};
    Matriz c{1, 3, {30, 36, 42}};

    EXPECT_EQ(a * b, c);
}

