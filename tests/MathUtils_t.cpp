#include <gtest/gtest.h>
#include <Matriz.hpp>
#include <MathUtils.hpp>
#include <limits>

TEST(MathUtils, Identity)
{
    const auto a = MatrizI(3);
    Matriz b{3, 3, {1, 2, 3, 4, 5, 6, 7, 8, 9}};

    EXPECT_EQ(a * b, b);
}

TEST(MathUtils, MatrizZeros)
{
    const auto a = MatrizZeros(3, 3);
    Matriz b{3, 3, {1, 2, 3, 4, 5, 6, 7, 8, 9}};

    EXPECT_EQ(a * b, a);
}

TEST(MathUtils, Inversa)
{
    Matriz a{3, 3, 
        {1, 5, 7, 
        5, 7, 9,
        2, 4, 5}
    };
    Matriz b = CalcInvMatriz(a);

    const auto err = (b*a) - MatrizI(3);

    MostrarMatriz(b);
    MostrarMatriz(a);
    MostrarMatriz(err);

    for(size_t i = 0; i < 3; i++)
    {
        for(size_t j = 0; j < 3; j++)
        {
            EXPECT_NEAR(err(i, j), 0.0, std::numeric_limits<decltype(err(i,j))>::epsilon()*10);
        }
    }
}
