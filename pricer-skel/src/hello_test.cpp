#include "pnl/pnl_random.h"
#include "pnl/pnl_vector.h"
#include "OptionBasket.cpp"
#include <gtest/gtest.h>

namespace {

// Demonstrate some basic assertions.
TEST(HelloTest, BasicAssertions) {
  // Expect two strings not to be equal.
  EXPECT_STRNE("hello", "world");
  // Expect equality.
  EXPECT_EQ(7 * 6, 42);
}

TEST(HelloTest, MoreAssertions) {
    EXPECT_FALSE(false);
    EXPECT_TRUE(true);
}

TEST(OptionBasket, Entrecote) {
    PnlVect *G = pnl_vect_create(3);
    G->array[0] = 0.5;
    G->array[1] = 0.4;
    G->array[2] = 0.1;
    OptionBasket opt = OptionBasket(1.0, 1, 1, G, 1.0);
    PnlMat *M = pnl_mat_create(1,3);
    M->array[0] = 6.0;
    M->array[1] = 3.0;
    M->array[2] = 2.0;
    EXPECT_DOUBLE_EQ(opt.payoff(M), 0.0);
}
}