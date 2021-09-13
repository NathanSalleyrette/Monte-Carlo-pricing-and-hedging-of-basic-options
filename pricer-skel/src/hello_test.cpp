#include <gtest/gtest.h>
// #include "pnl/pnl_random.h"
// #include "pnl/pnl_vector.h"
#include "OptionBasket.cpp"

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

// TEST(OptionBasket, Payoff) {
//     PnlVect *G = pnl_vect_new();
//     OptionBasket opt = OptionBasket(1.0, 1, 1, G, 10.0);
//     PnlMat *M = pnl_mat_new();
//     EXPECT_DOUBLE_EQ(opt.payoff(M), 0.0);
// }