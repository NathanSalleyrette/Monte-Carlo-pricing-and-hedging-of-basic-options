#include "pnl/pnl_random.h"
#include "pnl/pnl_vector.h"
#include "OptionBasket.cpp"
#include "MonteCarlo.cpp"
#include "BlackScholesModel.cpp"
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

TEST(OptionBasket, CalculPrix) {
  // Initialisation du randomisateur
    PnlRng *rng = pnl_rng_create(PNL_RNG_MERSENNE);
    pnl_rng_sseed(rng, 0);

    PnlVect *G = pnl_vect_create(3);
    G->array[0] = 0.5;
    G->array[1] = 0.4;
    G->array[2] = 0.1;
    OptionBasket opt = OptionBasket(1.0, 252, G->size, G, 10.0);

      // Initialisation Objet BlackScholes
    PnlVect *Sigma = pnl_vect_create_from_list(3, 0.5, 0.5, 0.5);
    PnlVect *Spot = pnl_vect_create_from_list(3, 8.0, 12.0, 15.0);
    BlackScholesModel bl = BlackScholesModel(3, 0.02, 0.15, Sigma, Spot);

    // Initialisation Objet MonteCarlo
    double prix;
    double std_dev;
    MonteCarlo mtc = MonteCarlo(&bl, &opt, rng, 1.0, 1000);
    mtc.price(prix, std_dev);
    
    
    EXPECT_NE(prix, 0.0);
    EXPECT_NE(std_dev, 0.0);
    EXPECT_DOUBLE_EQ(prix, 1.5070789066119803);
    EXPECT_DOUBLE_EQ(std_dev, 0.17162142732337202);

    pnl_vect_free(&G);
    pnl_vect_free(&Sigma);
    pnl_vect_free(&Spot);
    pnl_rng_free(&rng);


}
}