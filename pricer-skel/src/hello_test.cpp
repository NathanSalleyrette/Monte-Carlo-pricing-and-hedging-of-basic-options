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

TEST(OptionBasket, Payoff) {
    double nbSample = 1;
    PnlVect *G = pnl_vect_create(3);
    G->array[0] = 0.5;
    G->array[1] = 0.4;
    G->array[2] = 0.1;
    OptionBasket opt = OptionBasket(1.0, nbSample, G->size, G, 1.0);
    PnlMat *M = pnl_mat_create(G->size,nbSample+1);
    M->array[0] = 6.0;
    M->array[1] = 3.0;
    M->array[2] = 2.0;
    M->array[3] = 6.0;
    M->array[4] = 3.0;
    M->array[5] = 2.0;
    EXPECT_DOUBLE_EQ(opt.payoff(M), 3.1000000000000005);

    // On libère la mémoire
    pnl_vect_free(&G);
    pnl_mat_free(&M);
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

    PnlVect *sumShift = pnl_vect_create(3);
    PnlVect *sumShiftSquare = pnl_vect_create(3);

    // Initialisation Objet MonteCarlo
    double prix;
    double std_dev;
    MonteCarlo mtc = MonteCarlo(&bl, &opt, rng, 1.0, 1000, sumShift, sumShiftSquare);
    mtc.price(prix, std_dev);
    
    PnlVect *deltas = pnl_vect_create(3);
    PnlVect *std_dev_d = pnl_vect_create(3);

    mtc.delta(deltas, std_dev_d); 
    pnl_vect_print(deltas);
    // On teste
    EXPECT_NE(prix, 0.0);
    EXPECT_NE(std_dev, 0.0);
    EXPECT_DOUBLE_EQ(prix, 1.5070789066119803);
    EXPECT_DOUBLE_EQ(std_dev, 0.17162142732337202);

    // pnl_rng_sseed(rng, 1);
    // mtc.price(prix, std_dev);
    // EXPECT_DOUBLE_EQ(prix, 1.5070789066119803);
    // EXPECT_DOUBLE_EQ(std_dev, 0.17162142732337202);



    // On libère 
    pnl_vect_free(&G);
    pnl_vect_free(&Sigma);
    pnl_vect_free(&Spot);
    pnl_rng_free(&rng);


}

TEST(OptionBasket, CalculPrixInstantT) {
  // Initialisation du randomisateur
    PnlRng *rng = pnl_rng_create(PNL_RNG_MERSENNE);
    pnl_rng_sseed(rng, 0);

    PnlVect *G = pnl_vect_create(3);
    G->array[0] = 0.5;
    G->array[1] = 0.4;
    G->array[2] = 0.1;
    OptionBasket opt = OptionBasket(1.0, 2, G->size, G, 10.0);

      // Initialisation Objet BlackScholes
    PnlVect *Sigma = pnl_vect_create_from_list(3, 0.5, 0.5, 0.5);
    PnlVect *Spot = pnl_vect_create_from_list(3, 8.0, 12.0, 15.0);
    PnlMat *Past = pnl_mat_create_from_list(3,2, 8.0, 9.0, 12.0, 11.0, 15.0, 15.5);
    BlackScholesModel bl = BlackScholesModel(3, 0.02, 0.15, Sigma, Spot);

    PnlVect *sumShift = pnl_vect_create(3);
    PnlVect *sumShiftSquare = pnl_vect_create(3);

    // Initialisation Objet MonteCarlo
    double prix = 0.0;
    double std_dev = 0.0;
    MonteCarlo mtc = MonteCarlo(&bl, &opt, rng, 1.0, 1000, sumShift, sumShiftSquare);
    mtc.price(Past, 0.3, prix, std_dev);
    
    PnlVect *deltas = pnl_vect_create(3);
    PnlVect *std_dev_d = pnl_vect_create(3);

    mtc.delta(deltas, std_dev_d); 
    pnl_vect_print(deltas);
    // On teste
    EXPECT_NE(prix, 0.0);
    EXPECT_NE(std_dev, 0.0);
    EXPECT_DOUBLE_EQ(prix, 1.5070789066119803);
    EXPECT_DOUBLE_EQ(std_dev, 0.17162142732337202);

    // pnl_rng_sseed(rng, 1);
    // mtc.price(prix, std_dev);
    // EXPECT_DOUBLE_EQ(prix, 1.5070789066119803);
    // EXPECT_DOUBLE_EQ(std_dev, 0.17162142732337202);



    // On libère 
    pnl_vect_free(&G);
    pnl_vect_free(&Sigma);
    pnl_vect_free(&Spot);
    pnl_rng_free(&rng);


}
}