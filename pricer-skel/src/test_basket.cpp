#include "pnl/pnl_random.h"
#include "pnl/pnl_vector.h"
#include "OptionBasket.cpp"
#include "OptionAsian.cpp"
#include "MonteCarlo.cpp"
#include "BlackScholesModel.cpp"
#include <gtest/gtest.h>

namespace {



TEST(CalculPrixInstant0, OptionBasket1) {
  // Initialisation du randomisateur
    PnlRng *rng = pnl_rng_create(PNL_RNG_MERSENNE);
    pnl_rng_sseed(rng, time(NULL));

    PnlVect *G = pnl_vect_create_from_scalar(40, 0.025);

    OptionBasket opt = OptionBasket(3, 1, G->size, G, 100.0);

      // Initialisation Objet BlackScholes
    PnlVect *Sigma = pnl_vect_create_from_scalar(40, 0.2);
    PnlVect *Spot = pnl_vect_create_from_scalar(40, 100.0);
    //PnlMat *Past = pnl_mat_create_from_list(3,2, 8.0, 9.0, 12.0, 11.0, 15.0, 15.5);
    BlackScholesModel bl = BlackScholesModel(40, 0.04879, 0.0, Sigma, Spot);

    // Initialisation Objet MonteCarlo
    double prix = 0.0;
    double std_dev = 0.0;
    MonteCarlo mtc = MonteCarlo(&bl, &opt, rng, 0.000001, 50000);
    mtc.price(prix, std_dev);
    
    PnlVect *deltas = pnl_vect_create(G->size);
    PnlVect *std_dev_d = pnl_vect_create(G->size);

    //mtc.delta(deltas, std_dev_d); 
    // On teste
    EXPECT_NE(prix, 0.0);
    EXPECT_NE(std_dev, 0.0);
    EXPECT_DOUBLE_EQ(prix, 13.6098);
    EXPECT_DOUBLE_EQ(std_dev, 0.0300791);

    // pnl_rng_sseed(rng, 1);
    // mtc.price(prix, std_dev);
    // EXPECT_DOUBLE_EQ(prix, 1.5070789066119803);
    // EXPECT_DOUBLE_EQ(std_dev, 0.17162142732337202);



    // On libère 
    pnl_vect_free(&G);
    pnl_vect_free(&Sigma);
    pnl_vect_free(&Spot);
    pnl_vect_free(&deltas);
    pnl_vect_free(&std_dev_d);

    pnl_rng_free(&rng);
    bl.~BlackScholesModel();
    mtc.~MonteCarlo();



}


TEST(CalculPrixInstant0, OptionBasket2) {
  // Initialisation du randomisateur
    PnlRng *rng = pnl_rng_create(PNL_RNG_MERSENNE);
    pnl_rng_sseed(rng, time(NULL));

    PnlVect *G = pnl_vect_create_from_scalar(40, 0.025);

    OptionBasket opt = OptionBasket(1, 1, G->size, G, 100.0);

      // Initialisation Objet BlackScholes
    PnlVect *Sigma = pnl_vect_create_from_scalar(40, 0.2);
    PnlVect *Spot = pnl_vect_create_from_scalar(40, 100.0);
    //PnlMat *Past = pnl_mat_create_from_list(3,2, 8.0, 9.0, 12.0, 11.0, 15.0, 15.5);
    BlackScholesModel bl = BlackScholesModel(40, 0.04879, 0.7, Sigma, Spot);

    // Initialisation Objet MonteCarlo
    double prix = 0.0;
    double std_dev = 0.0;
    MonteCarlo mtc = MonteCarlo(&bl, &opt, rng, 0.000001, 50000);
    mtc.price(prix, std_dev);
    
    PnlVect *deltas = pnl_vect_create(G->size);
    PnlVect *std_dev_d = pnl_vect_create(G->size);

    //mtc.delta(deltas, std_dev_d); 
    // On teste
    EXPECT_NE(prix, 0.0);
    EXPECT_NE(std_dev, 0.0);
    EXPECT_DOUBLE_EQ(prix, 9.239495);
    EXPECT_DOUBLE_EQ(std_dev, 0.055287);

    // pnl_rng_sseed(rng, 1);
    // mtc.price(prix, std_dev);
    // EXPECT_DOUBLE_EQ(prix, 1.5070789066119803);
    // EXPECT_DOUBLE_EQ(std_dev, 0.17162142732337202);



    // On libère 
    pnl_vect_free(&G);
    pnl_vect_free(&Sigma);
    pnl_vect_free(&Spot);
    pnl_vect_free(&deltas);
    pnl_vect_free(&std_dev_d);


    pnl_rng_free(&rng);
    bl.~BlackScholesModel();
    mtc.~MonteCarlo();



}


TEST(CalculPrixInstant0, OptionBasket3) {
  // Initialisation du randomisateur
    PnlRng *rng = pnl_rng_create(PNL_RNG_MERSENNE);
    pnl_rng_sseed(rng, time(NULL));

    PnlVect *G = pnl_vect_create_from_scalar(40, 0.025);

    OptionBasket opt = OptionBasket(3, 1, G->size, G, 100.0);

      // Initialisation Objet BlackScholes
    PnlVect *Sigma = pnl_vect_create_from_scalar(40, 0.2);
    PnlVect *Spot = pnl_vect_create_from_scalar(40, 100.0);
    //PnlMat *Past = pnl_mat_create_from_list(3,2, 8.0, 9.0, 12.0, 11.0, 15.0, 15.5);
    BlackScholesModel bl = BlackScholesModel(40, 0.04879, 0.0, Sigma, Spot);

    // Initialisation Objet MonteCarlo
    double prix = 0.0;
    double std_dev = 0.0;
    MonteCarlo mtc = MonteCarlo(&bl, &opt, rng, 0.000001, 50000);
    mtc.price(prix, std_dev);
    
    PnlVect *deltas = pnl_vect_create(40);
    PnlVect *std_dev_d = pnl_vect_create(40);

    //mtc.delta(deltas, std_dev_d); 
    // On teste
    EXPECT_NE(prix, 0.0);
    EXPECT_NE(std_dev, 0.0);
    EXPECT_DOUBLE_EQ(prix, 13.631533);
    EXPECT_DOUBLE_EQ(std_dev, 0.025073);

    // pnl_rng_sseed(rng, 1);
    // mtc.price(prix, std_dev);
    // EXPECT_DOUBLE_EQ(prix, 1.5070789066119803);
    // EXPECT_DOUBLE_EQ(std_dev, 0.17162142732337202);



    // On libère 
    pnl_vect_free(&G);
    pnl_vect_free(&Sigma);
    pnl_vect_free(&Spot);
    pnl_vect_free(&deltas);
    pnl_vect_free(&std_dev_d);



    pnl_rng_free(&rng);
    bl.~BlackScholesModel();
    mtc.~MonteCarlo();


}


}