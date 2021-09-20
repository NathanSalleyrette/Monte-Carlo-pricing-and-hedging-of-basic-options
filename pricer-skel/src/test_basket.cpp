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
    MonteCarlo mtc = MonteCarlo(&bl, &opt, rng, 0.1, 50000);
    mtc.price(prix, std_dev);
    
    PnlVect *deltas = pnl_vect_create(G->size);
    PnlVect *std_dev_d = pnl_vect_create(G->size);

    mtc.delta(deltas, std_dev_d); 
    // On teste
    EXPECT_NE(prix, 0.0);
    EXPECT_NE(std_dev, 0.0);
    // EXPECT_NEAR(prix, 13.603843, 1.96 * 0.025111); (vient du .data)
    EXPECT_NEAR(prix, 13.6098, 1.96 * 0.0252271);
    // EXPECT_DOUBLE_EQ(prix, 13.6098);
    // EXPECT_DOUBLE_EQ(std_dev, 0.0300791);

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
    MonteCarlo mtc = MonteCarlo(&bl, &opt, rng, 0.1, 50000);
    mtc.price(prix, std_dev);
    
    PnlVect *deltas = pnl_vect_create(G->size);
    PnlVect *std_dev_d = pnl_vect_create(G->size);

    mtc.delta(deltas, std_dev_d); 
    // On teste
    EXPECT_NE(prix, 0.0);
    EXPECT_NE(std_dev, 0.0);
    // EXPECT_NEAR(prix, 9.239495, 1.96 * 0.055287); (vient du .data)
    EXPECT_NEAR(prix, 9.21359, 1.96 * 0.0555551);
    // EXPECT_DOUBLE_EQ(prix, 9.239495);
    // EXPECT_DOUBLE_EQ(std_dev, 0.055287);

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


TEST(CalculPrixInstantnonull, OptionBasket1) {
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
    double prix0 = 0.0;
    double std_dev0 = 0.0;
    MonteCarlo mtc = MonteCarlo(&bl, &opt, rng, 0.1, 50000);
    mtc.price(prix0, std_dev0);
    
    PnlVect *deltas = pnl_vect_create(G->size);
    PnlVect *std_dev_d = pnl_vect_create(G->size);

    mtc.delta(deltas, std_dev_d);
    // decalage vers la gauche
    OptionBasket opt2 = OptionBasket(1, 1, G->size, G, 100.0);
    PnlMat *past = pnl_mat_create_from_scalar(G->size, 2, 100.0);
    double prix = 0.0;
    double std_dev = 0.0;

    mtc.price(past, 2.0, prix, std_dev);

    double prix2 = 0.0;
    double std_dev2 = 0.0;
    MonteCarlo mtc2 = MonteCarlo(&bl, &opt2, rng, 0.1, 50000);

    mtc2.price(prix2, std_dev2);



    // On teste
    EXPECT_NE(prix, prix2);
    EXPECT_NE(std_dev, std_dev2);
    EXPECT_NEAR(prix, prix2, 1.96 * std_dev2);
    // EXPECT_NEAR(prix, 13.6098, 0.0300791);


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

TEST(CalculPrixInstantnonull, OptionBasketNouvelle) {
  // Initialisation du randomisateur
    PnlRng *rng = pnl_rng_create(PNL_RNG_MERSENNE);
    pnl_rng_sseed(rng, time(NULL));

    PnlVect *G = pnl_vect_create_from_scalar(40, 0.025);

    OptionBasket opt = OptionBasket(4, 4, G->size, G, 100.0);

      // Initialisation Objet BlackScholes
    PnlVect *Sigma = pnl_vect_create_from_scalar(40, 0.2);
    PnlVect *Spot = pnl_vect_create_from_scalar(40, 100.0);
    //PnlMat *Past = pnl_mat_create_from_list(3,2, 8.0, 9.0, 12.0, 11.0, 15.0, 15.5);
    BlackScholesModel bl = BlackScholesModel(40, 0.04879, 0.0, Sigma, Spot);

    // Initialisation Objet MonteCarlo
    double prix0 = 0.0;
    double std_dev0 = 0.0;
    MonteCarlo mtc = MonteCarlo(&bl, &opt, rng, 0.1, 50000);
    mtc.price(prix0, std_dev0);
    
    PnlVect *deltas = pnl_vect_create(G->size);
    PnlVect *std_dev_d = pnl_vect_create(G->size);

    mtc.delta(deltas, std_dev_d);
    // decalage vers la gauche

    double prix = 0.0;
    double std_dev = 0.0;
    PnlMat *past = pnl_mat_create_from_scalar(G->size, 2, 100.0);

    mtc.price(past, 0.5, prix, std_dev);

    OptionBasket opt2 = OptionBasket(3.5, 3, G->size, G, 100.0);
    double prix2 = 0.0;
    double std_dev2 = 0.0;
    MonteCarlo mtc2 = MonteCarlo(&bl, &opt2, rng, 0.1, 50000);

    mtc2.price(prix2, std_dev2);



    // On teste
    EXPECT_NE(prix, prix2);
    EXPECT_NE(std_dev, std_dev2);
    EXPECT_NEAR(prix, prix2, 1.96 * std_dev2);
    // EXPECT_NEAR(prix, 13.6098, 0.0300791);


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

TEST(CalculPrixInstanttdicretise, OptionBasketNouvelle) {
  // Initialisation du randomisateur
    PnlRng *rng = pnl_rng_create(PNL_RNG_MERSENNE);
    pnl_rng_sseed(rng, time(NULL));

    PnlVect *G = pnl_vect_create_from_scalar(40, 0.025);

    OptionBasket opt = OptionBasket(4, 4, G->size, G, 100.0);

      // Initialisation Objet BlackScholes
    PnlVect *Sigma = pnl_vect_create_from_scalar(40, 0.2);
    PnlVect *Spot = pnl_vect_create_from_scalar(40, 100.0);
    //PnlMat *Past = pnl_mat_create_from_list(3,2, 8.0, 9.0, 12.0, 11.0, 15.0, 15.5);
    BlackScholesModel bl = BlackScholesModel(40, 0.04879, 0.0, Sigma, Spot);

    // Initialisation Objet MonteCarlo
    double prix0 = 0.0;
    double std_dev0 = 0.0;
    MonteCarlo mtc = MonteCarlo(&bl, &opt, rng, 0.1, 50000);
    mtc.price(prix0, std_dev0);
    
    PnlVect *deltas = pnl_vect_create(G->size);
    PnlVect *std_dev_d = pnl_vect_create(G->size);

    mtc.delta(deltas, std_dev_d);
    // decalage vers la gauche

    double prix = 0.0;
    double std_dev = 0.0;
    PnlMat *past = pnl_mat_create_from_scalar(G->size, 2, 100.0);

    mtc.price(past, 1, prix, std_dev);

    OptionBasket opt2 = OptionBasket(3, 3, G->size, G, 100.0);
    double prix2 = 0.0;
    double std_dev2 = 0.0;
    MonteCarlo mtc2 = MonteCarlo(&bl, &opt2, rng, 0.1, 50000);

    mtc2.price(prix2, std_dev2);



    // On teste
    EXPECT_NE(prix, prix2);
    EXPECT_NE(std_dev, std_dev2);
    EXPECT_NEAR(prix, prix2, 1.96 * std_dev2);
    // EXPECT_NEAR(prix, 13.6098, 0.0300791);


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