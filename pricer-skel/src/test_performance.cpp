#include "pnl/pnl_vector.h"
#include "OptionAsian.cpp"
#include "OptionPerformance.cpp"
#include "MonteCarlo.cpp"
#include "BlackScholesModel.cpp"
#include <gtest/gtest.h>

namespace {


TEST(CalculPrixInstant0, OptionPerformance) {
  // Initialisation du randomisateur
    PnlRng *rng = pnl_rng_create(PNL_RNG_MERSENNE);
    pnl_rng_sseed(rng, time(NULL));

    PnlVect *G = pnl_vect_create_from_scalar(5, 0.2);

    OptionPerformance opt = OptionPerformance(2.0, 12, G->size, G, 100.0);

      // Initialisation Objet BlackScholes
    PnlVect *Sigma = pnl_vect_create_from_scalar(5, 0.2);
    PnlVect *Spot = pnl_vect_create_from_scalar(5, 100.0);
    //PnlMat *Past = pnl_mat_create_from_list(3,2, 8.0, 9.0, 12.0, 11.0, 15.0, 15.5);
    BlackScholesModel bl = BlackScholesModel(40, 0.03, 0.5, Sigma, Spot);

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
    EXPECT_NEAR(prix, 1.2573, 0.000589736);
    // EXPECT_DOUBLE_EQ(prix, 1.2573);
    // EXPECT_DOUBLE_EQ(std_dev, 0.000589736);

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
