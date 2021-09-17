
#include "pnl/pnl_random.h"
#include "pnl/pnl_vector.h"
#include "OptionBasket.cpp"
#include "OptionAsian.cpp"
#include "MonteCarlo.cpp"
#include "BlackScholesModel.cpp"
#include <gtest/gtest.h>

namespace {


TEST(OptionAsian1, CalculAuTemps0) {
  // Initialisation du randomisateur
    PnlRng *rng = pnl_rng_create(PNL_RNG_MERSENNE);
    pnl_rng_sseed(rng, time(NULL));

    PnlVect *G = pnl_vect_create(2);
    G->array[0] = 0.5;
    G->array[1] = 0.5;
    OptionAsian opt = OptionAsian(1.5, 150, G->size, G, 100);

      // Initialisation Objet BlackScholes
    PnlVect *Sigma = pnl_vect_create_from_list(2, 0.2, 0.2);
    PnlVect *Spot = pnl_vect_create_from_list(2, 100.0, 100.0);
    BlackScholesModel bl = BlackScholesModel(2, 0.02, 0.0, Sigma, Spot);

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
    EXPECT_NEAR(prix, 4.7208, 0.0300791);
    // EXPECT_DOUBLE_EQ(prix, 4.7208);
    // EXPECT_DOUBLE_EQ(std_dev, 0.0300791);

    pnl_rng_free(&rng);


}
}