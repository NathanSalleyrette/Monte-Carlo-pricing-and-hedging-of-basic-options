#include "pnl/pnl_random.h"
#include "pnl/pnl_vector.h"
#include "OptionBasket.cpp"
#include "OptionAsian.cpp"
#include "MonteCarlo.cpp"
#include "BlackScholesModel.cpp"
#include <gtest/gtest.h>

namespace {



TEST(PortfolioValue, basket_2) {
    // Initialisation du randomisateur
    PnlRng *rng = pnl_rng_create(PNL_RNG_MERSENNE);
    pnl_rng_sseed(rng, time(NULL));

    PnlVect *G = pnl_vect_create_from_scalar(2, 0.5);

    OptionBasket opt = OptionBasket(1, 1, G->size, G, 100.0);

      // Initialisation Objet BlackScholes
    PnlVect *Sigma = pnl_vect_create_from_scalar(G->size, 0.2);
    PnlVect *Spot = pnl_vect_create_from_scalar(G->size, 100.0);
    //PnlMat *Past = pnl_mat_create_from_list(3,2, 8.0, 9.0, 12.0, 11.0, 15.0, 15.5);
    BlackScholesModel bl = BlackScholesModel(G->size, 0.04879, 0.0, Sigma, Spot);
    PnlMat *RealPath = pnl_mat_create_from_file("../data-hedge/basket_2d_market.dat");
    //PnlMat *RealPath = pnl_mat_transpose(Path);

    // Initialisation Objet MonteCarlo
    double prix = 0.0;
    double std_dev = 0.0;
    double error = 0;
    MonteCarlo mtc = MonteCarlo(&bl, &opt, rng, 0.1, 50000);
    mtc.price(prix, std_dev);
    PnlVect *deltas = pnl_vect_create(G->size);
    PnlVect *std_dev_d = pnl_vect_create(G->size);

    mtc.delta(deltas, std_dev_d); 

    mtc.profitAndLoss(RealPath, error, prix, std_dev);
    cout << "error    " << error << endl;
    cout << "prix     " << prix << endl;
    cout << "std_dev  " << std_dev << endl;


    double fin = 0.0;
}


TEST(PortfolioValue, basket) {
    // Initialisation du randomisateur
    PnlRng *rng = pnl_rng_create(PNL_RNG_MERSENNE);
    pnl_rng_sseed(rng, time(NULL));

    PnlVect *G = pnl_vect_create_from_scalar(5, 0.2);

    OptionBasket opt = OptionBasket(1, 1, G->size, G, 100.0);

      // Initialisation Objet BlackScholes
    PnlVect *Sigma = pnl_vect_create_from_scalar(G->size, 0.2);
    PnlVect *Spot = pnl_vect_create_from_scalar(G->size, 100.0);
    //PnlMat *Past = pnl_mat_create_from_list(3,2, 8.0, 9.0, 12.0, 11.0, 15.0, 15.5);
    BlackScholesModel bl = BlackScholesModel(G->size, 0.04879, 0.0, Sigma, Spot);
    PnlMat *RealPath = pnl_mat_create_from_file("../data-hedge/basket_market.dat");
    //PnlMat *RealPath = pnl_mat_transpose(Path);

    // Initialisation Objet MonteCarlo
    double prix = 0.0;
    double std_dev = 0.0;
    double error = 0;
    MonteCarlo mtc = MonteCarlo(&bl, &opt, rng, 0.1, 50000);
    mtc.price(prix, std_dev);
    PnlVect *deltas = pnl_vect_create(G->size);
    PnlVect *std_dev_d = pnl_vect_create(G->size);

    mtc.delta(deltas, std_dev_d); 

    mtc.profitAndLoss(RealPath, error, prix, std_dev);
    cout << "error    " << error << endl;
    cout << "prix     " << prix << endl;
    cout << "std_dev  " << std_dev << endl;
    double fin = 0.0;
}

TEST(PortfolioValue, call_market) {
    // Initialisation du randomisateur
    PnlRng *rng = pnl_rng_create(PNL_RNG_MERSENNE);
    pnl_rng_sseed(rng, time(NULL));

    PnlVect *G = pnl_vect_create_from_scalar(1, 1);

    OptionBasket opt = OptionBasket(1, 1, G->size, G, 100.0);

      // Initialisation Objet BlackScholes
    PnlVect *Sigma = pnl_vect_create_from_scalar(G->size, 0.2);
    PnlVect *Spot = pnl_vect_create_from_scalar(G->size, 100.0);
    //PnlMat *Past = pnl_mat_create_from_list(3,2, 8.0, 9.0, 12.0, 11.0, 15.0, 15.5);
    BlackScholesModel bl = BlackScholesModel(G->size, 0.04879, 0.0, Sigma, Spot);
    PnlMat *RealPath = pnl_mat_create_from_file("../data-hedge/call_market.dat");
    //PnlMat *RealPath = pnl_mat_transpose(Path);

    // Initialisation Objet MonteCarlo
    double prix = 0.0;
    double std_dev = 0.0;
    double error = 0;
    MonteCarlo mtc = MonteCarlo(&bl, &opt, rng, 0.1, 50000);
    mtc.price(prix, std_dev);
    PnlVect *deltas = pnl_vect_create(G->size);
    PnlVect *std_dev_d = pnl_vect_create(G->size);

    mtc.delta(deltas, std_dev_d); 

    mtc.profitAndLoss(RealPath, error, prix, std_dev);
    cout << "error    " << error << endl;
    cout << "prix     " << prix << endl;
    cout << "std_dev  " << std_dev << endl;
    double fin = 0.0;
}

}