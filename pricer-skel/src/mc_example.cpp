#include <iostream>
#include <ctime>
#include "pnl/pnl_random.h"
#include "pnl/pnl_vector.h"
#include "OptionBasket.cpp"
#include "MonteCarlo.cpp"
#include "BlackScholesModel.cpp"
#include "OptionAsian.cpp"
#include "OptionPerformance.cpp"

using namespace std;

int main()
{
    PnlVect *G0 = pnl_vect_create(3);
    G0->array[0] = 0.5;
    G0->array[1] = 0.4;
    G0->array[2] = 0.1;
    OptionPerformance opt = OptionPerformance(1.0, 1, 3, G0, 10.0);
    PnlMat *M0 = pnl_mat_create(3,2) ;
    M0->array[0] = 6.0;
    M0->array[1] = 3.0;
    M0->array[2] = 2.0;
    M0->array[3] = 4.0;
    M0->array[4] = 5.0;
    M0->array[5] = 1.0;
    cout << opt.payoff(M0) << endl;
    PnlVect *G = pnl_vect_new();
    PnlRng *rng = pnl_rng_create(PNL_RNG_MERSENNE);
    int M = 1E5;
    int dim = 2;
    pnl_rng_sseed(rng, 0);

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

    

    double acc = 0., var = 0;

    for (int i = 0; i < M; i++)
    {
        pnl_vect_rng_normal(G, dim, rng);
        double tmp = pnl_vect_norm_two(G);
        acc += tmp;
        var += tmp * tmp;
    }

    acc /= M;
    var = var / M - acc * acc;

    cout << "E[||G||_2] = " << acc << endl;
    cout << "Var(||G||_2) = " << var << endl;

    pnl_vect_free(&G);
    pnl_rng_free(&rng);
    return 0;
}
