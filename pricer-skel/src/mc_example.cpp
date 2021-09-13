#include <iostream>
#include <ctime>
#include "pnl/pnl_random.h"
#include "pnl/pnl_vector.h"
#include "OptionBasket.cpp"

using namespace std;

int main()
{
    PnlVect *G0 = pnl_vect_create(3);
    G0->array[0] = 0.5;
    G0->array[1] = 0.4;
    G0->array[2] = 0.1;
    OptionBasket opt = OptionBasket(1.0, 1, 1, G0, 2.0);
    PnlMat *M0 = pnl_mat_create(1,3) ;
    M0->array[0] = 6.0;
    M0->array[1] = 3.0;
    M0->array[2] = 2.0;
    cout << opt.payoff(M0) << endl;
    PnlVect *G = pnl_vect_new();
    PnlRng *rng = pnl_rng_create(PNL_RNG_MERSENNE);
    int M = 1E5;
    int dim = 2;
    pnl_rng_sseed(rng, time(NULL));

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
