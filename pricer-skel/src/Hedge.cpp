#include <iostream>
#include <ctime>
#include "HedgingResults.cpp"
#include "OptionBasket.cpp"
#include "OptionPerformance.cpp"
#include "OptionAsian.cpp"
#include "MonteCarlo.cpp"
#include "BlackScholesModel.cpp"
#include "3rdparty/jlparser/parser.cpp"
using namespace std;

int main(int argc, char **argv)
{
    double T, r, c, strike;
    PnlVect *spot, *sigma, *divid, *pCoef;
    string type;
    int size, timeStep, nbDates;
    size_t n_samples;
    PnlRng *rng = pnl_rng_create(PNL_RNG_MERSENNE);
    pnl_rng_sseed(rng, time(NULL));
    char *infile = argv[2];
    Param *P = new Parser(infile);
    P->extract("option type", type);
    P->extract("maturity", T);
    P->extract("option size", size);
    P->extract("timeStep number", timeStep);
    P->extract("spot", spot, size);
    P->extract("volatility", sigma, size);
    P->extract("interest rate", r);
    P->extract("correlation", c);
    P->extract("payoff coefficients", pCoef, size);
    P->extract("hedging dates number", nbDates);
    char *marketInfile = argv[1];
    
    PnlMat *realPath = pnl_mat_create_from_file(marketInfile);
    //PnlMat *transposeRealPath = pnl_mat_transpose(realpath)
    // for (int i = 0; i < transposeRealPath->nm){
    //     transposeRealPath->array[i] = 
    // }
    
    if (P->extract("dividend rate", divid, size, true) == false)
    {
        divid = pnl_vect_create_from_zero(size);
    }

    P->extract("sample number", n_samples);
    // PnlVect *delta = pnl_vect_create(size);
    // PnlVect *delta_std_dev = pnl_vect_create(size);
    double prix = 0.0;
    double prix_std_dev = 0.0;
    double error = 0.0;
    if (type == "basket"){
        P->extract("strike", strike);
        OptionBasket opt = OptionBasket(T, timeStep, size, pCoef, strike);
        BlackScholesModel bl = BlackScholesModel(size, r, c, sigma, spot);
        MonteCarlo mtc = MonteCarlo(&bl, &opt, rng, 0.01, n_samples);
        mtc.profitAndLoss(realPath, error, prix, prix_std_dev);
        // mtc.price(prix, prix_std_dev); 
        // mtc.delta(delta,delta_std_dev);       
    }
    else if (type == "asian"){
        P->extract("strike", strike);
        OptionAsian opt = OptionAsian(T, timeStep, size, pCoef, strike);
        BlackScholesModel bl = BlackScholesModel(size, r, c, sigma, spot);
        MonteCarlo mtc = MonteCarlo(&bl, &opt, rng, 0.01, n_samples);
        mtc.profitAndLoss(realPath, error, prix, prix_std_dev);
        // mtc.price(prix, prix_std_dev); 
        // mtc.delta(delta,delta_std_dev); 
    }
    else{
        OptionPerformance opt = OptionPerformance(T, timeStep, size, pCoef, spot->array[0]);
        BlackScholesModel bl = BlackScholesModel(size, r, c, sigma, spot);
        MonteCarlo mtc = MonteCarlo(&bl, &opt, rng, 0.01, n_samples);
        mtc.profitAndLoss(realPath, error, prix, prix_std_dev);
        // mtc.price(prix, prix_std_dev); 
        // mtc.delta(delta,delta_std_dev); 
    }
    HedgingResults res = HedgingResults(prix, prix_std_dev, error);
    std::cout << res << std::endl;
    
    pnl_vect_free(&spot);
    pnl_vect_free(&sigma);
    pnl_vect_free(&divid);
    delete P;
    exit(0);
}
