#pragma once
#include <string>
#include "pnl/pnl_vector.h"
#include "pnl/pnl_matrix.h"
#include "Option.hpp"
class OptionBasket : public Option
{
private:
    const PnlVect *const weights;
    
    //string[] underlyingShares;
public:

    OptionBasket(double T_, int nbTimeSteps_, int size_, const PnlVect *const weights, double strike)
        :Option(T_, nbTimeSteps_, size_, strike),
        weights(weights)
    { }

    ~OptionBasket() { }

    OptionBasket(const OptionBasket &OptionBasket)
        :Option(OptionBasket.T_, OptionBasket.nbTimeSteps_, OptionBasket.size_, OptionBasket.strike),
        weights(OptionBasket.weights)
    { }

    double payoff(const PnlMat *path) override;
};


