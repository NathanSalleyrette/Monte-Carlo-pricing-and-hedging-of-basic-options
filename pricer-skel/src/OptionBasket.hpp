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

    OptionBasket(double T_, int nbTimeSteps_, int size_, const PnlVect *const weights)
        :Option(T_, nbTimeSteps_, size_),
        weights(weights)
        //,underlyingShares(underlyingShares)
    { }

    ~OptionBasket() { }

    OptionBasket(const OptionBasket &OptionBasket)
        :Option(OptionBasket.T_, OptionBasket.nbTimeSteps_, OptionBasket.size_),
        weights(OptionBasket.weights)
        //,underlyingShares(OptionBasket.underlyingShares)
    { }

    double payoff(const PnlMat *path) {
        return 0.0;
    }
};


