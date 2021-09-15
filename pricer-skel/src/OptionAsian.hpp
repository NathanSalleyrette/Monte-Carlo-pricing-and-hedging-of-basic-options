#pragma once
#include <string>
#include "pnl/pnl_vector.h"
#include "pnl/pnl_matrix.h"
#include "Option.hpp"
class OptionAsian : public Option
{
private:
    const PnlVect *const weights;
    
    //string[] underlyingShares;
public:

    OptionAsian(double T_, int nbTimeSteps_, int size_, const PnlVect *const weights, double strike)
        :Option(T_, nbTimeSteps_, size_, strike),
        weights(weights)
    { }

    ~OptionAsian() { }

    OptionAsian(const OptionAsian &OptionAsian)
        :Option(OptionAsian.T_, OptionAsian.nbTimeSteps_, OptionAsian.size_, OptionAsian.strike),
        weights(OptionAsian.weights)
    { }

    double payoff(const PnlMat *path) override;
};