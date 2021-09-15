#pragma once
#include <string>
#include "pnl/pnl_vector.h"
#include "pnl/pnl_matrix.h"
#include "Option.hpp"
class OptionPerformance : public Option
{
private:
    const PnlVect *const weights;
    
    //string[] underlyingShares;
public:

    OptionPerformance(double T_, int nbTimeSteps_, int size_, const PnlVect *const weights, double strike)
        :Option(T_, nbTimeSteps_, size_, strike),
        weights(weights)
    { }

    ~OptionPerformance() { }

    OptionPerformance(const OptionPerformance &OptionPerformance)
        :Option(OptionPerformance.T_, OptionPerformance.nbTimeSteps_, OptionPerformance.size_, OptionPerformance.strike),
        weights(OptionPerformance.weights)
    { }

    double payoff(const PnlMat *path) override;
};