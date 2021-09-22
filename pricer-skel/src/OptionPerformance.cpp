#include "OptionPerformance.hpp"
#include <iostream>
using namespace std;

double OptionPerformance :: payoff(const PnlMat *path){
    double payoff = 1;
    
    for(int i = 1; i <= this->nbTimeSteps_; ++i){
        double num = 0.0;
        double denom = 0.0;
        double val = 0.0;
        for(int j = 0; j < this->weights->size ; ++j){
            num += GET(this->weights, j) * MGET(path, i, j);
            denom += GET(this->weights, j) * MGET(path, i-1, j);
        }
        val = max((num / denom) - 1, 0.0);
        payoff += val;
    }

    return payoff;
}