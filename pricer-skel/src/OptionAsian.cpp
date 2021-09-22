#include "OptionAsian.hpp"
#include <iostream>
using namespace std;

double OptionAsian:: payoff(const PnlMat *path){
    double payoff = 0;
    for(int i = 0; i < this->weights->size ; ++i){
        double lil_sum = 0;
        for(int j = 0; j <= this->nbTimeSteps_ ; ++j){
            lil_sum += MGET(path,i, j);
        }
        payoff += GET(this->weights, i) * lil_sum / (this->nbTimeSteps_ + 1);
    }

    return std::max(payoff - this->strike, 0.0);
}