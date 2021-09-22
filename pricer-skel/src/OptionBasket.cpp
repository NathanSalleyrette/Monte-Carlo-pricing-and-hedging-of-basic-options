#include "OptionBasket.hpp"
#include <iostream>
using namespace std;

double OptionBasket :: payoff(const PnlMat *path){

    double payoff = 0;

    for(int i = 0; i < this->weights->size ; ++i){    
        payoff += GET(this->weights, i) * MGET(path, path->m - 1, i);
    }

    return  std::max(payoff - this->strike, 0.0);
}