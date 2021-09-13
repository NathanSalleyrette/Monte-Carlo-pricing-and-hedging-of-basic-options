#include "OptionBasket.hpp"

void OptionBasket :: payoff(const PnlMat *path){
    double payoff = 0;
    
    for(int i = 0; i < path.m; i++){
        for(int j = 0; j < path.n; j++){
            payoff += this->weights[j] * path[j,i];
        }
    }

    return max(payoff - this.strike, 0);

}