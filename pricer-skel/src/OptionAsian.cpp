#include "OptionAsian.hpp"
#include <iostream>
using namespace std;

double OptionAsian:: payoff(const PnlMat *path){
    double payoff = 0;
    
    for(int i = 0; i < this->weights->size ; ++i){
        //cout << this->weights->array[i] << " " << path->array[i, this->nbTimeSteps_ ] <<  endl;
        // payoff = somme (poids de l'actif * somme(strikes aux N+1 dates) / N+1)
        double lil_sum = 0;
        for(int j = 0; j <= this->nbTimeSteps_ ; ++j){
            lil_sum += path->array[i*path->n + j];
        }
        payoff += this->weights->array[i] * lil_sum / (this->nbTimeSteps_ + 1);
    }

    return  std::max(payoff - this->strike, 0.0);

}