#include "OptionPerformance.hpp"
#include <iostream>
using namespace std;

double OptionPerformance :: payoff(const PnlMat *path){
    double payoff = 1;
    
    for(int i = 1; i <= this->nbTimeSteps_; ++i){
        //cout << this->weights->array[i] << " " << path->array[i, this->nbTimeSteps_ ] <<  endl;
        // payoff = somme (poids de l'actif * la valeur de l'actif à la maturité)
        double num = 0.0;
        double denom = 0.0;
        double val = 0.0;
        for(int j = 0; j < this->weights->size ; ++j){
            
            num += GET(this->weights, j) * MGET(path, j, i);
            denom += GET(this->weights, j) * MGET(path, j, i - 1);
            
            //num += this->weights->array[j] * path->array[j*path->n +i];
            //denom += this->weights->array[j] * path->array[j*path->n +i-1];
        }
        val = max((num / denom) - 1, 0.0);
        payoff += val;
        //cout << payoff << endl;
    }

    return payoff;

}