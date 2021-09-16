#include "OptionBasket.hpp"
#include <iostream>
using namespace std;

double OptionBasket :: payoff(const PnlMat *path){
    double payoff = 0;
    
    for(int i = 0; i < this->weights->size ; ++i){
        //cout << this->weights->array[i] << " " << path->array[i, this->nbTimeSteps_ ] <<  endl;
        // payoff = somme (poids de l'actif * la valeur de l'actif à la maturité) 
       
        payoff += this->weights->array[i] * path->array[(i+1)*path->n - 1];
        //payoff += GET(this->weights,i) * MGET(path, i, path->n - 1);

        //cout << payoff << endl;
    }

    return  std::max(payoff - this->strike, 0.0);

}