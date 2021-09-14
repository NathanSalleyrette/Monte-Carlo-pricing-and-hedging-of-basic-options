#include "OptionBasket.hpp"
using namespace std;

double OptionBasket :: payoff(const PnlMat *path){
    double payoff = 0;
    
    for(int i = 0; i < this->weights->size ; ++i){
        cout << this->weights->array[i] << " " << path->array[i, this->nbTimeSteps_ ] <<  endl;
        payoff += this->weights->array[i] * path->array[i, this->nbTimeSteps_ ];
        cout << payoff << endl;
    }

    return  std::max(payoff - this->strike, 0.0);

}