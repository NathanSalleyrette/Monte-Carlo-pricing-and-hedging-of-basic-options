#include "BlackScholesModel.hpp"

void BlackScholesModel :: asset(PnlMat *path, double T, int nbTimeSteps, PnlRng *rng){
    
    double step = T/(double)nbTimeSteps;

    for(int s = 0; s < path->m ; s++ ){
        for(int i = 0; i < path->n; i ++){
            double t = i*step;
            double bt = sqrt(t) * pnl_rng_normal(rng);
            path->array[i, s] = this->spot_[s] * exp((this->r_ * - this->sigma_[s])/2.0 * t + this->sigma_[s]*bt);
            cout << path->array[i, s] << endl; 
        }
    }

}