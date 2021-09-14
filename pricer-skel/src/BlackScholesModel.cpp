#include "BlackScholesModel.hpp"

void BlackScholesModel :: asset(PnlMat *path, double T, int nbTimeSteps, PnlRng *rng){

    PnlMat *L = pnl_mat_create(path->m, path->m);

    for(int i = 0; i < path->m* path->m ; i++){
        L->array[i] = this->rho_;
    }

    for(int i = 0; i < path->m; i++){
        L->array[i*(path->m + 1)] = 1;
        path->array[i] = this->spot_->array[i];
    }
    pnl_mat_chol(L);
    double step = T/(double)nbTimeSteps;

    for(int i = 0; i < path->n ; i++ ){
        PnlVect *G = pnl_vect_create(path->m);
        pnl_vect_rng_normal(G, path->m, rng);
        for(int s = 0; s < path->m; s++){
            double t = i*step;
            double LG = 0;
            for(int j = 0 ; j < path->m ; j++){
                LG += L->array[s*path->m +j] * G->array[j];
            } 
            path->array[(i+1)*path->m + s] = path->array[i*path->m + s] * exp((this->r_ * - this->sigma_->array[s] * this->sigma_->array[s] )/2.0 * step + this->sigma_->array[s]*sqrt(step)*LG  );
            cout << path->array[i*path->m + s] << endl; 
        }
    }

}