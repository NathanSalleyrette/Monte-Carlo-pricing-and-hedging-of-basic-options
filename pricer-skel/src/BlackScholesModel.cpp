#include "BlackScholesModel.hpp"

void BlackScholesModel :: asset(PnlMat *path, double T, int nbTimeSteps, PnlRng *rng){

    PnlMat *L = pnl_mat_create(path->m, path->m);

    for(int i = 0; i < path->m* path->m ; i++){
        L->array[i] = this->rho_;
    }

    for(int i = 0; i < path->m; i++){
        L->array[i*(path->m + 1)] = 1;
        path->array[i*path->n] = this->spot_->array[i];
    }
    pnl_mat_chol(L);
    double step = T/(double)nbTimeSteps;

    for(int i = 0; i < path->n-1 ; i++ ){
        PnlVect *G = pnl_vect_create(path->m);
        pnl_vect_rng_normal(G, path->m, rng);
        for(int s = 0; s < path->m; s++){
            double LG = 0;
            for(int j = 0 ; j < path->m ; j++){
                LG += L->array[s*path->m +j] * G->array[j];
            } 
            path->array[s*path->n + i+1] = path->array[s*path->n + i] * exp((this->r_ - this->sigma_->array[s] * this->sigma_->array[s] )/2.0 * step + this->sigma_->array[s]*sqrt(step)*LG  );
        }
        //On libère G et L a la sortie des for
        pnl_vect_free(&G);
    }
    pnl_mat_free(&L);

}

void BlackScholesModel :: shiftAsset(PnlMat *shift_path, const PnlMat *path, int d, double h, double t, double timestep){
    
    for(int i = 0; i < path->n ; i++){
        if(i*timestep <= t){
            shift_path->array[d*path->n + i] = path->array[d*path->n + i];
        }else {
            shift_path->array[d*path->n + i] = (1.0+h)*path->array[d*path->n + i];
        }
    }

}