#include "BlackScholesModel.hpp"

void BlackScholesModel :: asset(PnlMat *path, double T, int nbTimeSteps, PnlRng *rng){

    PnlMat *L = pnl_mat_create_from_scalar(path->m, path->m, this->rho_);
    PnlVect *LignL = pnl_vect_create(path->m);
    for(int i = 0; i < path->m; i++){
        pnl_mat_set_diag(L, 1.0, i);
    }

    pnl_mat_set_col(path, this->spot_, 0);
    pnl_mat_chol(L);
    double step = T/(double)nbTimeSteps;

    PnlMat *Gmat = pnl_mat_create(nbTimeSteps, path->m);
    pnl_mat_rng_normal(Gmat, nbTimeSteps, path->m, rng);

    PnlVect *LignG = pnl_vect_create(path->m);

    for(int j = 0; j < path->m; ++j){
        pnl_mat_get_row(LignL, L, j);
        for(int i = 0; i < path->n - 1; ++i){
            double LG = 0.0;
            pnl_mat_get_row(LignG, Gmat, i);
            LG = pnl_vect_scalar_prod(LignL, LignG);
            
            double Sj = GET(this->sigma_, j);
            MLET(path, j,i+1) = MGET(path,j , i) * exp((this->r_ - Sj * Sj /2.0) * step + Sj*sqrt(step)*LG  );


        }
    }

    // for(int i = 0; i < path->n-1 ; i++ ){
    //     PnlVect *G = pnl_vect_create(path->m);
    //     pnl_vect_rng_normal(G, path->m, rng);
    //     for(int s = 0; s < path->m; s++){
    //         double LG = 0;
    //         for(int j = 0 ; j < path->m ; j++){
    //             LG += L->array[s*path->m +j] * G->array[j];
    //         } 
    //         path->array[s*path->n + i+1] = path->array[s*path->n + i] * exp((this->r_ - this->sigma_->array[s] * this->sigma_->array[s] )/2.0 * step + this->sigma_->array[s]*sqrt(step)*LG  );
    //     }
    //     //On libère G et L a la sortie des for
    //     pnl_vect_free(&G);
    // }
    // pnl_mat_free(&L);

}

void BlackScholesModel :: asset(PnlMat *path, double t, double T, int nbTimeSteps, PnlRng *rng, const PnlMat *past) {
    int simuremains = path->n - past->n + 1;
    PnlMat *L = pnl_mat_create_from_scalar(path->m, path->m, this->rho_);
    
    PnlVect *LignL = pnl_vect_create(path->m);
    
    for(int i = 0; i < path->m; i++){
        pnl_mat_set_diag(L, 1.0, i);
    }
    // Pas besoin de copier past dans path, c'est déjà fait
    pnl_mat_chol(L);
    double step = T/(double)nbTimeSteps;

    // Création de la matrice G randomisée
    PnlMat *GMat = pnl_mat_create(simuremains, path->m);
    pnl_mat_rng_normal(GMat, simuremains, path->m, rng);
    PnlVect *GLign = pnl_vect_create(path->m);

    // Pour chaque actif, on fait :
    for (int j = 0; j < path->m; j++) {
        
        PnlVect *S = pnl_vect_create(simuremains);
        double previous = 1;
        for(int s = 0; s < S->size; s++) {

            double LG = 0;
            pnl_mat_get_row(GLign, GMat, s);
            pnl_mat_get_row(LignL, L, j);
            LG = pnl_vect_scalar_prod(GLign, LignL);

            // Attention, il ne faut pas créer de nouveaux G pour chaque actif (à voir avec le prof)
            // PnlVect *G = pnl_vect_create(path->m);
            // pnl_vect_rng_normal(G, path->m, rng);
            // double LG = 0;

            // pnl_mat_get_row(LignL, L, j);
            // LG = pnl_vect_scalar_prod(LignL, G);
            
            if(s != 0) previous = GET(S, s-1);
            double Sj = GET(this->sigma_, j);
            LET(S, s) = previous * exp((this->r_ - Sj* Sj /2.0) * step + Sj*sqrt(step)*LG  );

        }
        pnl_vect_mult_scalar(S, past->array[(j+1)*past->n - 1]);

        for (int s = 0; s < S->size; s++) {
            path->array[(j+1)*path->n - simuremains + s] = S->array[s];
        }
    }



    // for(int i = 0; i < path->n-1 ; i++ ){
    //     PnlVect *G = pnl_vect_create(path->m);
    //     pnl_vect_rng_normal(G, path->m, rng);
    //     for(int s = 0; s < path->m; s++){
    //         double LG = 0;
    //         for(int j = 0 ; j < path->m ; j++){
    //             LG += L->array[s*path->m +j] * G->array[j];
    //         } 
    //     }
    //     //On libère G et L a la sortie des for
    //     pnl_vect_free(&G);
    // }
    // pnl_mat_free(&L);
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