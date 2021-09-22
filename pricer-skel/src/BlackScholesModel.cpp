#include "BlackScholesModel.hpp"

void BlackScholesModel :: asset(PnlMat *path, double T, int nbTimeSteps, PnlRng *rng){

    PnlMat *L = this->El;
    PnlVect *LignL = this->LignEl;
 
    pnl_mat_set_row(path, this->spot_, 0);

    double step = T/(double)nbTimeSteps;
    double sqrtstep = sqrt(step);
    double freeRate = this->r_;
    PnlMat *Gmat = pnl_mat_create(nbTimeSteps, path->n);
    pnl_mat_rng_normal(Gmat, nbTimeSteps, path->n, rng);

    PnlVect *LignG = pnl_vect_create(path->n);

    for(int j = 0; j < path->m-1; ++j){
        pnl_mat_get_row(LignG, Gmat, j);
        for(int i = 0; i < path->n; ++i){
            double LG = 0.0;
            pnl_mat_get_row(LignL, L, i);
            
            LG = pnl_vect_scalar_prod(LignL, LignG);
            
            double Si = GET(this->sigma_, i);
            
            MLET(path, j+1,i) = MGET(path,j , i) * exp((freeRate - Si * Si /2.0) * step + Si*sqrtstep*LG  );


        }
    }

    pnl_mat_free(&Gmat);
    pnl_vect_free(&LignG);


}

void BlackScholesModel :: asset(PnlMat *path, double t, double T, int nbTimeSteps, PnlRng *rng, const PnlMat *past) {
    double step = T/(double)nbTimeSteps;
    bool tisdiscretisation = std::fmod(t, step) == 0;
    double interstep = 0;
    int simuremains = 0;
    // On distingue le cas t est un temps de discrétisation ou non
    if (tisdiscretisation) simuremains = path->m - past->m;
    else simuremains = path->m - past->m + 1;
    
    
    PnlVect SEnd = pnl_vect_wrap_mat_row(past, past->m-1);
    PnlMat *L = this->El;
    PnlVect *LignL = this->LignEl;
    double freeRate = this->r_;

    // Création de la matrice G randomisée
    PnlMat *GMat = pnl_mat_create(simuremains, path->n);
    pnl_mat_rng_normal(GMat, simuremains, path->n, rng);
    PnlVect *GLign = pnl_vect_create(path->n);
    PnlVect *StPrevious = pnl_vect_create_from_scalar(path->n, 1);
    
    PnlVect *S = pnl_vect_create(path->n);

    // Pour chaque temps qui reste, on fait :
    for (int j = 0; j < simuremains; j++) {
        pnl_mat_get_row(GLign, GMat, j);
        // Si t est un temps de discrétisation alors les pas sont tout le temps égaux a step
        // sinon, le premier pas est la distance entre t et le prochain pas de discrétisation
        if (tisdiscretisation) interstep = step;
        else interstep = step - std::fmod(t, step);
        
        // Pour chaque actif
        for(int i = 0; i < path->n; i++) {

            double LG = 0;
            //PnlVect testL = pnl_vect_wrap_mat_row(L, j);
            pnl_mat_get_row(LignL, L, i);
            LG = pnl_vect_scalar_prod(GLign, LignL);
            //PnlVect testG = pnl_vect_wrap_mat_row(GMat, s);
            //LG = pnl_vect_scalar_prod(&testG, &testL);

            if(j!= 0)  {
                interstep = step;
            }
            double Si = GET(this->sigma_, i);
            LET(S, i) = GET(StPrevious, i) * exp((freeRate - Si* Si /2.0) * (interstep) + Si*sqrt(interstep)*LG  );
            LET(StPrevious, i) = GET(S, i);
            LET(S, i) = GET(&SEnd, i) * GET(S, i);


        }
        
        //pnl_vect_mult_scalar(S, si);
        
        for (int i = 0; i < S->size; i++) {
            MLET(path, path->m - simuremains + j,i) = GET(S, i);
        }
    }
    pnl_vect_free(&StPrevious);
    pnl_mat_free(&GMat);
    pnl_vect_free(&GLign);

}

void BlackScholesModel :: shiftAsset(PnlMat *shift_path, const PnlMat *path, int d, double h, double t, double timestep){
    //Nous n'utilisons pas MGET et GET car ils prennent plus de temps que ->array
    for(int i = 0; i < path->m ; i++){
        if(i*timestep < t){
            MLET(shift_path, i, d) = MGET(path, i, d);
        }else {
            MLET(shift_path, i, d) = (1.0+h)*MGET(path, i, d);
        }
    }
}


void BlackScholesModel :: simul_market(PnlMat *path, double T, int nbTimeSteps, PnlRng *rng){
    PnlMat *L = this->El;
    PnlVect *LignL = this->LignEl;
 
    pnl_mat_set_col(path, this->spot_, 0);

    double step = T/(double)nbTimeSteps;
    double sqrtstep = sqrt(step);

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
            MLET(path, j,i+1) = this->spot_->array[j] * exp((this->Trend->array[j] - Sj * Sj /2.0) * step * i + Sj*LG  );


        }
    }

    pnl_mat_free(&Gmat);
    pnl_vect_free(&LignG);
}

void BlackScholesModel :: setTrend(double trend) {
    this->Trend = pnl_vect_create_from_scalar(this->size_, trend);
}