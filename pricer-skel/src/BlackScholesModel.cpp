#include "BlackScholesModel.hpp"

void BlackScholesModel :: asset(PnlMat *path, double T, int nbTimeSteps, PnlRng *rng){

    PnlMat *L = this->El;
    PnlVect *LignL = this->LignEl;
 
    pnl_mat_set_col(path, this->spot_, 0);

    double step = T/(double)nbTimeSteps;
    double sqrtstep = sqrt(step);
    double freeRate = this->r_;
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
            
            MLET(path, j,i+1) = MGET(path,j , i) * exp((freeRate - Sj * Sj /2.0) * step + Sj*sqrtstep*LG  );


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
    if (tisdiscretisation) simuremains = path->n - past->n;
    else simuremains = path->n - past->n + 1;

    

    PnlMat *L = this->El;
    PnlVect *LignL = this->LignEl;
    double freeRate = this->r_;

    // Création de la matrice G randomisée
    PnlMat *GMat = pnl_mat_create(simuremains, path->m);
    pnl_mat_rng_normal(GMat, simuremains, path->m, rng);
    PnlVect *GLign = pnl_vect_create(path->m);

    // Pour chaque actif, on fait :
    for (int j = 0; j < path->m; j++) {
        
        PnlVect *S = pnl_vect_create(simuremains);
        double previous = 1.0;
        double sj = past->array[(j+1)*past->n - 1];
        //PnlVect testL = pnl_vect_wrap_mat_row(L, j);
        pnl_mat_get_row(LignL, L, j);


        // Si t est un temps de discrétisation alors les pas sont tout le temps égaux a step
        // sinon, le premier pas est la distance entre t et le prochain pas de discrétisation
        if (tisdiscretisation) interstep = step;
        else interstep = step - std::fmod(t, step);
        
        for(int s = 0; s < S->size; s++) {

            double LG = 0;
            pnl_mat_get_row(GLign, GMat, s);
            LG = pnl_vect_scalar_prod(GLign, LignL);

            //PnlVect testG = pnl_vect_wrap_mat_row(GMat, s);
            //LG = pnl_vect_scalar_prod(&testG, &testL);

            
            if(s != 0)  {
                previous = GET(S, s-1);
                interstep = step;
            }
            double Sj = GET(this->sigma_, j);
            LET(S, s) = previous * exp((freeRate - Sj* Sj /2.0) * (interstep) + Sj*sqrt(interstep)*LG  );

        }
        pnl_vect_mult_scalar(S, sj);

        for (int s = 0; s < S->size; s++) {
            path->array[(j+1)*path->n - simuremains + s] = S->array[s];
        }
    }

    pnl_mat_free(&GMat);
    pnl_vect_free(&GLign);    
}

void BlackScholesModel :: shiftAsset(PnlMat *shift_path, const PnlMat *path, int d, double h, double t, double timestep){
    

    //Nous n'utilisons pas MGET et GET car ils prennent plus de temps que ->array
    for(int i = 0; i < path->n ; i++){
        if(i*timestep < t){
            shift_path->array[d*path->n + i] = path->array[d*path->n + i];
        }else {
            shift_path->array[d*path->n + i] = (1.0+h)*path->array[d*path->n + i];
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