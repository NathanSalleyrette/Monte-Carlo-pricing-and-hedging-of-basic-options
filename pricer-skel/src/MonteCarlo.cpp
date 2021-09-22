#include "MonteCarlo.hpp"
#include "math.h"


void MonteCarlo::price(double &prix, double &std_dev){
    //On simule la trajectoire
    double sumpayoff = 0.0;
    double squaresum = 0.0;
    double respayoff = 0.0;
    double sigma = 0.0;
    double sample = this->nbSamples_;
    pnl_vect_set_zero(this->sumShift);
    pnl_vect_set_zero(this->sumShiftSquare);


    double timestep = this->opt_->T_ / this->opt_->nbTimeSteps_;
    PnlMat *path = pnl_mat_create(this->opt_->size_, this->opt_->nbTimeSteps_ + 1);
    for(int i = 0; i< this->nbSamples_; i++){
        this->mod_->asset(path, this->opt_->T_ , this->opt_->nbTimeSteps_, this->rng_);
        respayoff = this->opt_->payoff(path);
        sumpayoff += respayoff;
        squaresum += respayoff * respayoff;
        
        PnlMat *shiftpath = pnl_mat_copy(path);


        for(int d = 0; d < path->m; d++){
            // PnlMat *shiftpath = pnl_mat_create(this->opt_->size_, this->opt_->nbTimeSteps_+1);
            // pnl_mat_clone(shift_path, path);
            this->mod_->shiftAsset(shiftpath, path, d, this->fdStep_, 0.0, timestep);
            double payoffhplus = this->opt_->payoff(shiftpath);
            LET(this->sumShift,d) = GET(this->sumShift,d)+payoffhplus;
            
            this->mod_->shiftAsset(shiftpath, path, d, -this->fdStep_, 0.0, timestep);
            double payoffhmoins = this->opt_->payoff(shiftpath);
            LET(this->sumShift,d) = GET(this->sumShift,d)-payoffhmoins;
            LET(this->sumShiftSquare,d) = GET(this->sumShiftSquare, d) + (payoffhplus-payoffhmoins)*(payoffhplus-payoffhmoins);
            this->mod_->shiftAsset(shiftpath, path, d, 0.0, 0.0, timestep);

            
        }
        //On détruit les objets inutiles
        pnl_mat_free(&shiftpath);


    }
    pnl_mat_free(&path);
    //pnl_mat_free(&shiftpath);


    prix =  sumpayoff*exp(- this->mod_->r_ * this->opt_->T_) / sample;

    sigma = squaresum / sample - (sumpayoff / sample)*(sumpayoff / sample);
    sigma = sigma*exp(- 2 * this->mod_->r_ * this->opt_->T_);
    std_dev = sqrt(sigma) / sqrt(sample);
}


void MonteCarlo::price(const PnlMat *past, double t, double &prix, double &std_dev) {
    double price = 0.0;
    double squaresum = 0.0;
    double respayoff = 0.0;
    double sigma = 0.0;
    double sample = this->nbSamples_;
    double timestep = this->opt_->T_ / this->opt_->nbTimeSteps_;
    PnlMat *path = pnl_mat_create(this->opt_->size_, this->opt_->nbTimeSteps_ + 1);
    pnl_vect_set_zero(this->sumShift);
    pnl_vect_set_zero(this->sumShiftSquare);

    for(int i = 0; i< this->nbSamples_; i++){
        pnl_mat_set_subblock(path, past, 0, 0);
        this->mod_->asset(path, t, this->opt_->T_ , this->opt_->nbTimeSteps_, this->rng_, past);
        
        respayoff = this->opt_->payoff(path);
        price += respayoff;
        squaresum += respayoff * respayoff;

        PnlMat *shiftpath = pnl_mat_copy(path);
        
        // Pour tous les actifs
        for(int d = 0; d < path->m; d++){
            // PnlMat *shiftpath = pnl_mat_create(this->opt_->size_, this->opt_->nbTimeSteps_+1);
            // pnl_mat_clone(shift_path, path);
            this->mod_->shiftAsset(shiftpath, path, d, this->fdStep_, t, timestep);
            double payoffhplus = this->opt_->payoff(shiftpath);
            LET(this->sumShift,d) = GET(this->sumShift,d)+payoffhplus;
            
            this->mod_->shiftAsset(shiftpath, path, d, -this->fdStep_, t, timestep);
            double payoffhmoins = this->opt_->payoff(shiftpath);
            LET(this->sumShift,d) = GET(this->sumShift,d)-payoffhmoins;
            LET(this->sumShiftSquare,d) = GET(this->sumShiftSquare, d) + (payoffhplus-payoffhmoins)*(payoffhplus-payoffhmoins);
            this->mod_->shiftAsset(shiftpath, path, d, 0.0, t, timestep);

            
        }
    }

    prix =  price*exp(- this->mod_->r_ * (this->opt_->T_ - t)) / this->nbSamples_;
    sigma = squaresum / this->nbSamples_ - (price / this->nbSamples_)*(price / this->nbSamples_);
    if (t == this->opt_->T_) sigma = 0.0;
    sigma = sigma*exp(- 2 * this->mod_->r_ * (this->opt_->T_ - t));
    std_dev = sqrt(sigma) / sqrt(this->nbSamples_);
}


void MonteCarlo :: delta(PnlVect *delta, PnlVect *std_dev){

    for(int i = 0; i < this->sumShift->size; i++){
        LET(delta,i) = exp( - this->mod_->r_ * this->opt_->T_)/(this->nbSamples_ * 2 * this->fdStep_ * GET(this->mod_->spot_,i)) * GET(this->sumShift,i);
        
        double denom = 2*GET(this->mod_->spot_, i)*this->fdStep_;
        double sigma  = GET(this->sumShiftSquare,i) /(denom*denom*this->nbSamples_) - (GET(this->sumShift,i) /(denom*this->nbSamples_))*(GET(this->sumShift,i)/(denom*this->nbSamples_));
        sigma *= exp(-2 * this->mod_->r_ * this->opt_->T_);
        LET(std_dev,i) = sqrt(sigma) / sqrt(this->nbSamples_);
    }    

}


void MonteCarlo :: delta(const PnlMat *past, double t, PnlVect *delta, PnlVect *std_dev){
    for(int i = 0; i < this->sumShift->size; i++){
        LET(delta,i) = exp( - this->mod_->r_ * (this->opt_->T_ - t))/(this->nbSamples_ * 2 * this->fdStep_ * MGET( past ,i, past->n - 1)) * GET(this->sumShift,i);
        
        double denom = 2*GET(this->mod_->spot_, i)*this->fdStep_;
        double sigma  = GET(this->sumShiftSquare,i) /(denom*denom*this->nbSamples_) - (GET(this->sumShift,i) /(denom*this->nbSamples_))*(GET(this->sumShift,i)/(denom*this->nbSamples_));
        sigma *= exp(-2 * this->mod_->r_ * (this->opt_->T_ - t));
        LET(std_dev,i) = sqrt(sigma) / sqrt(this->nbSamples_);
    }
}

void MonteCarlo :: portfolioValue(PnlVect *delta, double &value, double &price){
    value = price;

    for(int i = 0; i < delta->size; i++){
        value -= delta->array[i] * this->mod_->spot_->array[i];
    }
}


void MonteCarlo :: portfolioValue(PnlVect *deltaNow, PnlVect *deltaLast, double &value, double &H, PnlVect *spots){
    value *= exp(this->mod_->r_ * this->opt_->T_/H);

    for(int i = 0; i < deltaNow->size; i++){
        value -= (deltaNow->array[i] - deltaLast->array[i]) * spots->array[i];
    }
}

void MonteCarlo :: profitAndLoss(PnlMat *pathReal, double &error, double &prix0, double &std_dev0){
    
    // Calcul de value initial
    
    double prix = 0.0;
    double std_dev = 0.0;
    double stepHedging = this->opt_->T_/(pathReal->n - 1);
    double paststep = this->opt_->T_ / this->opt_->nbTimeSteps_;
    int step = (int) (pathReal->n - 1)/this->opt_->nbTimeSteps_;
    double value = 0.0;
    double H = pathReal->n - 1;
    PnlVect *deltaNow = pnl_vect_create(pathReal->m);
    PnlVect *deltaLast = pnl_vect_create(pathReal->m);
    PnlVect *std_dev_delta = pnl_vect_create(pathReal->m);
    PnlVect *lignSpots = pnl_vect_create(pathReal->m);
    price(prix0, std_dev0);
    delta(deltaNow, std_dev_delta);
    portfolioValue(deltaNow, value, prix0);

    PnlMat *pastfull = pnl_mat_create(pathReal->m, this->opt_->nbTimeSteps_ + 1);
    PnlMat *past = pnl_mat_create(pathReal->m, 1);

    for(int i = 0; i <= this->opt_->nbTimeSteps_; i++ ){
        pnl_mat_get_col(lignSpots, pathReal, (int)(step*i));
        pnl_mat_set_col(pastfull, lignSpots,i);
    }

    
    // Prix en 0 déjà calculé plus haut, donc on commence au premier élément
    for(int i = 1 ; i < pathReal->n ;  i++){
        
        
        int end = (int) (i*stepHedging)/paststep;
        pnl_mat_extract_subblock(past, pastfull, 0, pathReal->m, 0, end+1);
        // Calcul des différentes valeurs de portefeuille
        price( past, i * stepHedging, prix, std_dev);
        deltaLast = pnl_vect_copy(deltaNow);
        delta(past,i * stepHedging ,deltaNow, std_dev_delta);
        pnl_mat_get_col(lignSpots,pathReal, i);
        portfolioValue(deltaNow,deltaLast,value, H , lignSpots);

        if (i == pathReal->n - 2) {
            double r = 0.0;
        }

    }

    //Calcul de l'erreur
    error = value;

    for(int i = 0; i < deltaNow->size   ; i++){
        error += deltaNow->array[i] * MGET(pathReal, i, pathReal->n-1);

    }

    error -= this->opt_->payoff(pastfull);
}