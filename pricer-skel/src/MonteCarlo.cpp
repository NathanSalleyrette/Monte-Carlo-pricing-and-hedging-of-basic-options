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
    PnlMat *path = pnl_mat_create(this->opt_->nbTimeSteps_ + 1, this->opt_->size_);
    for(int i = 0; i< this->nbSamples_; i++){
        this->mod_->asset(path, this->opt_->T_ , this->opt_->nbTimeSteps_, this->rng_);
        respayoff = this->opt_->payoff(path);
        sumpayoff += respayoff;
        squaresum += respayoff * respayoff;
        
        PnlMat *shiftpath = pnl_mat_copy(path);


        for(int d = 0; d < path->n; d++){
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
    PnlMat *path = pnl_mat_create(this->opt_->nbTimeSteps_ + 1, this->opt_->size_);
    pnl_vect_set_zero(this->sumShift);
    pnl_vect_set_zero(this->sumShiftSquare);

    for(int i = 0; i< sample; i++){
        pnl_mat_set_subblock(path, past, 0, 0);
        this->mod_->asset(path, t, this->opt_->T_ , this->opt_->nbTimeSteps_, this->rng_, past);
        
        respayoff = this->opt_->payoff(path);
        price += respayoff;
        squaresum += respayoff * respayoff;

        PnlMat *shiftpath = pnl_mat_copy(path);
        
        // Pour tous les actifs
        for(int d = 0; d < path->n; d++){

            this->mod_->shiftAsset(shiftpath, path, d, this->fdStep_, t, timestep);
            double payoffhplus = this->opt_->payoff(shiftpath);
            LET(this->sumShift,d) = GET(this->sumShift,d)+payoffhplus;
            
            this->mod_->shiftAsset(shiftpath, path, d, -this->fdStep_, t, timestep);
            double payoffhmoins = this->opt_->payoff(shiftpath);
            LET(this->sumShift,d) = GET(this->sumShift,d)-payoffhmoins;
            LET(this->sumShiftSquare,d) = GET(this->sumShiftSquare, d) + (payoffhplus-payoffhmoins)*(payoffhplus-payoffhmoins);
            this->mod_->shiftAsset(shiftpath, path, d, 0.0, t, timestep);

            
        }

        pnl_mat_free(&shiftpath);
    }

    double rTt = this->mod_->r_ * (this->opt_->T_ - t);
    prix =  price*exp(- rTt) / sample;
    sigma = squaresum / sample - (price / sample)*(price / sample);
    
    // Si t = maturité alors le payoff est calculé avec exactitude
    if (t == this->opt_->T_) sigma = 0.0;
    
    sigma = sigma*exp(- 2 * rTt);
    std_dev = sqrt(sigma) / sqrt(sample);

    pnl_mat_free(&path);
}


void MonteCarlo :: delta(PnlVect *delta, PnlVect *std_dev){
    // Calcul de la formule mathématique (4) donnée à la page 4 de l'énoncé.
    double sample = this->nbSamples_;
    double rT = this->mod_->r_ * this->opt_->T_;

    for(int i = 0; i < this->sumShift->size; i++){
        
        double denom = 2*GET(this->mod_->spot_, i)*this->fdStep_;
        double sumshifti = GET(this->sumShift, i);
        double sumshiftsquarei = GET(this->sumShiftSquare,i);

        LET(delta,i) = exp( - rT)/(sample * denom) * sumshifti;
        
        double sigma  = sumshiftsquarei /(denom*denom*sample) - (sumshifti /(denom*sample))*(sumshifti/(denom*sample));
        
        sigma *= exp(-2 * rT);
        
        LET(std_dev,i) = sqrt(sigma) / sqrt(sample);
    }    

}


void MonteCarlo :: delta(const PnlMat *past, double t, PnlVect *delta, PnlVect *std_dev){
    
    double sample = this->nbSamples_;
    double rTt = this->mod_->r_ * (this->opt_->T_ - t);

    for(int i = 0; i < this->sumShift->size; i++){
        
        double denom = 2*MGET( past ,past->m - 1, i)*this->fdStep_;
        double sumshifti = GET(this->sumShift, i);
        double sumshiftsquarei = GET(this->sumShiftSquare,i);

        LET(delta,i) = exp( - rTt)/(sample* denom) * sumshifti;
        
        double sigma  = sumshiftsquarei /(denom*denom*sample) - (sumshifti /(denom*sample))*(sumshifti/(denom*sample));
        
        sigma *= exp(-2 * rTt);
        
        LET(std_dev,i) = sqrt(sigma) / sqrt(sample);
    }
}

void MonteCarlo :: portfolioValue(PnlVect *delta, double &value, double &price){
    value = price;

    for(int i = 0; i < delta->size; i++){
        value -= GET(delta,i) * GET(this->mod_->spot_, i);
    }
}


void MonteCarlo :: portfolioValue(PnlVect *deltaNow, PnlVect *deltaLast, double &value, double &H, PnlVect *spots){
    value *= exp(this->mod_->r_ * this->opt_->T_/H);

    for(int i = 0; i < deltaNow->size; i++){
        value -= (GET(deltaNow,i) - GET(deltaLast,i)) * GET(spots,i);
    }
}

void MonteCarlo :: profitAndLoss(PnlMat *pathReal, double &error, double &prix0, double &std_dev0){
    
    // Calcul de value initial
    double prix = 0.0;
    double std_dev = 0.0;
    double stepHedging = this->opt_->T_/(pathReal->m - 1);
    double paststep = this->opt_->T_ / this->opt_->nbTimeSteps_;
    int step = (int) (pathReal->m - 1)/this->opt_->nbTimeSteps_;
    double value = 0.0;
    double H = pathReal->m - 1;
    PnlVect *deltaNow = pnl_vect_create(pathReal->n);
    PnlVect *deltaLast = pnl_vect_create(pathReal->n);
    PnlVect *std_dev_delta = pnl_vect_create(pathReal->n);
    PnlVect *lignSpots = pnl_vect_create(pathReal->n);
    price(prix0, std_dev0);
    delta(deltaNow, std_dev_delta);
    portfolioValue(deltaNow, value, prix0);

    PnlMat *pastfull = pnl_mat_create(this->opt_->nbTimeSteps_ + 1, pathReal->n);
    PnlMat *past = pnl_mat_create(1, pathReal->n);

    for(int i = 0; i <= this->opt_->nbTimeSteps_; i++ ){
        pnl_mat_get_row(lignSpots, pathReal, (int)(step*i));
        pnl_mat_set_row(pastfull, lignSpots,i);
    }

    
    // Prix en 0 déjà calculé plus haut, donc on commence au premier élément
    for(int i = 1 ; i < pathReal->m ;  i++){
        int end = (int) (i*stepHedging)/paststep;
        pnl_mat_extract_subblock(past, pastfull, 0, end+1, 0, pathReal->n);
        // Calcul des différentes valeurs de portefeuille
        price( past, i * stepHedging, prix, std_dev);
        deltaLast = pnl_vect_copy(deltaNow);
        delta(past,i * stepHedging ,deltaNow, std_dev_delta);
        pnl_mat_get_row(lignSpots,pathReal, i);
        portfolioValue(deltaNow,deltaLast,value, H , lignSpots);
    }

    //Calcul de l'erreur
    error = value;

    for(int i = 0; i < deltaNow->size   ; i++){
        error += deltaNow->array[i] * MGET(pathReal, pathReal->m-1, i);
    }

    error -= this->opt_->payoff(pastfull);

    pnl_vect_free(&deltaNow);
    pnl_vect_free(&deltaLast);
    pnl_vect_free(&std_dev_delta);
    pnl_vect_free(&lignSpots);
    pnl_mat_free(&pastfull);
    pnl_mat_free(&past);

}