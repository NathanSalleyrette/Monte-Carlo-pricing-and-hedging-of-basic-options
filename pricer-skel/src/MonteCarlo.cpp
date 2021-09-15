#include "MonteCarlo.hpp"
#include "math.h"


void MonteCarlo::price(double &prix, double &std_dev){
    //On simule la trajectoire
    double price = 0.0;
    double squaresum = 0.0;
    double respayoff = 0.0;
    double sigma = 0.0;
    double timestep = this->opt_->T_ / this->opt_->nbTimeSteps_;
    for(int i = 0; i< this->nbSamples_; i++){
        PnlMat *path = pnl_mat_create(this->opt_->size_ , this->opt_->nbTimeSteps_+1);
        this->mod_->asset(path, this->opt_->T_ , this->opt_->nbTimeSteps_, this->rng_);
        respayoff = this->opt_->payoff(path);
        price += respayoff;
        squaresum += respayoff * respayoff;
        PnlMat *shiftpath = pnl_mat_create(this->opt_->size_, this->opt_->nbTimeSteps_+1);
        for(int d = 0; d < path->n; d++){
            this->mod_->shiftAsset(shiftpath, path, d, this->fdStep_, 0.0, timestep);
            this->sumShift->array[d] += this->opt_->payoff(shiftpath);
            this->mod_->shiftAsset(shiftpath, path, d, -this->fdStep_, 0.0, timestep);
            this->sumShift->array[d] -= this->opt_->payoff(shiftpath);
            this->sumShiftSquare->array[d] = this->sumShift->array[d] * this->sumShift->array[d];
            
        }
        //On détruit les objets inutiles
        pnl_mat_free(&path);
        pnl_mat_free(&shiftpath);


    }
    prix =  price*exp(- this->mod_->r_ * this->opt_->T_) / this->nbSamples_;
    sigma = squaresum / this->nbSamples_ - (price / this->nbSamples_)*(price / this->nbSamples_);
    sigma = sigma*exp(- 2 * this->mod_->r_ * this->opt_->T_);
    std_dev = 1.96 * sqrt(sigma) / sqrt(this->nbSamples_);
}

void MonteCarlo :: delta(PnlVect *delta, PnlVect *std_dev){

    for(int i = 0; i < this->sumShift->size; i++){
        delta->array[i] = exp( - this->mod_->r_ * this->opt_->T_)/(this->nbSamples_ * 2 * this->fdStep_ * this->mod_->spot_->array[i]) * this->sumShift->array[i];
        double sigma  = this->sumShiftSquare->array[i]/this->nbSamples_ - (this->sumShiftSquare->array[i] /this->nbSamples_)*(this->sumShiftSquare->array[i]/ this->nbSamples_);
        sigma *= exp(-2 * this->mod_->r_ * this->opt_->T_);
        std_dev->array[i] = 1.96 * sqrt(sigma) / sqrt(this->nbSamples_);
    }    

}