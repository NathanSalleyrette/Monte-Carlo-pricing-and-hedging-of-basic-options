#include "MonteCarlo.hpp"


void MonteCarlo::price(double &prix, double &std_dev){
    //On simule la trajectoire
    int M = 1000;
    double price = 0.0;
    for(int i = 0; i< M; i++){
        PnlMat *path = pnl_mat_create(3, this->opt_->nbTimeSteps_+1);
        this->mod_->asset(path, this->opt_->T_ , this->opt_->nbTimeSteps_, this->rng_);
        price += this->opt_->payoff(path);
    }
    prix =  price*exp(- this->mod_->r_ * this->opt_->T_) / M;
    std_dev = 1.96 * 2.0;
}