#include "MonteCarlo.hpp"
#include "math.h"


void MonteCarlo::price(double &prix, double &std_dev){
    //On simule la trajectoire
    double price = 0.0;
    double squaresum = 0.0;
    double respayoff = 0.0;
    double sigma = 0.0;
    for(int i = 0; i< this->nbSamples_; i++){
        PnlMat *path = pnl_mat_create(3, this->opt_->nbTimeSteps_+1);
        this->mod_->asset(path, this->opt_->T_ , this->opt_->nbTimeSteps_, this->rng_);
        respayoff = this->opt_->payoff(path);
        price += respayoff;
        squaresum += respayoff * respayoff;
        
        //On détruit les objets inutiles
        pnl_mat_free(&path);


    }
    prix =  price*exp(- this->mod_->r_ * this->opt_->T_) / this->nbSamples_;
    sigma = squaresum / this->nbSamples_ - (price / this->nbSamples_)*(price / this->nbSamples_);
    sigma = sigma*exp(- 2 * this->mod_->r_ * this->opt_->T_);
    std_dev = 1.96 * sqrt(sigma) / sqrt(this->nbSamples_);
}