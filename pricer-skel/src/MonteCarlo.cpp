#pragma once

#include "Option.hpp"
#include "BlackScholesModel.hpp"
#include "pnl/pnl_random.h"
#include "MonteCarlo.hpp"


void MonteCarlo::price(double &prix, double &std_dev){
    //On simule la trajectoire
    int M = 1000;
    double price = 0.0;
    for(int i < 0; i< M; i++){
        PnlMat *path = pnl_mat_new();
        this->mod_->asset(path, this->opt_->T_ , this->nbTimeSteps, this->rng);
        price += this->opt_->payoff(path);
    }
    *price =  price*exp(- this->mod_->r_ * this->opt_->T_) / M;
    *std_dev = 1.96 * 2;
}