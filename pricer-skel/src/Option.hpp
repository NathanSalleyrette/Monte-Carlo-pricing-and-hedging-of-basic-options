#pragma once

#include "pnl/pnl_vector.h"
#include "pnl/pnl_matrix.h"

/// \brief Classe Option abstraite
class Option
{
public:
    double T_; /// maturité
    int nbTimeSteps_; /// nombre de pas de temps de discrétisation
    int size_; /// dimension du modèle, redondant avec BlackScholesModel::size_
    double strike;
    Option(double T_, int nbTimeSteps_, int size_, double strike)
        :T_(T_)
        ,nbTimeSteps_(nbTimeSteps_)
        ,size_(size_)
        ,strike(strike)
    { }

    ~Option() { }

    Option(const Option &Opt)
        :T_(Opt.T_)
        ,nbTimeSteps_(Opt.nbTimeSteps_)
        ,size_(Opt.size_)
        ,strike(Opt.strike)
    { }



    /**
     * Calcule la valeur du payoff sur la trajectoire
     *
     * @param[in] path est une matrice de taille (N+1) x d
     * contenant une trajectoire du modèle telle que créée
     * par la fonction asset. (N = nbTimeSteps_) (d = nombre d'actifs = size_) 
     * @return phi(trajectoire)
     */
    virtual double payoff(const PnlMat *path) = 0;
};


