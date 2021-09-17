#pragma once

#include "Option.hpp"
#include "BlackScholesModel.hpp"
#include "pnl/pnl_random.h"
class MonteCarlo
{
public:
    BlackScholesModel *mod_; /*! pointeur vers le modèle */
    Option *opt_; /*! pointeur sur l'option */
    PnlRng *rng_; /*! pointeur sur le générateur */
    double fdStep_; /*! pas de différence finie */
    int nbSamples_; /*! nombre de tirages Monte Carlo */
    PnlVect *sumShift;
    PnlVect *sumShiftSquare; 


    MonteCarlo(BlackScholesModel *mod_, Option *opt_, PnlRng *rng_, double fdStep_, int nbSamples_)
        : mod_(mod_),
        opt_(opt_),
        rng_(rng_),
        fdStep_(fdStep_),
        nbSamples_(nbSamples_)

    {
        this->sumShift = pnl_vect_create_from_zero(opt_->size_);
        this->sumShiftSquare = pnl_vect_create_from_zero(opt_->size_);
    }

    ~MonteCarlo() {
        pnl_vect_free(&sumShift);
        pnl_vect_free(&sumShiftSquare);
    }
    
    MonteCarlo(const MonteCarlo &MonteCarlo)
        : mod_(MonteCarlo.mod_),
        opt_(MonteCarlo.opt_),
        rng_(MonteCarlo.rng_),
        fdStep_(MonteCarlo.fdStep_),
        nbSamples_(MonteCarlo.nbSamples_)
    {}

    /**
     * Calcule le prix de l'option à la date 0
     *
     * @param[out] prix valeur de l'estimateur Monte Carlo
     * @param[out] ic écart type de l'estimateur
     */
    void price(double &prix, double &std_dev);

    /**
     * Calcule le prix de l'option à la date t
     *
     * @param[in]  past contient la trajectoire du sous-jacent
     * jusqu'à l'instant t
     * @param[in] t date à laquelle le calcul est fait
     * @param[out] prix contient le prix
     * @param[out] std_dev contient l'écart type de l'estimateur
     */
    void price(const PnlMat *past, double t, double &prix, double &std_dev);

    /**
     * Calcule le delta de l'option à la date t
     *
     * @param[in] past contient la trajectoire du sous-jacent
     * jusqu'à l'instant t
     * @param[in] t date à laquelle le calcul est fait
     * @param[out] delta contient le vecteur de delta
     * @param[out] std_dev contient l'écart type de l'estimateur
     */
    void delta(const PnlMat *past, double t, PnlVect *delta, PnlVect *std_dev);

    /**
     * Calcule le delta de l'option à la date 0
     *
     * @param[in] t date à laquelle le calcul est fait
     * @param[out] delta contient le vecteur de delta
     * @param[out] std_dev contient l'écart type de l'estimateur
     */
    void delta(PnlVect *delta, PnlVect *std_dev);


    /**
     * Calcul du portefeuille de couverture à la date 0
     * @param[in] delta contient le vecteur de delta en 0
     * @param[out] value contient la valeur du portefeuille à la date 0
     */
    double portfolioValue(PnlVect *delta);

    /**
     * Calcul du portefeuille à une date t > 0
     * @param[in] delta contient le vecteur des delta pour la trajectoire de marché
     * @param[in] path contient la trajactoire de marché
     * @param[out] value valeur du portefeuille à l'instant t
     */
    double portfolioValue(PnlMat *delta, PnlMat *path);


};


