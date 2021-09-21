#pragma once

#include "pnl/pnl_random.h"
#include "pnl/pnl_vector.h"
#include "pnl/pnl_matrix.h"

/// \brief Modèle de Black Scholes
class BlackScholesModel
{
public:
    int size_; /// nombre d'actifs du modèle
    double r_; /// taux d'intérêt
    double rho_; /// paramètre de corrélation
    PnlVect *sigma_; /// vecteur de volatilités
    PnlVect *spot_; /// valeurs initiales des sous-jacents
    PnlMat *El; /// Matrice triangulaire inférieure
    PnlVect *LignEl; /// Ligne de El
    PnlVect *Trend;

    BlackScholesModel(int size_, double r_, double rho_,PnlVect *sigma_, PnlVect *spot_)
        : size_(size_)
        , r_(r_)
        , rho_(rho_)
        , sigma_(sigma_)
        , spot_(spot_)
    { 
        this->El = pnl_mat_create_from_scalar(size_, size_, rho_);
        this->LignEl = pnl_vect_create(size_);
        for(int i = 0; i < size_; i++){
            pnl_mat_set_diag(this->El, 1.0, i);
        }

        pnl_mat_chol(this->El);
    }

    ~BlackScholesModel() { 
        pnl_mat_free(&El);
        pnl_vect_free(&LignEl);
    }

    BlackScholesModel(const BlackScholesModel &BlackScholesModel)
        : size_(BlackScholesModel.size_)
        , r_(BlackScholesModel.r_)
        , rho_(BlackScholesModel.rho_)
        , sigma_(BlackScholesModel.sigma_)
        , spot_(BlackScholesModel.spot_)
    { }

    /**
     * Génère une trajectoire du modèle et la stocke dans path
     *
     * @param[out] path contient une trajectoire du modèle.
     * C'est une matrice de taille (nbTimeSteps+1) x d
     * @param[in] T  maturité
     * @param[in] nbTimeSteps nombre de dates de constatation
     */
    void asset(PnlMat *path, double T, int nbTimeSteps, PnlRng *rng);
    

    /**
     * Calcule une trajectoire du modèle connaissant le
     * passé jusqu' à la date t
     *
     * @param[out] path  contient une trajectoire du sous-jacent
     * donnée jusqu'à l'instant t par la matrice past
     * @param[in] t date jusqu'à laquelle on connait la trajectoire.
     * t n'est pas forcément une date de discrétisation
     * @param[in] nbTimeSteps nombre de pas de constatation
     * @param[in] T date jusqu'à laquelle on simule la trajectoire
     * @param[in] past trajectoire réalisée jusqu'a la date t
     */
    void asset(PnlMat *path, double t, double T, int nbTimeSteps, PnlRng *rng, const PnlMat *past);

    /**
     * Shift d'une trajectoire du sous-jacent
     *
     * @param[in]  path contient en input la trajectoire
     * du sous-jacent
     * @param[out] shift_path contient la trajectoire path
     * dont la composante d a été shiftée par (1+h)
     * à partir de la date t.
     * @param[in] t date à partir de laquelle on shift
     * @param[in] h pas de différences finies
     * @param[in] d indice du sous-jacent à shifter
     * @param[in] timestep pas de constatation du sous-jacent
     */
    void shiftAsset(PnlMat *shift_path, const PnlMat *path, int d, double h, double t, double timestep);


    void simul_market(PnlMat *path, double T, int nbTimeSteps, PnlRng *rng);

    void setTrend(double trend);

};


