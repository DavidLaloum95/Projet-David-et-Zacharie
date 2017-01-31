///////////////////////////////////////////////////////////////////////////////////////
////             BlackScholes.cpp : Résultat théorique attendus                   /////
///////////////////////////////////////////////////////////////////////////////////////


#include "BlackScholes.hpp"
#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>


// Densité de la loi normale

double norm_pdf(const double& x) {
    return (1.0/(pow(2*M_PI,0.5)))*exp(-0.5*x*x);
}

// Fonction de répartition de la loi normale

double norm_cdf(const double& x) {
    double k = 1.0/(1.0 + 0.2316419*x);
    double k_sum = k*(0.319381530 + k*(-0.356563782 + k*(1.781477937 + k*(-1.821255978 + 1.330274429*k))));
    
    if (x >= 0.0) {
        return (1.0 - (1.0/(pow(2*M_PI,0.5)))*exp(-0.5*x*x) * k_sum);
    } else {
        return 1.0 - norm_cdf(-x);
    }
    
}

// calcul de d_j, pour j dans {1,2}.

double d_j(const int& j, const double& S, const double& K, const double& r, const double& v, const double& T) {
    return (log(S/K) + (r + (pow(-1,j-1))*0.5*v*v)*T)/(v*(pow(T,0.5)));
}



// Calcul  d'une option vanille Européenne put price, basé sur
// prix de l'action S, strike K, taux d'interet r, volatilité sigma
// et le temps de maturité T

double put_price(const double& S, const double& K, const double& r, const double& v, const double& T) {
    return -S*norm_cdf(-d_j(1, S, K, r, v, T))+K*exp(-r*T) * norm_cdf(-d_j(2, S, K, r, v, T));
}

