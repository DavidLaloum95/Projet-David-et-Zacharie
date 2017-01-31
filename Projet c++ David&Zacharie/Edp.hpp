///////////////////////////////////////////////////////////////////////////////////////
////                            Fichier EDP.h                                     /////
///////////////////////////////////////////////////////////////////////////////////////

#ifndef Edp_hpp
#define Edp_hpp

#include <stdio.h>
#include <iostream>
using namespace std;
#include "Matrice.hpp"

Matrice EulerExplicite_Put(double Smax, double K, double Sigma, double r, double M, double N, double T);

Matrice EulerCentré_Put(double Smax, double K, double Sigma, double r, double M, double N, double T);

Matrice EulerImplicite_Put(double Smax, double K, double Sigma, double r, double M, double N, double T);

Matrice Crank_Nicholson(double Smax, double K, double Sigma, double r, double M, double N, double T);




#endif /* Edp_hpp */
