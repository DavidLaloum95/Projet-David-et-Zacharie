///////////////////////////////////////////////////////////////////////////////////////
////  Projet C++ = Pricing de produit dérivé par EDP : Modèle de Black & Scholes  /////
///////////////////////////////////////////////////////////////////////////////////////

#define _USE_MATH_DEFINES

#include <cmath>
#include <iostream>
using namespace std;
#include "Matrice.hpp"
#include "BlackScholes.hpp"
#include "Edp.hpp"


int main(int argc, const char * argv[]) {
    
    
    // CONDITIONS INITIALES :
    
    
    const int M=19;                              // Maillage de L'espace
    const int N=20;                              // Maillage du temps
    const int Smax = 200;                        // Valeur de l'action est bornée
    const int K=100;                             // Strike ou prix d'exercice de l'option
    const int T=1;                               // Maturité de l'option
    const double r=0.1;                          // Taux d'interêt
    const double Sigma = 0.2;                    // Volatilité de l'action
    const double h= (double)Smax/((double)M+1);  // Pas de l'espace
    
    
    
    
    cout << "Valeur du Put par le schéma d'Euler explicite :"<<endl;
    EulerExplicite_Put(Smax, K, Sigma, r,  M, N, T).print();
    cout << "Valeur du Put par le schéma d'Euler centré :"<<endl;
    EulerCentré_Put(Smax, K, Sigma, r,  M, N, T).print();
    cout << "Valeur du Put par le schéma d'Euler implicite :"<<endl;
    EulerImplicite_Put(Smax, K, Sigma, r,  M, N, T).print();
    cout << "Valeur du Put par le schéma de Crank Nicholson :"<<endl;
    Crank_Nicholson(Smax, K, Sigma, r,  M, N, T).print();
    
    
    // Affichage des valeurs théoriques
    
    cout << "Valeur du Put théorique : "<<endl;
    
    for(int i=0; i<M+1;i++){
        double put = put_price(h*i, K, r, Sigma, T);
        cout <<put<<endl;
        
    }
    
    
    
    
    
    
    
    
    return EXIT_SUCCESS;
}











