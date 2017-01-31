///////////////////////////////////////////////////////////////////////////////////////
////                            Fichier Matrice.h                                 /////
///////////////////////////////////////////////////////////////////////////////////////


#ifndef Matrice_hpp
#define Matrice_hpp
#include <iostream>
using namespace std;
#include <stdio.h>
#include <cmath>


class Matrice
{
    
    private :
    int _nl;                                       // Nombre de ligne
    int _nc;                                       // Nombre de Colonne
    double ** _matrice;                            // Tableau dynamique bidimensionnel
    
    
public:
    Matrice(int,int);                              // Constructeur
    ~Matrice();                                    // Destructeur
    double get(int , int );                        // Obtenir le coefficient i,j
    void set(int ,int , double );                  // Donner une valeur au coefficient i,j
    void print();                                  // Affiche la matrice
    void add(Matrice&, Matrice&);                  // Somme de deux matrices
    void produit(Matrice&,Matrice&);               // Produit de deux matrices
    void prodConst(double);                        // Produit d'une matrice par une constante
    
    double dimCol();                               // Nombre de colonne de la matrice
    double dimLig();                               // Nombre de Ligne de la matrice
    
    void equal(Matrice& m);                        // Fait une copie de la matrice
    
    Matrice getCol(int j);                         // Obtenir la colonne j de la matrice
    void setCol(int j,int k, Matrice& m);          // Remplacer la colonne j de la matrice
    // par la colonne i de la matrice m
    
    Matrice getLig(int i);                         // Obtenir la Ligne i de la matrice
    
    
    void echange_ligne( int i, int j);             // Echange la ligne i et la ligne j : Li<->Lj
    void transvection(int i, int j, double mu);    // Li<-Li+mu*Lj
    double pivot_partiel( int j0);                 // Trouver i tel que |Mi,j0| soit maximale
    void inverse(Matrice& m);                      // Inverse une matrice carré
    
};




#endif /* Matrice_hpp */
