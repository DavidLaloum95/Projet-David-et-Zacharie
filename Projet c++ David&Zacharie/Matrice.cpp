///////////////////////////////////////////////////////////////////////////////////////
/// Matrice.cpp : Construction de la classe matrice et des opérateur spécifiques   /////
///////////////////////////////////////////////////////////////////////////////////////


#include "Matrice.hpp"
#include <iostream>
using namespace std;


// Constructeur d'une Matrice

Matrice:: Matrice(int n, int p){
    _nl=n;
    _nc=p;
    _matrice =  new double * [n];
    for (int i = 0; i < n; i++){
        
        _matrice[i] = new double [p];
    }
    for (int i = 0; i < n; i++){
        for (int j = 0; j < p; j++){
            
            _matrice[i][j]=0;
        }
    }
    
}


// Obtenir le nombdre de colonnes :

double Matrice::dimCol(){
    return _nc;
}


// Obtenir le nombre de ligne :

double Matrice::dimLig(){
    return _nl;
}


// Obtenir le coefficient i,j :

double Matrice:: get(int i, int j){
    return _matrice[i][j];
}

// Donner une valeur au coeficient i,j :

void Matrice::set(int k,int l, double val){
    _matrice[k][l]= val;
}

// Affichage de la matrice :

void Matrice:: print(){
    for (int i = 0; i < _nl; i++){
        for (int j = 0; j <_nc; j++){
            cout << " "<< get(i,j)<<"  ";
        }
        cout <<"\n";
    }
    cout <<"\n";
}


// Programme pour addition de deux matrices de même tailles :

void  Matrice::add(Matrice& a, Matrice& b){
    for (int i = 0; i < _nl; i++){
        for (int j = 0; j <_nc; j++){
            _matrice[i][j]=a.get(i,j) + b.get(i,j);
        }
        
    }
}


// Multiplication de deux matrices de tailles adaptées :

void Matrice::produit(Matrice& a, Matrice& b){
    
    if (a.dimCol() != b.dimLig() ){
        cout << "problème dimension";
        
    }
    else {
        for (int i = 0; i < a.dimLig(); i++){
            for (int j = 0; j < b.dimCol(); j++){
                double res=0;
                
                for (int k = 0; k<b.dimLig();k++){
                    
                    res += a.get(i,k) * b.get(k,j);
                    
                }
                
                set(i,j,res);
            }
        }
    }
    
    _nl=a.dimLig();
    _nc=b.dimCol();
}


// Multiplier une matrice par une constante

void Matrice::prodConst(double a){
    
    for (int i = 0; i < _nl; i++){
        for (int j = 0; j < _nc; j++){
            
            set(i,j,get(i,j)*a);
        }
    }
}

// creer une copie de matrice

void Matrice::equal(Matrice &m){
    for (int i = 0; i < _nl; i++){
        for (int j = 0; j < _nc; j++){
            
            set(i,j,m.get(i,j));
        }
    }
}

// Obtenir la colonne j de la matrice

Matrice Matrice::getCol(int j){
    Matrice res(_nl,1);
    for(int i=0;i<_nl;i++){
        res.set(i,0,get(i,j));
    }
    return res;
    
}

// Remplir la colonne j de ma matrice par la colonne k de la matrice m

void Matrice::setCol(int j, int k ,Matrice& m){
    
    for (int i = 0; i < _nl; i++){
        
        set(i,j,m.get(i,k));
    }
}

// Obtenir la ligne i de la matrice

Matrice Matrice::getLig(int i){
    Matrice res(1,_nc);
    for(int j=0;j<_nc;j++){
        res.set(0,j,get(i,j));
    }
    return res;
    
}

// Echanger la ligne i avec la ligne j

void Matrice::echange_ligne(int i, int j){
    Matrice m(_nl,_nc);
    for (int i = 0; i < _nl; i++){
        for (int j = 0; j < _nc; j++){
            
            m.set(i,j,get(i,j));
        }
    }
    
    
    for(int n=0;n<_nl;n++){
        for (int p=0;p<_nc;p++){
            if( n != i and n!=j ){
                set(n,p,m.get(n,p));
            }
            if(n==i){
                set(i,p,m.get(j,p));
            }
            if(n==j){
                set(j,p,m.get(i,p));
            }
        }
    }
    
}


// Transvection : Li <- Li+mu*Lj

void Matrice::transvection( int i, int j, double mu){
    
    for(int n=0;n<_nl;n++){
        for (int p=0;p<_nc;p++){
            if( n != i){
                set(n,p,get(n,p));
            }
            if(n==i){
                set(i,p,get(i,p)+mu*get(j,p));
            }
            
        }
    }
    
}

// Trouver i tel que le coefficient |Mi,j0| soit maximale

double Matrice::pivot_partiel( int j0){
    
    int imax=j0;
    for (int i= j0+1;i<_nl;i++){
        if(abs(get(i,j0))>abs(get(imax,j0))){
            imax=i;
            
        }
    }
    return imax;
}


// Inverse de la matrice par la méthode du pivot de gauss

void Matrice::inverse(Matrice& m){
    Matrice X(_nl,_nc);
    
    
    for(int i=0; i<X.dimLig();i++){
        
        X.set(i,i,1);
    }
    
    for(int j=0; j<_nl-1 ;j++){
        int imax = pivot_partiel(j);
        echange_ligne(imax,j);
        X.echange_ligne(imax,j);
        
        for(int i=j+1; i<_nl;i++){
            X.transvection(i,j,-(get(i,j)/get(j,j)));
            transvection(i,j,-(get(i,j)/get(j,j)));
            
            
        }
    }
    
    for(int j=_nl-1; j>=0;j--){
        for(int i=j-1;i>=0;i--){
            X.transvection(i,j,-(get(i,j)/get(j,j)));
            transvection(i,j,-(get(i,j)/get(j,j)));
        }
    }
    for(int i=0;i<_nl;i++){
        for(int j=0;j<_nc;j++){
            X.set(i,j,X.get(i,j)/get(i,i));
            set(i,j,get(i,j)/get(i,i));
        }
    }
    
    for (int i = 0; i < _nl; i++){
        for (int j = 0; j < _nc; j++){
            
            set(i,j,X.get(i,j));
        }
    }
}






// Destructeur d'une matrice //

Matrice::~Matrice()
{
    for(int i = 0;i < _nl;++i) {
        
        delete  [] _matrice[i];
        
    }
    
    delete [] _matrice;
    
    
}


//////////////////////////////////////////////////////////////////////////


