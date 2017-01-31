///////////////////////////////////////////////////////////////////////////////////////
////       EDP.cpp: Implémentation de l'algorithme de résolution d'une EDP        /////
///////////////////////////////////////////////////////////////////////////////////////


#include "Edp.hpp"
#include <iostream>
using namespace std;
#include <stdio.h>
#include "Matrice.hpp"


// Premier Algorithme
// Résolution par le schéma d'Euler explicite décentré à droite


Matrice EulerExplicite_Put(double Smax, double K, double Sigma, double r, double M, double N, double T){
    
    
    
    // CREATION DES MATRICES ET DES VECTEURS NECESSAIRES
    // On les remplit de 0 puis, on en changera les valeurs
    
    
    
    Matrice s(M+2,1);      // Maillage de l'esapce
    Matrice temps(N+2,1);  // Maillage du temps
    Matrice A2(M+1,M+1);   // Matrice A2 définie comme sur le rapport
    Matrice A1(M+1,M+1);   // Matrice A1 définie comme sur le rapport
    Matrice P(M+1,N+1);    // Matrice P donc les colonnes sont les P^n
    Matrice Id(M+1,M+1);   // Matrice Identité
    Matrice rId(M+1,M+1);  // Matrice r*Id
    Matrice A(M+1,M+1);    // Matrice A définie comme ID-A avec A définie comme sur le rapport (cf plus loin)
    
    
    const double h= (double)Smax/((double)M+1);  // Pas de l'espace
    const double deltaT= (double)T/((double)N);  // Pas du temps
    
    
    // On remplit le vecteur s en utilisant le maillage d'espace suivant : sj=j*h pour 0<=j<=M+1 :
    
    for (int i=0; i<s.dimLig();i++){
        s.set(i,0,i*h);
    }
    
    // On remplit la matrice temps en utilisant le maillage de temps suivant : tempsj=j*deltaT pour 0<=j<=N+1 :
    
    for (int i=0; i<temps.dimLig();i++){
        temps.set(i,0,i*deltaT);
    }
    
    
    
    
    // CONSTRUCTION DE LA MATRICE Id ET rId :
    
    for(int i=0; i<Id.dimLig();i++){
        
        Id.set(i,i,1);
        rId.set(i,i,r);
        
    }
    
    
    // CONSTRUCTION DE LA MATRICE A1
    
    // On remplit la diagonale de A1 par -rsj :
    
    for(int i=0; i<A1.dimLig();i++){
        
        A1.set(i,i,s.get(i,0)*(-r));
        
    }
    
    // On remplit la sur-diagonale de A1  par rsj:
    
    for (int i = 0; i<A1.dimLig()-1;i++){
        A1.set(i,i+1,s.get(i,0)*(r));
    }
    
    // A1 = A1*1/h :
    
    A1.prodConst((double)1/(2*h));
    
    
    // CONSTRUCTION DE LA MATRICE A2
    
    // on pose aj=(sigma^2*sj^2)/2
    // On remplit la diagonale de A2 par -2aj:
    
    for(int i=0; i<A2.dimLig();i++){
        
        A2.set(i,i,(s.get(i,0)*s.get(i,0)* (Sigma*Sigma)));
        
    }
    
    
    // On remplit la sur-diagonale de A2 par aj:
    
    
    for (int i = 0; i<A2.dimLig()-1;i++){
        A2.set(i,i+1,-(s.get(i,0)*s.get(i,0)* (Sigma*Sigma)/2));
    }
    
    // On remplit la sous-diagonale de A2 par aj:
    
    for (int i = 1; i<A2.dimLig();i++){
        A2.set(i,i-1,-(s.get(i,0)*s.get(i,0)* (Sigma*Sigma)/2));
    }
    
    // A2 = A2 * 1/h^2
    
    A2.prodConst((double) 1/(h*h));
    
    
    
    // CONSTRUCTION DE LA MATRICE A
    
    // On remplit la matrice A = Id-DeltaT*(A2+A1+rId):
    
    A.add(A2,A1);
    A.add(A,rId);
    A.prodConst(-deltaT);
    A.add(A,Id);
    
    
    
    // CONSTRUCTION DE LA MATRICE P
    
    // On remplit la Matrice P en utilisant les conditions initiales :
    
    for (int i=0; i<P.dimLig();i++){
        if (K-s.get(i,0) >= 0){
            P.set(i,0,(K-s.get(i,0)));
        }
        else {
            P.set(i,0,0);
            
        }}
    
    
    // MATRICE TEMPORAIRE TOT
    
    Matrice TOT(M+1,M+1);
    
    // On remplit la matrice P colonne par colonne
    // en utilisant : Pn+1=A*Pn
    
    for(int n=0; n<P.dimCol()-1;n++){
        Matrice res(P.dimLig(),1);
        for(int i=0;i<P.dimLig();i++){
            res.set(i,0,P.get(i,n));
        }
        
        TOT.produit(A,res);
        
        for (int i = 0; i < P.dimLig(); i++){
            
            P.set(i,n+1,TOT.get(i,0));
        }
    }
    
    // On retourne le prix du put pour les différents sj
    
    
    return P.getCol(P.dimCol()-1);
    
}





////////////////////////////////////////////////////////////////////////////////////


// Deuxième Algorithme
// Résolution par le schéma d'Euler Explicite centré.


Matrice EulerCentré_Put(double Smax, double K, double Sigma, double r, double M, double N, double T){
    
    // CREATION DES MATRICES ET DES VECTEURS NECESSAIRES
    // On les remplit de 0 puis, on en changera les valeurs
    
    
    
    Matrice s(M+2,1);      // Maillage de l'esapce
    Matrice temps(N+2,1);  // Maillage du temps
    Matrice A2(M+1,M+1);   // Matrice A2 définie comme sur le rapport
    Matrice A1(M+1,M+1);   // Matrice A1 définie comme sur le rapport
    Matrice P(M+1,N+1);    // Matrice P donc les colonnes sont les P^n
    Matrice Id(M+1,M+1);   // Matrice Identité
    Matrice rId(M+1,M+1);  // Matrice r*Id
    Matrice A(M+1,M+1);    // Matrice A définie comme ID-A avec A définie comme sur le rapport (cf plus loin)
    
    
    const double h= (double)Smax/((double)M+1);  // Pas de l'espace
    const double deltaT= (double)T/((double)N);  // Pas du temps
    
    
    // On remplit le vecteur s en utilisant le maillage d'espace suivant : sj=j*h pour 0<=j<=M+1 :
    
    for (int i=0; i<s.dimLig();i++){
        s.set(i,0,i*h);
    }
    
    // On remplit la matrice temps en utilisant le maillage de temps suivant : tempsj=j*deltaT pour 0<=j<=N+1 :
    
    for (int i=0; i<temps.dimLig();i++){
        temps.set(i,0,i*deltaT);
    }
    
    
    // CONSTRUCTION DE LA MATRICE Id ET rId
    
    // On remplit la matrice Id et rId :
    
    for(int i=0; i<Id.dimLig();i++){
        
        Id.set(i,i,1);
        rId.set(i,i,r);
        
    }
    
    
    
    // CONSTRUCTION DE LA MATRICE A1
    
    // On remplit la SOUS diagonale de A1 (a la différence du schéma décentré à droite):
    
    for(int i=1; i<A1.dimLig();i++){
        
        A1.set(i,i-1,s.get(i,0)*(-r));
        
    }
    
    // On remplit la sur-diagonale de A1 :
    
    for (int i = 0; i<A1.dimLig()-1;i++){
        A1.set(i,i+1,s.get(i,0)*(r));
    }
    
    // A1 = A1*1/2h :
    
    A1.prodConst((double)1/(2*h));
    
    
    
    // CONSTRUCTION DE LA MATRICE A2
    
    // On remplit la diagonale de A2 :
    
    for(int i=0; i<A2.dimLig();i++){
        
        A2.set(i,i,(s.get(i,0)*s.get(i,0)* (Sigma*Sigma)));
        
    }
    
    
    // On remplit la sur-diagonale de A2 :
    
    
    for (int i = 0; i<A2.dimLig()-1;i++){
        A2.set(i,i+1,-(s.get(i,0)*s.get(i,0)* (Sigma*Sigma)/2));
    }
    
    // On remplit la sous-diagonale de A2 :
    
    for (int i = 1; i<A2.dimLig();i++){
        A2.set(i,i-1,-(s.get(i,0)*s.get(i,0)* (Sigma*Sigma)/2));
    }
    
    // A2 = A2 * 1/h^2
    
    A2.prodConst((double) 1/(h*h));
    
    
    
    // CONSTRUCTION DE LA MATRICE A
    
    // On remplit la matrice A :
    
    A.add(A2,A1);
    A.add(A,rId);
    A.prodConst(-deltaT);
    A.add(A,Id);
    
    
    
    // CONSTRUCTION DE LA MATRICE P
    
    // On remplit la Matrice P en utilisant les conditions initiales :
    
    for (int i=0; i<P.dimLig();i++){
        if (K-s.get(i,0) >= 0){
            P.set(i,0,(K-s.get(i,0)));
        }
        else {
            P.set(i,0,0);
            
        }
    }
    
    // MATRICE TEMPORAIRE TOT
    
    Matrice TOT(M+1,M+1);
    
    // On remplit la matrice P colonne par colonne
    // en utilisant : Pn+1=A*Pn
    
    for(int n=0; n<P.dimCol()-1;n++){
        Matrice res(P.dimLig(),1);
        for(int i=0;i<P.dimLig();i++){
            res.set(i,0,P.get(i,n));
        }
        
        TOT.produit(A,res);
        
        for (int i = 0; i < P.dimLig(); i++){
            
            P.set(i,n+1,TOT.get(i,0));
        }
    }
    
    // On retourne le prix du put pour les différents sj
    
    return P.getCol(P.dimCol()-1);
    
}


////////////////////////////////////////////////////////////////////////////////////


// Troisième Algorithme
// Résolution par le schéma d'Euler Implicite.




Matrice EulerImplicite_Put(double Smax, double K, double Sigma, double r, double M, double N, double T){
    
    // CREATION DES MATRICES ET DES VECTEURS NECESSAIRES
    // On les remplit de 0 puis, on en changera les valeurs
    
    
    
    Matrice s(M+2,1);      // Maillage de l'esapce
    Matrice temps(N+2,1);  // Maillage du temps
    Matrice A2(M+1,M+1);   // Matrice A2 définie comme sur le rapport
    Matrice A1(M+1,M+1);   // Matrice A1 définie comme sur le rapport
    Matrice P(M+1,N);      // Matrice P donc les colonnes sont les P^n
    Matrice Id(M+1,M+1);   // Matrice Identité
    Matrice rId(M+1,M+1);  // Matrice r*Id
    Matrice A(M+1,M+1);    // Matrice A définie comme sur le rapport
    
    
    const double h= (double)Smax/((double)M+1);  // Pas de l'espace
    const double deltaT= (double)T/((double)N);  // Pas du temps
    
    
    // On remplit le vecteur s en utilisant le maillage d'espace suivant : sj=j*h pour 0<=j<=M+1 :
    
    for (int i=0; i<s.dimLig();i++){
        s.set(i,0,i*h);
    }
    
    // On remplit la matrice temps en utilisant le maillage de temps suivant : tempsj=j*deltaT pour 0<=j<=N+1 :
    
    for (int i=0; i<temps.dimLig();i++){
        temps.set(i,0,i*deltaT);
    }
    
    
    // CONSTRUCTION DE LA MATRICE Id ET rId
    
    // On remplit la matrice Id et rId :
    
    for(int i=0; i<Id.dimLig();i++){
        
        Id.set(i,i,1);
        rId.set(i,i,r);
        
    }
    
    
    // CONSTRUCTION DE LA MATRICE A1
    
    // On remplit la SOUS diagonale de A1 :
    for(int i=1; i<A1.dimLig();i++){
        
        A1.set(i,i-1,s.get(i,0)*(-r));
        
    }
    
    // On remplit la sur-diagonale de A1 :
    
    for (int i = 0; i<A1.dimLig()-1;i++){
        A1.set(i,i+1,s.get(i,0)*(r));
    }
    
    // A1 = A1*1/2h :
    
    A1.prodConst((double)1/(2*h));
    
    
    
    // CONSTRUCTION DE LA MATRICE A2
    
    // On remplit la diagonale de A2 :
    
    for(int i=0; i<A2.dimLig();i++){
        
        A2.set(i,i,(s.get(i,0)*s.get(i,0)* (Sigma*Sigma)));
        
    }
    
    
    // On remplit la sur-diagonale de A2 :
    
    
    for (int i = 0; i<A2.dimLig()-1;i++){
        A2.set(i,i+1,-(s.get(i,0)*s.get(i,0)* (Sigma*Sigma)/2));
    }
    
    // On remplit la sous-diagonale de A2 :
    
    for (int i = 1; i<A2.dimLig();i++){
        A2.set(i,i-1,-(s.get(i,0)*s.get(i,0)* (Sigma*Sigma)/2));
    }
    
    // A2 = A2 * 1/h^2
    
    A2.prodConst((double) 1/(h*h));
    
    
    
    // CONSTRUCTION DE LA MATRICE A
    
    // On remplit la matrice A = Id + deltaT*( A2+A1+rId):
    
    A.add(A2,A1);
    A.add(A,rId);
    A.prodConst(deltaT);
    A.add(A,Id);
    
    // On inverse A: A^-1
    
    A.inverse(A);
    
    
    
    
    // CONSTRUCTION DE LA MATRICE P
    
    // On remplit la Matrice P en utilisant les conditions initiales :
    
    for (int i=0; i<P.dimLig();i++){
        if (K-s.get(i,0) >= 0){
            P.set(i,0,(K-s.get(i,0)));
        }
        else {
            P.set(i,0,0);
            
        }
        
    }
    
    // MATRICE TEMPORAIRE TOT
    
    Matrice TOT(M+1,M+1);
    
    
    // On remplit la matrice P colonne par colonne
    // en utilisant : Pn+1=(A^-1)*Pn
    
    
    for(int n=0; n<P.dimCol()-1;n++){
        Matrice res(P.dimLig(),1);
        for(int i=0;i<P.dimLig();i++){
            res.set(i,0,P.get(i,n));
        }
        
        TOT.produit(A,res);
        
        for (int i = 0; i < P.dimLig(); i++){
            
            P.set(i,n+1,TOT.get(i,0));
        }
        
    }
    
    
    // On retourne le prix du put pour les différents sj
    
    return P.getCol(P.dimCol()-1);
    
}


////////////////////////////////////////////////////////////////////////////////////


//  Quatrième Algorithme
// Résolution par le schéma de Crank Nicholson.





Matrice Crank_Nicholson(double Smax, double K, double Sigma, double r, double M, double N, double T){
    
    // CREATION DES MATRICES ET DES VECTEURS NECESSAIRES
    // On les remplit de 0 puis, on en changera les valeurs
    
    
    
    Matrice s(M+2,1);      // Maillage de l'esapce
    Matrice temps(N+2,1);  // Maillage du temps
    Matrice A2I(M+1,M+1);  // Matrice A2 définie comme pour Euler Implicite
    Matrice A2(M+1,M+1);   // Matrice A2 définie comme pour Euler Explicite
    Matrice A1(M+1,M+1);   // Matrice A1 définie comme sur le rapport
    Matrice P(M+1,N);      // Matrice P donc les colonnes sont les P^n
    Matrice Id(M+1,M+1);   // Matrice Identité
    Matrice rId(M+1,M+1);  // Matrice r*Id
    Matrice A(M+1,M+1);    // Matrice A définie comme sur le rapport
    
    
    const double h= (double)Smax/((double)M+1);  // Pas de l'espace
    const double deltaT= (double)T/((double)N);  // Pas du temps
    
    
    // On remplit le vecteur s en utilisant le maillage d'espace suivant : sj=j*h pour 0<=j<=M+1 :
    
    for (int i=0; i<s.dimLig();i++){
        s.set(i,0,i*h);
    }
    
    // On remplit la matrice temps en utilisant le maillage de temps suivant : tempsj=j*deltaT pour 0<=j<=N+1 :
    
    for (int i=0; i<temps.dimLig();i++){
        temps.set(i,0,i*deltaT);
    }
    
    
    // CONSTRUCTION DE LA MATRICE Id ET rId
    
    // On remplit la matrice Id et rId :
    
    for(int i=0; i<Id.dimLig();i++){
        
        Id.set(i,i,1);
        rId.set(i,i,r);
        
    }
    
    
    // CONSTRUCTION DE LA MATRICE A1
    
    // On remplit la SOUS diagonale de A1 :
    for(int i=1; i<A1.dimLig();i++){
        
        A1.set(i,i-1,s.get(i,0)*(-r));
        
    }
    
    // On remplit la sur-diagonale de A1 :
    
    for (int i = 0; i<A1.dimLig()-1;i++){
        A1.set(i,i+1,s.get(i,0)*(r));
    }
    
    // A1 = A1*1/2h :
    
    A1.prodConst((double)1/(2*h));
    
    
    
    // CONSTRUCTION DE LA MATRICE A2I
    
    // On remplit la diagonale de A2I :
    
    for(int i=0; i<A2I.dimLig();i++){
        
        A2I.set(i,i,(s.get(i,0)*s.get(i,0)* (Sigma*Sigma)));
        
    }
    
    
    // On remplit la sur-diagonale de A2I :
    
    
    for (int i = 0; i<A2I.dimLig()-1;i++){
        A2I.set(i,i+1,-(s.get(i,0)*s.get(i,0)* (Sigma*Sigma)/2));
    }
    
    // On remplit la sous-diagonale de A2I :
    
    for (int i = 1; i<A2I.dimLig();i++){
        A2I.set(i,i-1,-(s.get(i,0)*s.get(i,0)* (Sigma*Sigma)/2));
    }
    
    // A2I = A2 * 1/h^2
    
    A2I.prodConst((double) 1/(h*h));
    
    
    // CONSTRUCTION DE LA MATRICE A2
    
    // On remplit la diagonale de A2 :
    
    for(int i=0; i<A2.dimLig();i++){
        
        A2.set(i,i,(s.get(i,0)*s.get(i,0)* (Sigma*Sigma)));
        
    }
    
    
    // On remplit la sur-diagonale de A2 :
    
    
    for (int i = 0; i<A2.dimLig()-1;i++){
        A2.set(i,i+1,-(s.get(i,0)*s.get(i,0)* (Sigma*Sigma)/2));
    }
    
    // On remplit la sous-diagonale de A2 :
    
    for (int i = 1; i<A2.dimLig();i++){
        A2.set(i,i-1,-(s.get(i,0)*s.get(i,0)* (Sigma*Sigma)/2));
    }
    
    // A2 = A2 * 1/h^2
    
    A2.prodConst((double) 1/(h*h));
    
    
    
    
    
    // CONSTRUCTION DE LA MATRICE A
    
    // On remplit la matrice A = Id + deltaT*(A2/2+A1+rId):
    A2I.prodConst((double) 1 /(double)2);
    A.add(A2I,A1);
    A.add(A,rId);
    A.prodConst(deltaT);
    A.add(A,Id);
    
    // On inverse A: A^-1
    
    A.inverse(A);
    
    // On multiplie A par Id - DeltatT(A1/2))
    A1.prodConst(-deltaT/2);
    A1.add(A1,Id);
    A.produit(A,A1);
    
    
    
    
    // CONSTRUCTION DE LA MATRICE P
    
    // On remplit la Matrice P en utilisant les conditions initiales :
    
    for (int i=0; i<P.dimLig();i++){
        if (K-s.get(i,0) >= 0){
            P.set(i,0,(K-s.get(i,0)));
        }
        else {
            P.set(i,0,0);
            
        }
        
    }
    
    // MATRICE TEMPORAIRE TOT
    
    Matrice TOT(M+1,M+1);
    
    
    // On remplit la matrice P colonne par colonne
    // en utilisant : Pn+1=(A^-1)*Pn
    
    
    for(int n=0; n<P.dimCol()-1;n++){
        Matrice res(P.dimLig(),1);
        for(int i=0;i<P.dimLig();i++){
            res.set(i,0,P.get(i,n));
        }
        
        TOT.produit(A,res);
        
        for (int i = 0; i < P.dimLig(); i++){
            
            P.set(i,n+1,TOT.get(i,0));
        }
        
    }
    
    
    // On retourne le prix du put pour les différents sj
    
    return P.getCol(P.dimCol()-1);
}



