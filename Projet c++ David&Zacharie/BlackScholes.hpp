///////////////////////////////////////////////////////////////////////////////////////
////                            BlackScholes.h                                    /////
///////////////////////////////////////////////////////////////////////////////////////


#ifndef BlackScholes_hpp
#define BlackScholes_hpp
#define _USE_MATH_DEFINES

#include <iostream>
#include <cmath>
#include <stdio.h>

double norm_pdf(const double& x);
double norm_cdf(const double& x);
double d_j(const int& j, const double& S, const double& K, const double& r, const double& v, const double& T);
double put_price(const double& S, const double& K, const double& r, const double& v, const double& T);
#endif /* BlackScholes_hpp */
