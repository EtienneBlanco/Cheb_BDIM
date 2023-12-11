/*======================================================\
|                                                       |
|          Functions used in this resolution            |
|	                                                    |
\======================================================*/


// On this file are gathered all the definitions of the mathematical functions used.

#ifndef CHEB_FUNCTIONS_H
#define CHEB_FUNCTIONS_H

#include <math.h>

//////////////////////////////////////////////////
//prototype of the used functions               //
double cheb_T(int n, double x);                 //
double cheb_int(int i);                         //
double y_(double x, double a, double b);        //
double y_inv(double x, double a, double b);     //
double y_log(double x, double a, double b);     //
double y_log_inv(double x, double a, double b); //
double node(int N, int i);                      //
double K_gg(double x);                          //
double K_simpl(double x);                       //
double K_gq(double x);                          //
double K_qq(double x);                          //
double K_qg(double x);                          //
double D_t0_e(double x, double eps);            //
double D_t0_a(double x, double tau);            //
//////////////////////////////////////////////////


//return the value of the Chebyshev polynomial T_n on x
double cheb_T(int n, double x)
{
    double T;
    if (x<-1)
    {
        cout << "\nError cheb_T cannot be evaluated on " << x;
        return 0;
    }
    else if (x>1)
    {
        cout << "\nError cheb_T cannot be evaluated on " << x;
        return 0;
    }
    T = cos(n*acos(x));

    return T;
}

//Return the integral of T_i (over [-1,1])
double cheb_int(int i)
{
    double y;
    if (i==1)
        //y = .5;
        y = 0;
    else
        //y = (i*sin(i*M_PI/2)-1)/double(i*i-1);
        y = (pow(-1., double(i))+1)/double(1-i*i);
    return y;
}

//roots of the Chebyshev polynomials
double node(int N, int i){
    return cos(M_PI/N*(i+.5));}

//linear bijection going from [a,b] to [-1,1]
double y_(double x, double a, double b){
    return (2*x-b-a)/(b-a);}

//inverse of y_
double y_inv(double x, double a, double b){
    return 0.5*((b-a)*x+b+a);}

//logarithmic bijection going from [a,b] to [-1,1]
double y_log(double x, double a, double b){
    return 1+2*log(x/b)/log(b/a);}

//inverse of y_
double y_log_inv(double x, double a, double b){
    return b*pow(b/a, .5*(x-1));}


//Kernel of the evolution equation for the pure gluon case
double K_gonly(double x){
    return pow(1-x+x*x,2.5)/pow(x*(1-x),1.5);}

//Simplified kernel of the evolution equation (in the pure gluons case)
double K_simpl(double x){
    return 1./pow(x*(1-x),1.5);}

//Kernel of the evolution equation (gluon to gluon part)
double K_gg(double x){
    return pow(C_A,1.5)*pow(1-x+x*x,2.5)/pow(x*(1-x),1.5);}

//Kernel of the evolution equation (gluon to quark part)
double K_gq(double x){
    return C_F/2.*(1+pow(1-x,2))/x*sqrt(((1-x)*C_A+x*x*C_F)/(x*(1-x)));}

//Kernel of the evolution equation (quark to quark part)
double K_qq(double x){
    return C_F/2.*(1+x*x)/(1-x)*sqrt((x*C_A+pow(1-x,2)*C_F)/(x*(1-x)));}

//Kernel of the evolution equation (quark to gluon part)
double K_qg(double x){
    return Nf*T_F*(x*x+pow(1-x,2))*sqrt((C_F-x*(1-x)*C_A)/(x*(1-x)));}


//Possible initial solution at t_0 : a Gaussian picked in 1 to mimic a delta function
double D_t0_e(double x, double eps){
    return sqrt(2./M_PI)/eps*exp(-pow((1-x)/eps,2.)/2.);}

//Possible initial solution at t_0 : the analytical solution of the simplified BDIM equations
double D_t0_a(double x, double tau){
    return tau/sqrt(x*pow(1-x,3.))*exp(-M_PI*tau*tau/(1-x));}


#endif // CHEB_FUNCTIONS_H
