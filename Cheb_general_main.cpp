/*=============================================================================\
|                                                                              |
|  Resolution of the integrated BDIM equations (accounting or not for quarks)  |
|	                                                                           |
\=============================================================================*/



//The resolution of the evolution equations is done by expanding the distribution functions in
//Chebyshev polynomials and the time differential equation is solved by a simple Euler method.
//
//The result is given as a grid in x and t. The point in x are the nodes of the Chebyshev polynomials.


#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
using namespace std;

const long double M_PI = 3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679;
const int Nc = 3; //Number of colors
const int Nf = 3; //Number of active quark flavors
const double C_A = Nc;
const double C_F = (Nc*Nc-1)/2./Nc;
const double T_F = .5;

//Associated files
#include "Cheb_functions.h" //file containing functions related to Chebyshev polynomials and also to the kernel and the starting FF
#include "Cheb_parameters.h" //file containing the parameters classS
#include "Cheb_method_g_S_NS.h" //file containing functions related to the resolution method for the system of equations on g, S and NS FF
#include "Cheb_method_g.h" //file containing functions related to the resolution method for the pure gluon case

void parameters_from_file(int *N, double *t0, double *dt, int *n_t, int *n_t_w, double *eps);
void cheb_write_grid_header(int N, double eps, ofstream &grid);
int main_g(parameters param);
int main_g_S_NS(parameters param);

int main(void)
{
    //Parameters reading
    ////////////////////
    parameters param;
    param.from_file("Parameters");

    int evo = param.evo;

    if (evo == 0)
        main_g(param);
    else if (evo == 1)
        main_g_S_NS(param);
    else
    {
        cout << "\n# /!\\ Wrong parameter 'evo' : " << evo << endl;
        cout << "# Should be set to 0 or 1 : " << endl;
        cout << "# Program stopped" << endl;
        return 1;
    }
    return 0;
}

int main_g(parameters param)
{
    //Parameters initialization
    ///////////////////////////
    int Nx = param.Nx;
    int n_t = param.Nt;
    int n_t_w = param.Nt_w;
    double t0 = param.t0;
    double dt = param.dt;
    double eps = param.eps;

    //Memory allocation
    ///////////////////
    cout << "\n-Memory allocation";
    double *D = new double[Nx]; //Fragmentation function D_i,t=D(x_i,t_0+dt*t)
    double **T = new double*[Nx]; //Chebyshev polynomial i evaluated on node j
    double *X = new double[Nx]; //Images of the roots of T_N
    double **S = new double*[Nx]; //Matrix S

    for (int i = 0; i < Nx; ++i)
    {
        T[i] = new double[Nx];
        S[i] = new double[Nx];
    }

    cout << "\n*Memory allocated*";

    //Initialization (time independent matrix calculation)
    ////////////////
    cheb_initilization_g(X, T, S, param);
    cout << "\n*Initialization done*";

    //Grid files preparation
    ////////////////////////
    ofstream D_grid;
    string filename = param.gridname;

    param.write();
    D_grid.open(filename+"gonly", ios::out);
    cout << "\n\nGrid file created";

    //Euler method
    //////////////
    cout << "\n\nEuler method launched";
    cheb_Euler_g(X, T, D, S, D_grid, param);
    cout << "\n*Equation solved*";

    //Memory liberation
    ///////////////////
    D_grid.close();
    cout << "\n\n*Grid file closed*";
    for (int i = 0; i < Nx; ++i)
        delete[] S[i];
    delete[] S;
    delete[] D;
    delete[] X;

    cout << "\n\n*Freed memory*\n";


    return 0;
}

int main_g_S_NS(parameters param)
{
    //Parameters initialization
    ///////////////////////////
    int Nx = param.Nx;
    int n_t = param.Nt;
    int n_t_w = param.Nt_w;
    double t0 = param.t0;
    double dt = param.dt;
    double eps = param.eps;

    if (n_t_w>n_t)
        n_t_w = n_t;

    //Memory allocation
    ///////////////////
    cout << "\n-Memory allocation";
    double *D_g = new double[Nx]; //Gluon fragmentation function D_i,t=D(x_i,t_0+dt*t)
    double *D_S = new double[Nx]; //Singlet fragmentation function D_i,t=D(x_i,t_0+dt*t)
    double *D_NS = new double[Nx]; //Non-singlet fragmentation function D_i,t=D(x_i,t_0+dt*t)
    double **T = new double*[Nx]; //Chebyshev polynomial i evaluated on node j
    double *X = new double[Nx]; //Images of the roots of T_N
    double **S_gg = new double*[Nx]; //Some sum
    double **S_gS = new double*[Nx]; //Some sum
    double **S_SS = new double*[Nx]; //Some sum
    double **S_Sg = new double*[Nx]; //Some sum
    double **S_NSNS = new double*[Nx]; //Some sum

    for (int i = 0; i < Nx; ++i)
    {
        T[i] = new double[Nx];
        S_gg[i] = new double[Nx];
        S_gS[i] = new double[Nx];
        S_SS[i] = new double[Nx];
        S_Sg[i] = new double[Nx];
        S_NSNS[i] = new double[Nx];
    }
    cout << "\n*Memory allocated*";

    //Initialization (time independent matrix calculation)
    ////////////////
    cheb_initilization(X, T, S_gg, S_gS, S_SS, S_Sg, S_NSNS, param);
    cout << "\n*Initialization done*";

    //Grid files preparation
    ////////////////////////
    ofstream g_grid, S_grid, NS_grid;
    string filename = param.gridname;

    param.write();
    g_grid.open (filename+"g", ios::out);
    S_grid.open (filename+"S", ios::out);
    NS_grid.open (filename+"NS", ios::out);
    cout << "\n\nGrid files created";

    //Euler method
    //////////////
    cout << "\n\nEuler method launched";
    cheb_Euler(X, T, D_g, D_S, D_NS, S_gg, S_gS, S_SS, S_Sg, S_NSNS, g_grid, S_grid, NS_grid, param);
    cout << "\n*Equation solved*";

    //Memory liberation
    ///////////////////
    g_grid.close();
    S_grid.close();
    NS_grid.close();
    cout << "\n\n*Grid files closed*";

    for (int i = 0; i < Nx; ++i)
    {
        delete[] S_gg[i];
        delete[] S_gS[i];
        delete[] S_SS[i];
        delete[] S_Sg[i];
        delete[] S_NSNS[i];
    }
    delete[] S_gg;
    delete[] S_gS;
    delete[] S_SS;
    delete[] S_Sg;
    delete[] S_NSNS;
    delete[] D_g;
    delete[] D_S;
    delete[] D_NS;
    delete[] X;

    cout << "\n\n*Freed memory*";


    return 0;
}
