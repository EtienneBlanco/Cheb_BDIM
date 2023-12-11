/*=============================================================================\
|                                                                              |
|  Resolution of the integrated BDIM equations (accounting or not for quarks)  |
|	                                                                           |
\=============================================================================*/



//The resolution of the evolution equations is done by expanding the distribution functions in
//Chebyshev polynomials and the time differential equation is solved by a simple Euler method.
//
//The result is given as a grid in x and t. The points in x are the nodes of the Chebyshev polynomials.


#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include  <string.h>
using namespace std;

const long double M_PI = 3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679;
const int Nc = 3; //Number of colors
const int Nf = 3; //Number of active quark flavors
//Casimir operators
const double C_A = Nc;
const double C_F = (Nc*Nc-1)/2./Nc;
const double T_F = .5;

//Associated files
#include "Cheb_functions.h" //file containing functions related to Chebyshev polynomials and also to the kernel and the starting FF
#include "Cheb_parameters.h" //file containing the parameters classS
#include "Cheb_method_g_S_NS.h" //file containing functions related to the resolution method for the system of equations on g, S and NS FF
#include "Cheb_method_g.h" //file containing functions related to the resolution method for the pure gluon case

int program_init(string parametersfile); //reading of the parameters (with the choice of the BDIM equation to solve)
int program_g(parameters param);
int program_g_S_NS(parameters param);
int help(void); //help display

//argument management
int main(int argc, char** argv)
{
    if (argc == 1)
        //launch the program with default value for the name of the parameters file
        program_init("Parameters");
    else if (argc == 2)
    {
        //if (argv[1] == "help")
        if (strcmp(argv[1], "help")==0)
            //display help
            help();
        else
        {
            //launch the program with specified name of the parameters file
            return program_init(argv[1]);
        }
    }
    else
    {
        //wrong input display
        cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
        cout << "!!  /!\\ Wrong argument(s)  !!" << endl;
        cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;

        cout << "\nThis program should be launched with either :" << endl;
        cout << "-no argument (parameters are initialized from the file 'Parameters') :" << endl;
        cout << "\t" << argv[0] << endl;
        cout << "-the name of the parameters file to be use in argument :" << endl;
        cout << "\t" << argv[0] << " 'parameterfile'" << endl;
        cout << "-help argument to display help :" << endl;
        cout << "\t" << argv[0] << " help" << endl;
    }
    return 1;
}

//help
int help(void)
{
    cout << "\n#////////////#" << endl;
    cout << "#//  Help  //#" << endl;
    cout << "#/////////////////////////////////////////////#" << endl;

    cout << "# This program calculates solution of the integrated BDIM equations either accounting for quarks in the cascade or not." << endl;
    cout << "# The fragmentation functions are evolved from an initial distribution, during a time tau." << endl;
    cout << "# The solutions are then given in terms of grids on the nodes of the Nth Chebyshev polynomial in x and using a regular spacing in time." << endl;
    cout << "#" << endl;
    cout << "# More details on the parameters are given in the file 'Parameters' and in the manual." << endl;
    cout << "#" << endl;
    cout << "# The results are given in either 2 files for gluon dominated cascades." << endl;
    cout << "# These files are the gluon fragmentation function D_g and the parameter file used to produce it, with suffixes as follow : " << endl;
    cout << "# -'gonly' for D_g as a 3 column grid (tau, x, D_g(x,tau))" << endl;
    cout << "# -'_parameters' for the parameter file (written as parameters appear in the console when they are read)" << endl;
    cout << "#" << endl;
    cout << "# In the case accounting for quarks in the cascades, 4 files are created." << endl;
    cout << "# These files are the gluon/color singlet/color non singlet fragmentation function D_g/D_S/D_NS and the parameter file used to produce them, with suffixes as follow : " << endl;
    cout << "# -'g' for D_g as a 3 column grid (tau, x, D_g(x,tau))" << endl;
    cout << "# -'S' for D_S as a 3 column grid (tau, x, D_S(x,tau))" << endl;
    cout << "# -'NS' for D_NS as a 3 column grid (tau, x, D_NS(x,tau))" << endl;
    cout << "# -'_parameters' for the parameter file (written as parameters appear in the console when they are read)" << endl;
    cout << "#" << endl;
    cout << "# The program should be launched either with the name of the parameters file in argument or without argument." << endl;
    cout << "# In this case, parameters are read from the file 'Parameters' if present or set to default values." << endl;

    cout << "#/////////////////////////////////////////////#\n" << endl;

    return 0;
}

int program_init(string parametersfile)
{
    //Parameters reading
    ////////////////////
    parameters param;
    param.from_file(parametersfile);

    int evo = param.evo; // type of cascade

    if ((evo == 0)||(evo == 2))
        program_g(param);
    else if (evo == 1)
        program_g_S_NS(param);
    else
    {
        cout << "\n# /!\\ Wrong parameter 'evo' : " << evo << endl;
        cout << "# Should be set to 0 or 1" << endl;
        cout << "# Program stopped" << endl;
        return 1;
    }
    return 0;
}

int program_g(parameters param)
{
    //Parameters initialization
    ///////////////////////////
    int Nx = param.Nx;
    int n_t = param.Nt;
    int n_t_w = param.Nt_w;

    if (n_t_w>n_t)
        n_t_w = n_t;

    //Memory allocation
    ///////////////////
    cout << "\n-Memory allocation" << endl;
    double *D = new double[Nx]; //Fragmentation function D_i,t=D(x_i,t_0+dt*t)
    double **T = new double*[Nx]; //Chebyshev polynomial i evaluated on node j
    double *X = new double[Nx]; //Images of the roots of T_N
    double **S = new double*[Nx]; //Matrix S

    for (int i = 0; i < Nx; ++i)
    {
        T[i] = new double[Nx];
        S[i] = new double[Nx];
    }

    cout << "*Memory allocated*" << endl;

    //Initialization (time independent matrix calculation)
    ////////////////
    cheb_initilization_g(X, T, S, param);
    cout << "*Initialization done*" << endl;

    //Grid files preparation
    ////////////////////////
    ofstream D_grid;
    string filename = param.gridname;

    param.write();
    D_grid.open(filename+"gonly", ios::out);
    cout << "\nGrid file created" << endl;

    //Euler method
    //////////////
    cout << "\nEuler method launched" << endl;
    cheb_Euler_g(X, T, D, S, D_grid, param);
    cout << "*Equation solved*" << endl;

    //Memory liberation
    ///////////////////
    D_grid.close();
    cout << "\n*Grid file closed*" << endl;
    for (int i = 0; i < Nx; ++i)
        delete[] S[i];
    delete[] S;
    delete[] D;
    delete[] X;

    cout << "\n*Freed memory*" << endl;


    return 0;
}

int program_g_S_NS(parameters param)
{
    //Parameters initialization
    ///////////////////////////
    int Nx = param.Nx;
    int n_t = param.Nt;
    int n_t_w = param.Nt_w;

    if (n_t_w>n_t)
        n_t_w = n_t;

    //Memory allocation
    ///////////////////
    cout << "\n-Memory allocation" << endl;
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
    cout << "\n*Memory allocated*" << endl;

    //Initialization (time independent matrix calculation)
    ////////////////
    cheb_initilization(X, T, S_gg, S_gS, S_SS, S_Sg, S_NSNS, param);
    cout << "\n*Initialization done*" << endl;

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
