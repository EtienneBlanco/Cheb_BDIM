/*======================================================\
|                                                       |
|   Resolution method based on Chebyshev polynomials    |
|	                                                    |
\======================================================*/


// On this file are gathered the functions directly related to the resolution method.

#ifndef CHEB_METHOD_G_S_NS_H
#define CHEB_METHOD_G_S_NS_H

#include <math.h>

////////////////////////////////////////
//prototype of the used functions     //
////////////////////////////////////////
void cheb_init_S_chebint(double *&X, double **&T, double **&S_gg, double **&S_gS, double **&S_SS, double **&S_Sg, double **&S_NSNS, int N);
void cheb_init_S_chebint_log(double *&X, double **&T, double **&S_gg, double **&S_gS, double **&S_SS, double **&S_Sg, double **&S_NSNS, int N, double eps);

void write_D(double *&X, double **&T, double *&D_g, double *&D_S, double *&D_NS, double t, ofstream &g_grid, ofstream &S_grid, ofstream &NS_grid, parameters param);
void cheb_Euler(double *&X, double **&T, double *&D_g, double *&D_S, double *&D_NS, double **&S_gg, double **&S_gS, double **&S_SS, double **&S_Sg, double **&S_NSNS, ofstream &g_grid, ofstream &S_grid, ofstream &NS_grid, parameters param);

double normalize_D(double *&X, double **&T, double D_g[], double D_S[], double D_NS[], parameters param);
double normalize_D(double *&X, double **&T, double D[], parameters param);
////////////////////////////////////////


//Function that initialize useful elements
void cheb_initilization(double *&X, double **&T, double **&S_gg, double **&S_gS, double **&S_SS, double **&S_Sg, double **&S_NSNS, parameters param)
{
    int N = param.Nx;
    double eps = param.eps;

    //Calculation of recurrent quantities
    for (int i = 0; i < N; ++i)
    {
        for (int j = 0; j < N; ++j)
        {
            // /!\ /!\ /!\
            //Can be improved by iterative definition of Chebyshev polynomials
            T[i][j] = cheb_T(i, node(N,j));
            S_gg[i][j] = 0; //initialization of S_gg_ij
            S_gS[i][j] = 0; //initialization of S_gS_ij
            S_SS[i][j] = 0; //initialization of S_SS_ij
            S_Sg[i][j] = 0; //initialization of S_Sg_ij
            S_NSNS[i][j] = 0; //initialization of S_NSNS_ij
        }
    }
    cout << "\n-Tij calculated" << endl;

    //Calculation of Sij : time independent matrix appearing in the differential equation
    if (eps==0)
        cheb_init_S_chebint(X, T, S_gg, S_gS, S_SS, S_Sg, S_NSNS, N); //with linear scale
    else
        cheb_init_S_chebint_log(X, T, S_gg, S_gS, S_SS, S_Sg, S_NSNS, N, eps); //with log scale
    cout << "-Sij calculated" << endl;
}

//Euler method to solve the differential equation, writing solutions
void cheb_Euler(double *&X, double **&T, double *&D_g, double *&D_S, double *&D_NS, double **&S_gg, double **&S_gS, double **&S_SS, double **&S_Sg, double **&S_NSNS, ofstream &g_grid, ofstream &S_grid, ofstream &NS_grid, parameters param)
{
    double dD_g, dD_S, dD_NS;

    int N = param.Nx;

    int n_t = param.Nt;
    int n_t_w = param.Nt_w;
    double t0 = param.t0;
    double dt = param.dt;

    double Ieps = param.Ieps;
    double C_g0 = param.C_g0;
    double C_S0 = param.C_S0;
    double C_NS0 = param.C_NS0;
    bool initgrid = param.initgrid;

    //Initialisation of the solution
    if (initgrid == 1) // initialisation with grid file
    {
        string initgridname = param.initgridname;
        ifstream initgrid_file;
        initgrid_file.open(initgridname);

        double D0x = 0;
        for (int k = 0; k < N; ++k)
        {
            initgrid_file >> D0x;

            D_g[k] = C_g0*D0x; // gluon starting distribution
            D_S[k] = C_S0*D0x; // Singlet starting distribution
            D_NS[k] = C_NS0*D0x; // Non-singlet starting distribution
        }
        initgrid_file.close();
    }
    else if (Ieps == 0) // initialisation with narrow Gaussians
    {
        for (int k = 0; k < N; ++k)
        {
            D_g[k] = C_g0*D_t0_a(X[k],t0); // gluon starting distribution
            D_S[k] = C_S0*D_t0_a(X[k],t0); // Singlet starting distribution
            D_NS[k] = C_NS0*D_t0_a(X[k],t0); // Non-singlet starting distribution
        }
    }
    else // initialisation analytical solution (of the simplified BDIM)
    {
        for (int k = 0; k < N; ++k)
        {
            D_g[k] = C_g0*D_t0_e(X[k], Ieps); // gluon starting distribution
            D_S[k] = C_S0*D_t0_e(X[k], Ieps); // Singlet starting distribution
            D_NS[k] = C_NS0*D_t0_e(X[k], Ieps); // Non-singlet starting distribution
        }
    }
    // Writing of the initial solution
    write_D(X, T, D_g, D_S, D_NS, t0, g_grid, S_grid, NS_grid, param);

    // Euler method
    ///////////////
    for (int t = 1; t < n_t; ++t)
    {
        for (int k = 0; k < N; ++k)
        {
            //Calcualation of the fragmentation functions derivatives
            dD_g = 0;
            dD_S = 0;
            dD_NS = 0;
            for (int j = 0; j < N; ++j )
            {
                dD_g += dt*S_gg[k][j]*D_g[j];
                dD_g += dt*S_gS[k][j]*D_S[j];

                dD_S += dt*S_SS[k][j]*D_S[j];
                dD_S += dt*S_Sg[k][j]*D_g[j];

                dD_NS += dt*S_NSNS[k][j]*D_NS[j];
            }
            //Application of the Euler method at proper time t0+t*dt and on x=x_k
            D_g[k] = D_g[k] + dD_g;
            D_S[k] = D_S[k] + dD_S;
            D_NS[k] = D_NS[k] + dD_NS;

            //Force positivity of the coefficients
            if (param.positivity)
            {
                if (D_g[k] < 0)
                    D_g[k] = 0;
                if (D_S[k] < 0)
                    D_S[k] = 0;
                if (D_NS[k] < 0)
                    D_NS[k] = 0;
            }
        }

        // Writing the grids
        ////////////////////
        if (t%(n_t/n_t_w)==0)
            write_D(X, T, D_g, D_S, D_NS, t0+t*dt, g_grid, S_grid, NS_grid, param);
    }
}

//calculation of Sij using a linear scale for the nodes x_k
//PS : use as few loops as possible (maybe still optimizable)
void cheb_init_S_chebint(double *&X, double **&T, double **&S_gg, double **&S_gS, double **&S_SS, double **&S_Sg, double **&S_NSNS, int N)
{
    //Allocation of memory for intermediate elements
    double *intT = new double[N]; //integral of Ti(z) on [-1,1]
    double **X0k = new double*[N]; //Images of the roots of T_N
    double **Xk1 = new double*[N]; //Images of the roots of T_N

    //Integrals to calculate, see manual for more details
    double **I_g1 = new double*[N];
    double *I_g2 = new double[N];
    double I_g3 = 0;
    double **I_g4 = new double*[N];
    double **I_S1 = new double*[N];
    double *I_S2 = new double[N];
    double **I_S3 = new double*[N];

    double I; //Intermediate calculation variable

    for (int i = 0; i < N; ++i)
        X[i] = y_inv(node(N, i),0,1);

    for (int i = 0; i < N; ++i)
    {
        I_g1[i] = new double[N];
        I_g2[i] = 0;
        I_g4[i] = new double[N];
        I_S1[i] = new double[N];
        I_S2[i] = 0;
        I_S3[i] = new double[N];
        X0k[i] = new double[N];
        Xk1[i] = new double[N];

        intT[i] = cheb_int(i);

        for (int j = 0; j < N; ++j)
        {
            I_g1[i][j] = 0;
            I_g4[i][j] = 0;
            I_S1[i][j] = 0;
            I_S3[i][j] = 0;

            X0k[i][j] = y_inv(node(N, i),0,X[j]);
            Xk1[i][j] = y_inv(node(N, i),X[j],1);
        }
    }
    cout << "-Intermediate elements calculated" << endl;

    for (int k = 0; k < N; ++k)
    {
        for (int i = 0; i < N; ++i)
        {
            for (int m = 0; m < N; ++m)
            {
                for (int l = 0; l < N; ++l)
                {
                    //Calculation of all the I_ integrals using Chebyshev expansion
                    I = 1./N*T[l][m]*.5*(1-X[k])*intT[l];
                    if (l==0)
                    {
                        I_g1[k][i] += I*K_gg(Xk1[m][k])*(sqrt(Xk1[m][k])*cheb_T(i,y_(X[k]/Xk1[m][k],0,1))-Xk1[m][k]*T[i][k]);
                        I_g4[k][i] += I*K_gq(Xk1[m][k])*sqrt(Xk1[m][k])*cheb_T(i,y_(X[k]/Xk1[m][k],0,1));
                        I_S1[k][i] += I*K_qq(Xk1[m][k])*(sqrt(Xk1[m][k])*cheb_T(i,y_(X[k]/Xk1[m][k],0,1))-T[i][k]);
                        I_S3[k][i] += I*K_qg(Xk1[m][k])*sqrt(Xk1[m][k])*cheb_T(i,y_(X[k]/Xk1[m][k],0,1));
                    }
                    else
                    {
                        I_g1[k][i] += 2.*I*K_gg(Xk1[m][k])*(sqrt(Xk1[m][k])*cheb_T(i,y_(X[k]/Xk1[m][k],0,1))-Xk1[m][k]*T[i][k]);
                        I_g4[k][i] += 2.*I*K_gq(Xk1[m][k])*sqrt(Xk1[m][k])*cheb_T(i,y_(X[k]/Xk1[m][k],0,1));
                        I_S1[k][i] += 2.*I*K_qq(Xk1[m][k])*(sqrt(Xk1[m][k])*cheb_T(i,y_(X[k]/Xk1[m][k],0,1))-T[i][k]);
                        I_S3[k][i] += 2.*I*K_qg(Xk1[m][k])*sqrt(Xk1[m][k])*cheb_T(i,y_(X[k]/Xk1[m][k],0,1));
                    }

                    if (i==0)
                    {
                        if (l==0)
                        {
                            //I_g2 += 1./N*K_gg(X[m])*X[m]*T[l][m]*.5*intT[l];
                            //I_g3 += 1./N*K_qg(X[m])*X[m]*T[l][m]*.5*intT[l];
                            //I_S2 += 1./N*K_qq(X[m])*T[l][m]*.5*intT[l];

                            I_g2[k] += 1./N*K_gg(X0k[m][k])*X0k[m][k]*T[l][m]*.5*X[k]*intT[l];
                            I_S2[k] += 1./N*K_qq(X0k[m][k])*T[l][m]*.5*X[k]*intT[l];
                            if (k==0)
                                I_g3 += 1./N*K_qg(X[m])*X[m]*T[l][m]*.5*intT[l];
                        }
                        else
                        {
                            //I_g2 += 2./N*K_gg(X[m])*X[m]*T[l][m]*.5*intT[l];
                            //I_g3 += 2./N*K_qg(X[m])*X[m]*T[l][m]*.5*intT[l];
                            //I_S2 += 2./N*K_qq(X[m])*T[l][m]*.5*intT[l];

                            I_g2[k] += 2./N*K_gg(X0k[m][k])*X0k[m][k]*T[l][m]*.5*X[k]*intT[l];
                            I_S2[k] += 2./N*K_qq(X0k[m][k])*T[l][m]*.5*X[k]*intT[l];
                            if (k==0)
                                I_g3 += 2./N*K_qg(X[m])*X[m]*T[l][m]*.5*intT[l];
                        }
                    }
                }
            }
            //Calculation of the S matrices
            for (int j = 0; j < N; ++j)
            {
                if (i==0)
                {
                    S_gg[k][j] += 1./N/sqrt(X[k])*T[i][j]*I_g1[k][i];
                    S_gS[k][j] += 1./N/sqrt(X[k])*T[i][j]*I_g4[k][i];
                    S_SS[k][j] += 1./N/sqrt(X[k])*T[i][j]*I_S1[k][i];
                    S_Sg[k][j] += 1./N/sqrt(X[k])*T[i][j]*I_S3[k][i];
                    S_NSNS[k][j] += 1./N/sqrt(X[k])*T[i][j]*I_S1[k][i];
                }
                else
                {
                    S_gg[k][j] += 2./N/sqrt(X[k])*T[i][j]*I_g1[k][i];
                    S_gS[k][j] += 2./N/sqrt(X[k])*T[i][j]*I_g4[k][i];
                    S_SS[k][j] += 2./N/sqrt(X[k])*T[i][j]*I_S1[k][i];
                    S_Sg[k][j] += 2./N/sqrt(X[k])*T[i][j]*I_S3[k][i];
                    S_NSNS[k][j] += 2./N/sqrt(X[k])*T[i][j]*I_S1[k][i];
                }
            }
        }
        S_gg[k][k] -= 1./sqrt(X[k])*(I_g2[k]+I_g3);
        S_SS[k][k] -= 1./sqrt(X[k])*I_S2[k];
        S_NSNS[k][k] -= 1./sqrt(X[k])*I_S2[k];
    }

    //Memory liberation of the intermediate matrices
    for (int i = 0; i < N; ++i)
    {
        delete[] I_g1[i];
        delete[] I_g4[i];
        delete[] I_S1[i];
        delete[] I_S3[i];
        delete[] X0k[i];
        delete[] Xk1[i];
    }
    delete[] I_g1;
    delete[] I_g2;
    delete[] I_g4;
    delete[] I_S1;
    delete[] I_S2;
    delete[] I_S3;
    delete[] X0k;
    delete[] Xk1;
    delete[] intT;
}

//calculation of Sij using Chebyshev using a log scale for the nodes x_k
void cheb_init_S_chebint_log(double *&X, double **&T, double **&S_gg, double **&S_gS, double **&S_SS, double **&S_Sg, double **&S_NSNS, int N, double eps)
{
    //Allocation of memory for intermediate elements
    double *intT = new double[N]; //integral of Ti(z) on [-1,1]
    double **X0k = new double*[N]; //Images of the roots of T_N
    double **Xk1 = new double*[N]; //Images of the roots of T_N

    //Integrals to calculate, see manual for more details
    double **I_g1 = new double*[N];
    double *I_g2 = new double[N];
    double I_g3 = 0;
    double **I_g4 = new double*[N];
    double **I_S1 = new double*[N];
    double *I_S2 = new double[N];
    double **I_S3 = new double*[N];

    double I; //Intermediate calculation variable

    for (int i = 0; i < N; ++i)
        X[i] = y_log_inv(node(N, i),eps,1);

    for (int i = 0; i < N; ++i)
    {
        I_g1[i] = new double[N];
        I_g2[i] = 0;
        I_g4[i] = new double[N];
        I_S1[i] = new double[N];
        I_S2[i] = 0;
        I_S3[i] = new double[N];
        X0k[i] = new double[N];
        Xk1[i] = new double[N];

        intT[i] = cheb_int(i);

        for (int j = 0; j < N; ++j)
        {
            I_g1[i][j] = 0;
            I_g4[i][j] = 0;
            I_S1[i][j] = 0;
            I_S3[i][j] = 0;

            X0k[i][j] = y_inv(node(N, i),0,X[j]);
            Xk1[i][j] = y_inv(node(N, i),X[j],1);
        }
    }
    cout << "-Intermediate elements calculated" << endl;

    for (int k = 0; k < N; ++k)
    {
        for (int i = 0; i < N; ++i)
        {
            for (int m = 0; m < N; ++m)
            {
                for (int l = 0; l < N; ++l)
                {
                    I = 1./N*T[l][m]*.5*(1-X[k])*intT[l];
                    //Calculation of all the I_ integrals using Chebyshev expansion
                    if (l==0)
                    {
                        I_g1[k][i] += I*K_gg(Xk1[m][k])*(sqrt(Xk1[m][k])*cheb_T(i,y_log(X[k]/Xk1[m][k],eps,1))-Xk1[m][k]*T[i][k]);
                        I_g4[k][i] += I*K_gq(Xk1[m][k])*sqrt(Xk1[m][k])*cheb_T(i,y_log(X[k]/Xk1[m][k],eps,1));
                        I_S1[k][i] += I*K_qq(Xk1[m][k])*(sqrt(Xk1[m][k])*cheb_T(i,y_log(X[k]/Xk1[m][k],eps,1))-T[i][k]);
                        I_S3[k][i] += I*K_qg(Xk1[m][k])*sqrt(Xk1[m][k])*cheb_T(i,y_log(X[k]/Xk1[m][k],eps,1));
                    }
                    else
                    {
                        I_g1[k][i] += 2.*I*K_gg(Xk1[m][k])*(sqrt(Xk1[m][k])*cheb_T(i,y_log(X[k]/Xk1[m][k],eps,1))-Xk1[m][k]*T[i][k]);
                        I_g4[k][i] += 2.*I*K_gq(Xk1[m][k])*sqrt(Xk1[m][k])*cheb_T(i,y_log(X[k]/Xk1[m][k],eps,1));
                        I_S1[k][i] += 2.*I*K_qq(Xk1[m][k])*(sqrt(Xk1[m][k])*cheb_T(i,y_log(X[k]/Xk1[m][k],eps,1))-T[i][k]);
                        I_S3[k][i] += 2.*I*K_qg(Xk1[m][k])*sqrt(Xk1[m][k])*cheb_T(i,y_log(X[k]/Xk1[m][k],eps,1));
                    }

                    if (i==0)
                    {
                        if (l==0)
                        {
                            //I_g2 += 1./N*K_gg(X[m])*X[m]*T[l][m]*.5*intT[l];
                            //I_g3 += 1./N*K_qg(X[m])*X[m]*T[l][m]*.5*intT[l];
                            //I_S2 += 1./N*K_qq(X[m])*T[l][m]*.5*intT[l];

                            I_g2[k] += 1./N*K_gg(X0k[m][k])*X0k[m][k]*T[l][m]*.5*X[k]*intT[l];
                            I_S2[k] += 1./N*K_qq(X0k[m][k])*T[l][m]*.5*X[k]*intT[l];
                            if (k==0)
                                I_g3 += 1./N*K_qg(X[m])*X[m]*T[l][m]*.5*intT[l];
                        }
                        else
                        {
                            //I_g2 += 2./N*K_gg(X[m])*X[m]*T[l][m]*.5*intT[l];
                            //I_g3 += 2./N*K_qg(X[m])*X[m]*T[l][m]*.5*intT[l];
                            //I_S2 += 2./N*K_qq(X[m])*T[l][m]*.5*intT[l];

                            I_g2[k] += 2./N*K_gg(X0k[m][k])*X0k[m][k]*T[l][m]*.5*X[k]*intT[l];
                            I_S2[k] += 2./N*K_qq(X0k[m][k])*T[l][m]*.5*X[k]*intT[l];
                            if (k==0)
                                I_g3 += 2./N*K_qg(X[m])*X[m]*T[l][m]*.5*intT[l];
                        }
                    }
                }
            }
            //Calculation of the S matrices
            for (int j = 0; j < N; ++j)
            {
                if (i==0)
                {
                    S_gg[k][j] += 1./N/sqrt(X[k])*T[i][j]*I_g1[k][i];
                    S_gS[k][j] += 1./N/sqrt(X[k])*T[i][j]*I_g4[k][i];
                    S_SS[k][j] += 1./N/sqrt(X[k])*T[i][j]*I_S1[k][i];
                    S_Sg[k][j] += 1./N/sqrt(X[k])*T[i][j]*I_S3[k][i];
                    S_NSNS[k][j] += 1./N/sqrt(X[k])*T[i][j]*I_S1[k][i];
                }
                else
                {
                    S_gg[k][j] += 2./N/sqrt(X[k])*T[i][j]*I_g1[k][i];
                    S_gS[k][j] += 2./N/sqrt(X[k])*T[i][j]*I_g4[k][i];
                    S_SS[k][j] += 2./N/sqrt(X[k])*T[i][j]*I_S1[k][i];
                    S_Sg[k][j] += 2./N/sqrt(X[k])*T[i][j]*I_S3[k][i];
                    S_NSNS[k][j] += 2./N/sqrt(X[k])*T[i][j]*I_S1[k][i];
                }
            }
        }
        S_gg[k][k] -= 1./sqrt(X[k])*(I_g2[k]+I_g3);
        S_SS[k][k] -= 1./sqrt(X[k])*I_S2[k];
        S_NSNS[k][k] -= 1./sqrt(X[k])*I_S2[k];
    }

    //Memory liberation of the intermediate matrices
    for (int i = 0; i < N; ++i)
    {
        delete[] I_g1[i];
        delete[] I_g4[i];
        delete[] I_S1[i];
        delete[] I_S3[i];
        delete[] X0k[i];
        delete[] Xk1[i];
    }
    delete[] I_g1;
    delete[] I_g2;
    delete[] I_g4;
    delete[] I_S1;
    delete[] I_S2;
    delete[] I_S3;
    delete[] X0k;
    delete[] Xk1;
    delete[] intT;
}

// Calculation of the integral of D_g+D_S+D_NS or (D_g+D_S+D_NS)/x over x
double normalize_D(double *&X, double **&T, double D_g[], double D_S[], double D_NS[], parameters param)
{
    int N = param.Nx;
    double eps = param.eps;
    int norm = param.norm;

    double intD = 0;
    double factCheb[N][N]; //Coefficient depending only on i,j

    for (int i = 0; i < N; ++i)
    {
        for (int j = 0; j < N; ++j)
        {
            if (eps == 0)
                factCheb[i][j] = T[i][j]*cheb_int(i)/N;
            else
                factCheb[i][j] = -T[i][j]*pow(eps,(1-node(N, j))/2.)*log(eps)*cheb_int(i)/N;
        }
    }

    if (norm == 1)
    {
        for (int i = 0; i < N; ++i)
        {
            for (int j = 0; j < N; ++j)
            {
                intD +=1.*D_g[j]*factCheb[i][j];
                intD +=1.*D_S[j]*factCheb[i][j];
                //intD +=1.*D_NS[j]*factCheb[i][j];
            }
            intD -=.5*D_g[i]*factCheb[0][i];
            intD -=.5*D_S[i]*factCheb[0][i];
            //intD -=.5*D_NS[i]*factCheb[0][i];
        }
    }
    else if (norm == 2)
    {
        for (int i = 0; i < N; ++i)
        {
            for (int j = 0; j < N; ++j)
            {
                intD +=1.*D_g[j]/X[j]*factCheb[i][j];
                intD +=1.*D_S[j]/X[j]*factCheb[i][j];
                //intD +=1.*D_NS[j]/X[j]*factCheb[i][j];
            }
            intD -=.5*D_g[i]/X[i]*factCheb[0][i];
            intD -=.5*D_S[i]/X[i]*factCheb[0][i];
            //intD -=.5*D_NS[i]/X[i]*factCheb[0][i];
        }
    }

    return intD;
}

// Add the lines "X[k] D[k]" to g_grid
void write_D(double *&X, double **&T, double *&D_g, double *&D_S, double *&D_NS, double t, ofstream &g_grid, ofstream &S_grid, ofstream &NS_grid, parameters param)
{
    double intD_g, intD_S, intD_NS;

    // Not normalize solutions
    if (param.norm == 0)
    {//Writing of the solutions on grids
        for (int k = 0; k < param.Nx; ++k)
        {
            g_grid  << t << "\t" << X[k] << "\t" << D_g[k] << endl;
            S_grid  << t << "\t" << X[k] << "\t" << D_S[k] << endl;
            NS_grid << t << "\t" << X[k] << "\t" << D_NS[k] << endl;
        }
    }
    // Normalized fragmentation functions
    else
    {
        if ((param.norm == 1)||(param.norm == 2))
        {
            //Calculating the integral of (D_g+D_S+D_NS) or (D_g+D_S+D_NS)/x depending on norm
            intD_g = normalize_D(X, T, D_g, D_S, D_NS, param);
            intD_S = intD_g;
            intD_NS = intD_g;
        }
        else if ((param.norm == 3)||(param.norm == 4))
        {
            //Calculating the integral of D_f or D_f/x depending on norm
            intD_g = normalize_D(X, T, D_g, param);
            intD_S = normalize_D(X, T, D_S, param);
            intD_NS = normalize_D(X, T, D_NS, param);

            //the distributions are normalized except when they are null
            if (intD_g == 0)
                intD_g = 1.;
            if (intD_S == 0)
                intD_S = 1.;
            if (intD_NS == 0)
                intD_NS = 1.;
        }

        if (param.printintsig)
            cout << "Time : " << t << "\t\tTotal integral : " << intD_g << endl;

        //Writing of the normalized solution on grids
        for (int k = 0; k < param.Nx; ++k)
        {
            g_grid  << t << "\t" << X[k] << "\t" << D_g[k]/intD_g << endl;
            S_grid  << t << "\t" << X[k] << "\t" << D_S[k]/intD_S << endl;
            NS_grid << t << "\t" << X[k] << "\t" << D_NS[k]/intD_NS << endl;
        }
    }
}

// Calculation of the integral of D or D/x over x
double normalize_D(double *&X, double **&T, double D[], parameters param)
{
    int N = param.Nx;
    double eps = param.eps;
    int norm = param.norm;

    double intD = 0;
    double factCheb[N][N]; //Coefficient depending only on i,j

    for (int i = 0; i < N; ++i)
    {
        for (int j = 0; j < N; ++j)
        {
            if (eps == 0)
                factCheb[i][j] = T[i][j]*cheb_int(i)/N;
            else
                factCheb[i][j] = -T[i][j]*pow(eps,(1-node(N, j))/2.)*log(eps)*cheb_int(i)/N;
        }
    }

    if (norm == 3)
    {
        for (int i = 0; i < N; ++i)
        {
            for (int j = 0; j < N; ++j)
                intD +=1.*D[j]*factCheb[i][j];
            intD -=.5*D[i]*factCheb[0][i];
        }
    }
    else if (norm == 4)
    {
        for (int i = 0; i < N; ++i)
        {
            for (int j = 0; j < N; ++j)
                intD +=1.*D[j]/X[j]*factCheb[i][j];
            intD -=.5*D[i]/X[i]*factCheb[0][i];
        }
    }

    return intD;
}

#endif // CHEB_METHOD_G_S_NS_H
