/*======================================================\
|                                                       |
|   Resolution method based on Chebyshev polynomials    |
|	                                                    |
\======================================================*/


// On this file are gathered the functions directly related to the resolution method.

#ifndef CHEB_METHOD_G_H
#define CHEB_METHOD_G_H

#include <math.h>

////////////////////////////////////////
//prototype of the used functions     //
////////////////////////////////////////
void cheb_init_S_chebint_g(double *&X, double **&T, double **&S_gg, int N);
void cheb_init_S_chebint_glog(double *&X, double **&T, double **&S_gg, int N, double eps);

void cheb_Euler_g(double *&X, double **&T, double *&D_g, double **&S_gg, ofstream &g_grid, parameters param);

double normalize_Dg(double *&X, double **&T, double D[], parameters param);
////////////////////////////////////////


//Function that initialize useful elements
void cheb_initilization_g(double *&X, double **&T, double **&S_gg, parameters param)
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
        }
    }
    cout << "-Tij calculated" << endl;

    //Calculation of Sij : time independent matrix appearing in the differential equation
    if (eps==0)
        cheb_init_S_chebint_g(X, T, S_gg, N); //with linear scale
    else
        cheb_init_S_chebint_glog(X, T, S_gg, N, eps); //with log scale
    cout << "-Sij calculated" << endl;
}

//Euler method to solve the differential equation, writing solutions
void cheb_Euler_g(double *&X, double **&T, double *&D_g, double **&S_gg, ofstream &g_grid, parameters param)
{
    double dD_g, intD_g;

    int norm = param.norm;

    int N = param.Nx;

    int n_t = param.Nt;
    int n_t_w = param.Nt_w;
    double t0 = param.t0;
    double dt = param.dt;

    double Ieps = param.Ieps;

    for (int t = 0; t < n_t; ++t)
    {
        // Euler method
        ///////////////
        for (int k = 0; k < N; ++k)
        {
            if (t==0)
                D_g[k] = D_t0_e(X[k], Ieps);
            else
            {//Calcualation of the fragmentation functions derivatives
                dD_g = 0;
                for (int j = 0; j < N; ++j )
                    dD_g += dt*S_gg[k][j]*D_g[j];
                //Application of the Euler method at proper time t0+t*dt and on x=x_k
                D_g[k] = D_g[k] + dD_g;

                //Force positivity of the coefficients
                if ((param.positivity) && ((D_g[k] < 0)))
                    D_g[k] = 0;
            }
        }

        // Writing the grids
        ////////////////////
        if (t%(n_t/n_t_w)==0)
        {
            // Not normalize solutions
            if (norm == 0)
            {//Writing of the solutions on grids
                for (int k = 0; k < N; ++k)
                    g_grid << t0+t*dt << "\t" << X[k] << "\t" << D_g[k] << endl;
            }
            // Normalized fragmentation functions
            else
            {

                intD_g = normalize_Dg(X, T, D_g, param);

                if (param.printintsig)
                    cout << "Time : " << t0+t*dt << "\t\tTotal integral : " << intD_g << endl;

                //Writing of the normalized solution on grids
                for (int k = 0; k < N; ++k)
                    g_grid << t0+t*dt << "\t" << X[k] << "\t" << D_g[k]/intD_g << endl;
            }
        }
    }
}

//calculation of Sij using a linear scale for the nodes x_k
//PS : use as few loops as possible (maybe still optimizable)
void cheb_init_S_chebint_g(double *&X, double **&T, double **&S_gg, int N)
{
    //Allocation of memory for intermediate elements
    double *intT = new double[N]; //integral of Ti(z) on [-1,1]
    double **X0k = new double*[N]; //Images of the roots of T_N
    double **Xk1 = new double*[N]; //Images of the roots of T_N

    //Integrals to calculate, see manual for more details
    double **I_g1 = new double*[N];
    double *I_g2 = new double[N];
    double I_g3 = 0;

    double I; //Intermediate calculation variable

    for (int i = 0; i < N; ++i)
        X[i] = y_inv(node(N, i),0,1);

    for (int i = 0; i < N; ++i)
    {
        I_g1[i] = new double[N];
        I_g2[i] = 0;
        X0k[i] = new double[N];
        Xk1[i] = new double[N];

        intT[i] = cheb_int(i);

        for (int j = 0; j < N; ++j)
        {
            I_g1[i][j] = 0;

            X0k[i][j] = y_inv(node(N, i),0,X[j]);
            Xk1[i][j] = y_inv(node(N, i),X[j],1);
        }
    }
    cout << "\n-Intermediate elements calculated";

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
                        I_g1[k][i] += I*K_gg(Xk1[m][k])*(sqrt(Xk1[m][k])*cheb_T(i,y_(X[k]/Xk1[m][k],0,1))-Xk1[m][k]*T[i][k]);
                    else
                        I_g1[k][i] += 2*I*K_gg(Xk1[m][k])*(sqrt(Xk1[m][k])*cheb_T(i,y_(X[k]/Xk1[m][k],0,1))-Xk1[m][k]*T[i][k]);

                    if (i==0)
                    {
                        if (l==0)
                        {
                            I_g2[k] += 1./N*K_gg(X0k[m][k])*X0k[m][k]*T[l][m]*.5*X[k]*intT[l];
                            if (k==0)
                                I_g3 += 1./N*K_qg(X[m])*X[m]*T[l][m]*.5*intT[l];
                        }
                        else
                        {
                            I_g2[k] += 2./N*K_gg(X0k[m][k])*X0k[m][k]*T[l][m]*.5*X[k]*intT[l];
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
                    S_gg[k][j] += 1./N/sqrt(X[k])*T[i][j]*I_g1[k][i];
                else
                    S_gg[k][j] += 2./N/sqrt(X[k])*T[i][j]*I_g1[k][i];
            }
        }
        S_gg[k][k] -= 1./sqrt(X[k])*(I_g2[k]+I_g3);
    }

    //Memory liberation of the intermediate matrices
    for (int i = 0; i < N; ++i)
    {
        delete[] I_g1[i];
        delete[] X0k[i];
        delete[] Xk1[i];
    }
    delete[] I_g1;
    delete[] I_g2;
    delete[] X0k;
    delete[] Xk1;
    delete[] intT;
}

//calculation of Sij using Chebyshev using a log scale for the nodes x_k
void cheb_init_S_chebint_glog(double *&X, double **&T, double **&S_gg, int N, double eps)
{
    //Allocation of memory for intermediate elements
    double *intT = new double[N]; //integral of Ti(z) on [-1,1]
    double **X0k = new double*[N]; //Images of the roots of T_N
    double **Xk1 = new double*[N]; //Images of the roots of T_N

    //Integrals to calculate, see manual for more details
    double **I_g1 = new double*[N];
    double *I_g2 = new double[N];
    double I_g3 = 0;

    double I; //Intermediate calculation variable

    for (int i = 0; i < N; ++i)
        X[i] = y_log_inv(node(N, i),eps,1);

    for (int i = 0; i < N; ++i)
    {
        I_g1[i] = new double[N];
        I_g2[i] = 0;
        X0k[i] = new double[N];
        Xk1[i] = new double[N];

        intT[i] = cheb_int(i);

        for (int j = 0; j < N; ++j)
        {
            I_g1[i][j] = 0;

            X0k[i][j] = y_inv(node(N, i),0,X[j]);
            Xk1[i][j] = y_inv(node(N, i),X[j],1);
        }
    }
    cout << "\n-Intermediate elements calculated";

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
                        I_g1[k][i] += I*K_gg(Xk1[m][k])*(sqrt(Xk1[m][k])*cheb_T(i,y_log(X[k]/Xk1[m][k],eps,1))-Xk1[m][k]*T[i][k]);
                    else
                        I_g1[k][i] += 2*I*K_gg(Xk1[m][k])*(sqrt(Xk1[m][k])*cheb_T(i,y_log(X[k]/Xk1[m][k],eps,1))-Xk1[m][k]*T[i][k]);

                    if (i==0)
                    {
                        if (l==0)
                        {
                            I_g2[k] += 1./N*K_gg(X0k[m][k])*X0k[m][k]*T[l][m]*.5*X[k]*intT[l];
                            if (k==0)
                                I_g3 += 1./N*K_qg(X[m])*X[m]*T[l][m]*.5*intT[l];
                        }
                        else
                        {
                            I_g2[k] += 2./N*K_gg(X0k[m][k])*X0k[m][k]*T[l][m]*.5*X[k]*intT[l];
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
                    S_gg[k][j] += 1./N/sqrt(X[k])*T[i][j]*I_g1[k][i];
                else
                    S_gg[k][j] += 2./N/sqrt(X[k])*T[i][j]*I_g1[k][i];
            }
        }
        S_gg[k][k] -= 1./sqrt(X[k])*(I_g2[k]+I_g3);
    }

    //Memory liberation of the intermediate matrices
    for (int i = 0; i < N; ++i)
    {
        delete[] I_g1[i];
        delete[] X0k[i];
        delete[] Xk1[i];
    }
    delete[] I_g1;
    delete[] I_g2;
    delete[] X0k;
    delete[] Xk1;
    delete[] intT;
}

// Calculation of the integral of D or D/x over x
double normalize_Dg(double *&X, double **&T, double D[], parameters param)
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

    if ((norm == 1) || (norm == 3))
    {
        for (int i = 0; i < N; ++i)
        {
            for (int j = 0; j < N; ++j)
                intD +=1.*D[j]*factCheb[i][j];
            intD -=.5*D[i]*factCheb[0][i];
        }
    }
    else if ((norm == 2) || (norm == 4))
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

#endif // CHEB_METHOD_G_H
