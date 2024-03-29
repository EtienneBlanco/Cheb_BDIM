/*=====================================================\
|                                                      |
|                      Parameters                      |
|	                                                   |
\=====================================================*/


// On this file, a parameters class is defined with some default values and routines to read them from a file.

#ifndef CHEB_PARAMETERS_H
#define CHEB_PARAMETERS_H

#include <string>
#include <sstream>
#include <math.h>

double paramstr2double(string str); //Function that find a double in a parameter line and return it
int paramstr2int(string str); //Function that find an integer in a parameter and return it
bool paramstr2bool(string str); //Function that find a boolean in a parameter and return it
string paramstr2str(string str); //Function that find a string in a parameter and return it

//Constants parameters used
class parameters
{
    public:
        //Result options
        int norm; //normalization method (0 : no normalization, 1 : normalization of D(x), 2 : normalization of N(x))
        bool positivity; // Consider or not negative values in distribution grids (if equal to 1, all negative values are set to 0)

        //Evolution used
        int evo; //type of evolution : pure gluon (0) or with quarks (1)

        //Chebyshev method
        int Nx; //Number of Chebychev nodes for the integration in x
        double eps; //low-x cut-off when using logarithmic bijection

        //Euler method
        int Nt; // Number point in time for the Euler method
        int Nt_w; // number of points in time used for the result grid
        double t0; // initial time of the evolution
        double dt; //  time step

        //Initial distributions
        double Ieps; // Width of the initial distribution : a Gaussian peaked on x=1 or the analytical solution to the simplified BDIM if set to 0
        double C_g0; // Initial gluon distribution (coefficient in front of the Gaussian)
        double C_S0; // Initial quark singlet distribution (coefficient in front of the Gaussian)
        double C_NS0; // Initial quark non-singlet distribution (coefficient in front of the Gaussian)
        bool initgrid; // Set to 0, the initial distribution is described by the previous parameters, set to 1 its uses a grid
        string initgridname; // name of the grid to use when initgrid is set to 1


        string gridname; //name to be used for the result grid

        //internal parameters
        //bool printig = false; // option for printing D(x, t) when put in grids
        bool printintsig = false; // option for printing the integral of D(x) when calculated
        bool paramsfromfile;

        void from_file(string filename); // initialize parameters from a parameter file
        void default_(); // initialize parameters with default values
        void show(); // Display parameters in the class
        void write(); // Write parameters in a file

        //Parameters trackers
        //(track if parameters are set from default values or from file)
        int N_param = 16; //number of initial parameters
        bool *paramfromfile = new bool[N_param]; //parameters that are set through file
        int Nparamfromfile = 0; //number of parameters set through parameter file
        string *paramfromfile_name = new string[N_param]; //name (string) of the parameters
};


//function reading the parameters file
void parameters::from_file(string filename)
{
    string line;

    ifstream parameters_file;
    parameters_file.open(filename);

    default_(); // Initialize the parameters with default values in case some are not defined in file

    if (parameters_file.is_open())
    {
        while (getline(parameters_file, line))
        {
            if (!(line[0]=='/' && line[1]=='/'))
            {

                if (line[0]=='e' && line[1]=='v' && line[2]=='o')
                {
                    evo = paramstr2int(line);
                    paramfromfile[2] = true;
                    Nparamfromfile++;
                    if (!(evo == 0) && !(evo == 1)&& !(evo == 2))
                    {
                        evo = 1;
                        paramfromfile[2] = false;
                        Nparamfromfile--;
                    }
                }

                if (line[0]=='n' && line[1]=='o' && line[2]=='r' && line[3]=='m')
                {
                    norm = paramstr2int(line);
                    paramfromfile[0] = true;
                    Nparamfromfile++;
                    if (!(norm == 0) && !(norm == 1) && !(norm == 2) && !(norm == 3) && !(norm == 4))
                    {
                        norm = 0;
                        paramfromfile[0] = true;
                        Nparamfromfile--;
                    }
                }
                if (line[0]=='p' && line[1]=='o' && line[2]=='s' && line[3]=='i' && line[4]=='t' && line[5]=='i' && line[6]=='v' && line[7]=='i' && line[8]=='t' && line[9]=='y')
                {
                    positivity = paramstr2bool(line);
                    paramfromfile[1] = true;
                    Nparamfromfile++;
                }

                if (line[0]=='N' && line[1]=='x')
                {
                    Nx = paramstr2int(line);
                    paramfromfile[3] = true;
                    Nparamfromfile++;
                    if (Nx == 0)
                    {
                        Nx = 1;
                        paramfromfile[3] = true;
                        Nparamfromfile--;
                    }
                }
                if (line[0]=='e' && line[1]=='p' && line[2]=='s')
                {
                    eps = paramstr2double(line);
                    paramfromfile[4] = true;
                    Nparamfromfile++;
                }

                if (line[0]=='t'&& line[1]=='a' && line[2]=='u' && line[3]=='0')
                {
                    t0 = paramstr2double(line);
                    paramfromfile[7] = true;
                    Nparamfromfile++;
                }
                if (line[0]=='d' && line[1]=='t' && line[2]=='a' && line[3]=='u')
                {
                    dt = paramstr2double(line);
                    paramfromfile[8] = true;
                    Nparamfromfile++;
                }
                if (line[0]=='N' && line[1]=='t' && line[2]=='a' && line[3]=='u')
                {
                    Nt = paramstr2double(line);
                    paramfromfile[5] = true;
                    Nparamfromfile++;
                }
                if (line[0]=='F' && line[1]=='t' && line[2]=='a' && line[3]=='u')
                {
                    double Ft = paramstr2double(line);
                    Nt = ceil((Ft-t0)/dt);
                    cout << "Ntau obtained using : Ftau = " << Ft << ", tau0 = "<< t0 << ", dtau = "<< dt << endl;
                    paramfromfile[5] = true;
                    Nparamfromfile++;
                }
                if (line[0]=='N' && line[1]=='t' && line[2]=='_' && line[3]=='w')
                {
                    Nt_w = paramstr2double(line);
                    paramfromfile[6] = true;
                    Nparamfromfile++;
                }

                if (line[0]=='I' && line[1]=='e' && line[2]=='p' && line[3]=='s')
                {
                    Ieps = paramstr2double(line);
                    paramfromfile[9] = true;
                    Nparamfromfile++;
                }
                if (line[0]=='C' && line[1]=='_' && line[2]=='g' && line[3]=='0')
                {
                    C_g0 = paramstr2double(line);
                    paramfromfile[10] = true;
                    Nparamfromfile++;
                }
                if (line[0]=='C' && line[1]=='_' && line[2]=='S' && line[3]=='0')
                {
                    C_S0 = paramstr2double(line);
                    paramfromfile[11] = true;
                    Nparamfromfile++;
                }
                if (line[0]=='C' && line[1]=='_' && line[2]=='N' && line[3]=='S' && line[4]=='0')
                {
                    C_NS0 = paramstr2double(line);
                    paramfromfile[12] = true;
                    Nparamfromfile++;
                }
                if (line[0]=='i' && line[1]=='n' && line[2]=='i' && line[3]=='t' && line[4]=='g' && line[5]=='r' && line[6]=='i' && line[7]=='d')
                {
                    initgrid = paramstr2bool(line);
                    paramfromfile[13] = true;
                    Nparamfromfile++;
                }
                if (line[0]=='i' && line[1]=='n' && line[2]=='i' && line[3]=='t' && line[4]=='g' && line[5]=='r' && line[6]=='i' && line[7]=='d' && line[8]=='n' && line[9]=='a' && line[10]=='m' && line[11]=='e')
                {
                    initgridname = paramstr2str(line);
                    ifstream initgrid_file;
                    initgrid_file.open(initgridname);
                    paramfromfile[14] = true;
                    Nparamfromfile++;

                    if (!(initgrid_file.is_open())) //verification that the file has been found
                    {
                        initgrid = 0; //if not, initial distribution is set from parameters
                        paramfromfile[14] = false;
                        Nparamfromfile--;
                    }
                    initgrid_file.close();
                }

                if (line[0]=='g' && line[1]=='r' && line[2]=='i' && line[3]=='d' && line[4]=='n' && line[5]=='a' && line[6]=='m' && line[7]=='e')
                {
                    gridname = paramstr2str(line);
                    paramfromfile[15] = true;
                    Nparamfromfile++;
                }
            }
        }
        paramsfromfile = true;
    }
    else
        paramfromfile = false;
    show(); //display the parameters that will be used
    write(); //save them in a file

    parameters_file.close();
}

//Default setting for the parameters
void parameters::default_()
{
    evo = 1;

    norm = 0;
    positivity = 0;

    Nx = 50;
    eps = 1e-4;

    Nt = 1e6;
    Nt_w = 1e3;
    t0 = 0;
    dt = 1e-6;

    Ieps = 1e-2;
    C_g0 = 1.;
    C_S0 = 0;
    C_NS0 = 0;
    initgrid = false;
    initgridname = "resultgrid_";

    gridname = "Results/D_cheb";


    //initialization of the trackers
    for (int i = 0; i < N_param; ++i)
        paramfromfile[i] = false;
    paramfromfile_name[0] = "norm";
    paramfromfile_name[1] = "positivity";
    paramfromfile_name[2] = "evo";
    paramfromfile_name[3] = "Nx";
    paramfromfile_name[4] = "eps";
    paramfromfile_name[5] = "Ntau";
    paramfromfile_name[6] = "Nt_w";
    paramfromfile_name[7] = "tau0";
    paramfromfile_name[8] = "dtau";
    paramfromfile_name[9] = "Ieps";
    paramfromfile_name[10] = "C_g0";
    paramfromfile_name[11] = "C_S0";
    paramfromfile_name[12] = "C_NS0";
    paramfromfile_name[13] = "initgrid";
    paramfromfile_name[14] = "initgridname";
    paramfromfile_name[15] = "gridname";
}

//Display the parameters
void parameters::show()
{
    cout << "\n#//////////////////#" << endl;
    cout << "#//  Parameters  //#" << endl;
    cout << "#/////////////////////////////////////////////#" << endl;
    if (paramfromfile)
        cout << "# Parameters are set to the following values : " << endl;
    else
    {
        cout << "# /!\\ 'Parameters' file not found" << endl;
        cout << "# Parameters are set to the following default values : " << endl;
    }
    cout << "#" << endl;

    cout << "#Evolution equations :" << endl;
    cout << "# -Type : ";
    if (evo == 0)
        cout << "Pure gluons" << endl;
    else if (evo == 1)
        cout << "Quarks and gluons" << endl;
    else if (evo == 2)
        cout << "Pure gluons with simplified kernel" << endl;

    cout << "#" << endl;

    cout << "#Results :" << endl;
    cout << "# -Normalization : ";
    if (norm == 0)
        cout << "none" << endl;
    else if (norm == 1)
        cout << "D(x) is normalized" << endl;
    else if (norm == 2)
        cout << "N(x) is normalized" << endl;
    else if (norm == 3)
        cout << "each D_f(x) are normalized" << endl;
    else if (norm == 4)
        cout << "each N_f(x) are normalized" << endl;

    cout << "# -Positivity : ";
    if (positivity)
        cout << "Yes" << endl;
    else
        cout << "No" << endl;

    cout << "#" << endl;

    cout << "#Chebyshev method :" << endl;
    cout << "# -Nx = " << Nx << endl;
    if (eps>0)
        cout << "# -logarithmic bijection y: [" << eps << ", 1] -> [-1, 1]" << endl;
    else if (eps==0)
        cout << "# -linear bijection y: [0, 1] -> [-1, 1]" << endl;

    cout << "#" << endl;

    cout << "#Euler method :" << endl;
    cout << "# -Nt = " << Nt << endl;
    cout << "# -Nt_w = " << Nt_w << endl;
    cout << "# -t0 = " << t0 << endl;
    cout << "# -dt = " << dt << endl;

    cout << "#" << endl;

    cout << "#Initial distributions :" << endl;
    if (initgrid)
        cout << "# Initial solution set from the grid : " << initgridname << endl;
    else
    {
        cout << "# -Ieps = " << Ieps << endl;
        if (evo == 0)
            cout << "# -C_g0 = " << "1" << endl;
        if (evo == 1)
        {
            cout << "# -C_g0 = " << C_g0 << endl;
            cout << "# -C_S0 = " << C_S0 << endl;
            cout << "# -C_NS0 = " << C_NS0 << endl;
        }
    }

    cout << "#" << endl;

    cout << "# -gridname = " << gridname << endl;
    cout << "#/////////////////////////////////////////////#" << endl;

    if ((paramfromfile)&&(Nparamfromfile<N_param))
    {
        cout << "\n\nNparamfromfile-N_param = " << Nparamfromfile-N_param+1 << endl;
        cout << "\n\n/!\\ The following parameters have been set with default value (either not found / not set / not recognized / no proper argument) : " << endl;
        cout << "-";
        int Nparamprint = 0;
        for (int i = 0; i < N_param; ++i)
        {
            if (!paramfromfile[i])
            {
                cout << paramfromfile_name[i];
                Nparamprint++;
                if (Nparamprint < N_param-Nparamfromfile+1)
                    cout << ", ";
                else
                    cout << endl << endl << endl;
            }
        }
    }
}

//Save the parameters in a file
void parameters::write()
{
    ofstream param_file;
    param_file.open (gridname+"_parameters", ios::out);

    param_file << "\n#//////////////////#" << endl;
    param_file << "#//  Parameters  //#" << endl;
    param_file << "#/////////////////////////////////////////////#" << endl;
    param_file << "# Parameters are set to the following values : " << endl;
    param_file << "#" << endl;

    param_file << "#Evolution equations :" << endl;
    param_file << "# -Type : ";
    if (evo == 0)
        param_file << "Pure gluons" << endl;
    else if (evo == 1)
        param_file << "Quarks and gluons" << endl;
    else if (evo == 2)
        param_file << "Pure gluons with simplified kernel" << endl;

    param_file << "#" << endl;

    param_file << "#Results :" << endl;
    param_file << "# -Normalization : ";
    if (norm == 0)
        param_file << "none" << endl;
    else if (norm == 1)
        param_file << "D(x) is normalized" << endl;
    else if (norm == 2)
        param_file << "N(x) is normalized" << endl;
    else if (norm == 3)
        param_file << "each D_f(x) are normalized" << endl;
    else if (norm == 4)
        param_file << "each N_f(x) are normalized" << endl;

    param_file << "# -Positivity : ";
    if (positivity)
        param_file << "Yes" << endl;
    else
        param_file << "No" << endl;

    param_file << "#" << endl;

    param_file << "#Chebyshev method :" << endl;
    param_file << "# -Nx = " << Nx << endl;
    param_file << "# -eps = " << eps << endl;

    param_file << "#" << endl;

    param_file << "#Euler method :" << endl;
    param_file << "# -Nt = " << Nt << endl;
    param_file << "# -Nt_w = " << Nt_w << endl;
    param_file << "# -t0 = " << t0 << endl;
    param_file << "# -dt = " << dt << endl;

    param_file << "#" << endl;

    param_file << "#Initial distributions :" << endl;
    if (initgrid)
        param_file << "# Initial solution set from the grid : " << initgridname << endl;
    else
    {
        param_file << "# -Ieps = " << Ieps << endl;
        if (evo == 0)
            param_file << "# -C_g0 = " << "1" << endl;
        if (evo == 1)
        {
            param_file << "# -C_g0 = " << C_g0 << endl;
            param_file << "# -C_S0 = " << C_S0 << endl;
            param_file << "# -C_NS0 = " << C_NS0 << endl;
        }
    }

    param_file << "#" << endl;

    param_file << "# -gridname = " << gridname << endl;
    param_file << "#/////////////////////////////////////////////#" << endl;


    param_file.close();
}

double paramstr2double(string str)
{
    istringstream iss(str);
    string paramname, sign, rest;
    double param;

    iss >> paramname >> sign >> param >> rest;
    return param;
}

int paramstr2int(string str)
{
    istringstream iss(str);
    string paramname, sign, rest;
    int param;

    iss >> paramname >> sign >> param >> rest;
    return param;
}

bool paramstr2bool(string str)
{
    istringstream iss(str);
    string paramname, sign, rest;
    bool param;

    iss >> paramname >> sign >> param >> rest;
    return param;
}

string paramstr2str(string str)
{
    istringstream iss(str);
    string paramname, sign, rest;
    string param;

    iss >> paramname >> sign >> param >> rest;
    return param;
}

#endif // CHEB_PARAMETERS_H
