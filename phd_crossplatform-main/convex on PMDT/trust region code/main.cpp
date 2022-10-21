#include <dace/dace.h>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <gsl/gsl_multiroots.h>
#include <fmt/printf.h>
#include <gsl/gsl_vector.h>
#include "RK78.h"
#include "constants.h"
#include <chrono>
#include "dynamics.h"
#include "osctomean.h"

using namespace std::chrono;

int main()
{
    parameters params = createParams();
    double cst[9];

    cst[3] = params.Tmax;
    cst[4] = params.mstart;
    cst[5] = params.mu_earth;
    cst[6] = params.J2;
    cst[7] = params.Re;

    if (params.method < 2)
    {
        cst[8] = 1.0;
        cout << "Method 1 running" << endl;
    }
    else
    {
        cst[8] = 2.0;
        cout << "Method 2 running " << endl;
    }

    typedef DACE::DA state_type; // define state type as DA
    DACE::DA::init(2, 9);        // initialise DACE for X order in N variables
    DACE::AlgebraicVector<state_type> x0(6), fx(6), dv(3), dx(3);
    DACE::AlgebraicMatrix<state_type> Jacobian(6, 9), xval(params.N, 9);
    DACE::DA derivative(1), result(1), b;
    DACE::AlgebraicMatrix<double> STMgamma(6, 9), Jbar(6, 9);
    double statesGuess[params.N + 1][9], xarray[9], constant[6], B, normJbar, time, maxval, eclipse;
    vector<double> zerostate(9);
    vector<state_type> zerostateDA(9);

    fstream iter, outputfile, trustRegion, dxmax;
    int i = 0;

    // obtain reference trajectory from text file.
    iter.open("iter.txt", ios::in);
    if (iter.is_open())
    {
        string tp;
        while (getline(iter, tp))
        {
            removeDupWord(tp, xarray);
            for (unsigned int j = 0; j < 9; j++)
                statesGuess[i][j] = xarray[j];
            i++;
        }
        iter.close();
    }

    outputfile.open("outputcpp.txt", fstream::out);
    trustRegion.open("trustRegion.txt", fstream::out);
    dxmax.open("dxmax.txt", fstream::out);

    for (int k = 0; k < params.N; k++)
    // for (int k = 0; k <1; k++)
    {
        // update time.
        time += params.dt;

        for (unsigned int i = 0; i < 6; i++)
        {
            x0[i] = statesGuess[k][i] + DACE::DA(i + 1);

            if (i < 3)
                cst[i] = statesGuess[k][i + 6];
        }

        if (cst[8] < 2)
        {
            fx = RK78(6, x0, 0.0, params.dt, cst, DynamicsCart);
        }
        else
        {

            if (params.eclipses > 0)
            {
                eclipse = 1.0; // 1.00001 - getEclipse(time-0.5*params.dt,  x0, 1000.0 ,1.0,params.LU, params.Rs, params.Re);
                // cout << "nu eclipse = " << eclipse<< endl;
            }
            else
            {
                eclipse = 1.0;
            }
            for (unsigned int i = 0; i < 3; i++)
            {
                dv[i] = eclipse * (cst[i] + DACE::DA(i + 7)); //
                x0[i + 3] += dv[i];
            }
            // cout << x0 << endl;
            fx = RK78(6, x0, 0.0, params.dt, cst, DynamicsCart);
        }
        
      
            
 
          
        
        

        // cout << fx<< endl;

        STMgamma = fx.linear();

        // cout <<STMgamma<<endl;

        for (unsigned int i = 0; i < 6; i++)
            constant[i] = DACE::cons(fx[i]);

        // output stm and constant term
        for (unsigned int j = 0; j < 6; j++) // Prints row of x
        {
            for (unsigned int i = 0; i < 9; i++)
                outputfile << STMgamma.at(j, i) << "\t";

            outputfile << constant[j] << "\t";
            outputfile << std::endl;
        }

        // calculate the jacobian and jbar norm
        normJbar = 0;
        for (unsigned int i = 0; i < 6; i++)
        {
            for (unsigned int j = 1; j < 10; j++)
            {
                Jacobian.at(i, j - 1) = fx[i].deriv(j);
                Jbar.at(i, j - 1) = fx[i].deriv(j).eval(zerostate);
                normJbar += abs(Jbar.at(i, j - 1));
            }
        }

        // cout << normJbar << endl;

        // analyse the dependency on each variable individually
        for (int l = 0; l < 9; l++)
        {
            zerostateDA = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
            zerostateDA[l] = DACE::DA(l + 1);
            b = 0;
            for (unsigned int i = 0; i < 6; i++)
            {
                for (unsigned int j = 1; j < 10; j++)
                {
                    // cout << Jacobian.at(i, j - 1).eval(zerostateDA) << endl;
                    // cout << Jacobian.at(i, j - 1).eval(zerostateDA)<< endl;
                    b += abs(Jacobian.at(i, j - 1).eval(zerostateDA).deriv(l + 1));
                }
            }

         
            // cout <<  "nu" << b/normJbar << endl;
            // calculate x that corresponds to the given nu
            // xval.at(k, l) = params.nu * normJbar / b;
            xval.at(k, l) = params.nu * normJbar / b; //
            // cout << (xval.at(k, l)) << endl;
            // cout << DACE::cons(xval.at(k, l)) << endl;
            dxmax << DACE::cons(xval.at(k, l)) << "\t";

            // cout << b << endl;
        }
        dxmax << std::endl;
    }
    outputfile.close();
    dxmax.close();

    // find the max variation in each element to provide the required nu.
    for (unsigned int i = 0; i < 9; i++)
    {
        maxval = 0.0;
        for (unsigned int k = 0; k < params.N; k++)
            if (DACE::cons(xval.at(k, i)) > (maxval))
                maxval = DACE::cons(xval.at(k, i));

        // cout << maxval<<endl;
        trustRegion << maxval << "\t";
        trustRegion << std::endl;
    }

    trustRegion.close();
}
