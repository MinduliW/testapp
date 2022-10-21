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

    // cout << params.tf << endl;
    typedef DACE::DA state_type; // define state type as DA
    DACE::DA::init(2, 9);        // initialise DACE for X order in N variables

    DACE::AlgebraicVector<state_type> x0(6), fx(6), dv(3);
    DACE::AlgebraicMatrix<state_type> Jacobian(6, 9), xval(params.N, 9);
    DACE::DA derivative(1), result(1), eclipse, coef;
    DACE::AlgebraicMatrix<double> STMgamma(6, 9), Jbar(6, 9), b(6, 9);

    // read and store the guess data.
    double statesGuess[params.N + 1][9], xarray[9], constant[6],  B, normJbar,time,maxval;

    vector<double> zerostate(9);

    vector<state_type> zerostate2(9);

    fstream myfile;
    int i = 0;

    myfile.open("iter.txt", ios::in); // open a file to perform read operation using file object
    if (myfile.is_open())
    {
        // checking whether the file is open
        string tp;
        while (getline(myfile, tp))
        { // read data from file object and put it into string
            removeDupWord(tp, xarray);
            for (unsigned int j = 0; j < 9; j++)
            {
                // cout << "J=" << j << endl;
                statesGuess[i][j] = xarray[j];
            }
            i++;
            // cout <<  i << endl;
        }
        myfile.close(); // close the file object.
    }

    fstream outputfile,nonlinIndex;
    //DACE::AlgebraicMatrix<double> b(6, 9);
    outputfile.open("outputcpp.txt", fstream::out);
    nonlinIndex.open("nonlinIndex.txt", fstream::out);
    
     for (int k = 0; k < params.N; k++)
    // for (int k = 0; k < 1; k++)
    {
      // update time. 
      time += params.dt;
        
        for (unsigned int i = 0; i < 6; i++)
        {
            x0[i] = statesGuess[k][i] + DACE::DA(i + 1);
            
            if (i < 3)
            {
                // params.u0[i] = statesGuess[k][i + 6];
                cst[i] = statesGuess[k][i + 6];
            }
        }
      
    //   cout << x0 << endl;
    //   cout << getEclipse(500.0, x0, 1000.0 , 1.0,params.LU, params.Rs, params.Re) << endl;
   
        if (cst[8] < 2)
        {

            fx = RK78(6, x0, 0.0, params.dt, cst, DynamicsCart);
        }
        else
        {
           //eclipse = getEclipse(time-0.5*params.dt,  x0, 100.0 , 0.999,params.LU, params.Rs, params.Re);
           //cout << "nu eclipse = " << eclipse<< endl;
                 
            for (unsigned int i = 0; i < 3; i++)
            {
                // cout << params.u0[i] << endl;

                dv[i] = cst[i] + DACE::DA(i + 7); //
                // cout << dv<< endl;
                if (params.eclipses >0){
                    // cout << dv<< endl;
                   // dv[i] =  dv[i]*(1.0000001-eclipse);
                     
               }
                x0[i + 3] = x0[i + 3] + dv[i];
            }
           //cout << dv<< endl;
           

            fx = RK78(6, x0, 0.0, params.dt, cst, DynamicsCart);
        }

        // convert osculating to mean.

        STMgamma = fx.linear();


        for (unsigned int i = 0; i < 6; i++)
        {
            constant[i] = DACE::cons(fx[i]);
            // cout << constant[i]<< endl;
        }

        // std::cout << "STM and gamma" << STMgamma << endl;

        for (unsigned int j = 0; j < 6; j++) // Prints row of x
        {
            for (unsigned int i = 0; i < 9; i++)
            {
                outputfile << STMgamma.at(j, i) << "\t";
            }

            outputfile << constant[j] << "\t";

            outputfile << std::endl;
        }


         normJbar = 0;
         DACE:: DA B;
        // // calculate the jacobian
        // for (unsigned int i = 0; i < 6; i++)
        // {
        //     for (unsigned int j = 1; j < 10; j++)
        //     {
        //         derivative = fx[i].deriv(j);
        //         Jacobian.at(i, j - 1) = derivative;
        //         B += (Jacobian.at(i, j - 1));

        //         Jbar.at(i, j - 1) = derivative.eval(zerostate);

        //         normJbar += abs(Jbar.at(i, j-1));
        //         //b.at(i, j - 1);
        //     }
           
        // }

          // calculate the jacobian
        for (unsigned int i = 0; i < 6; i++)
        {
            for (unsigned int j = 1; j < 10; j++)
            {
                derivative = fx[i].deriv(j);
                Jacobian.at(i, j - 1) = derivative;
                //B += (Jacobian.at(i, j - 1));

                Jbar.at(i, j - 1) = derivative.eval(zerostate);

                normJbar += abs(Jbar.at(i, j-1));
                //b.at(i, j - 1);
            }
           
        }
        
        // analyse the dependency on each variable individually
  
        for(int l = 0; l <9; l++){
            zerostate2 = {0.0, 0.0, 0.0, 0.0,0.0, 0.0,0.0, 0.0,0.0};
            zerostate2[l] = DA(l+1);
            coef = 0;
            for (unsigned int i = 0; i < 6; i++){
                for (unsigned int j = 1; j < 10; j++)
                {
                    derivative = Jacobian.at(i, j - 1);

                     // get the absolute of the first order coefficient 
                    coef +=   abs(derivative.eval(zerostate2).deriv(l+1));

                }
                }

                cout << DACE::cons(coef) <<endl;

                  //nu = coef/normJbar;
                  xval.at(k,l) = params.nu*normJbar/coef;

                 

                }

       

    }

    // find the max variation in each element to provide the required nu.
    for (unsigned int i =0; i< 9; i++){
        maxval = 0.0;
        for (unsigned int k = 0; k < params.N; k++){
            if (DACE::cons(xval.at(k,i))> (maxval)){
                maxval = DACE::cons(xval.at(k,i));
            }

        }
        
        nonlinIndex << maxval << "\t";
          
        nonlinIndex << std::endl;
    }

    outputfile.close();
    nonlinIndex.close();
     //cout << "xval" << xval<<endl;
}

       //    cout << B << endl;
        //     cout << nu << endl;
        //      cout << normJbar << endl;
     
        //B =0;
       
        

        // DACE::AlgebraicMatrix<double>  b(6, 9);
        // for (unsigned int i = 0; i < 6; i++)
        // {
        //     for (unsigned int j = 1; j < 10; j++)
        //     {
        //         result = Jacobian.at(i, j - 1); 
        //        // cout << "here" << endl;
        //         //cout << result << endl;
        //        // b.at(i, j - 1) = 0;
        //         // for (unsigned int p = 0; p < 9; p++)
        //         // {
        //         //     // cout << result.deriv(p + 1) << endl;
        //         //     // cout <<abs(DACE::cons(result.deriv(p + 1))) << endl;
        //         //     //    cout <<"here" << endl;
        //         //     b.at(i, j - 1) += abs(DACE::cons(result.deriv(p + 1)));
        //         // }
                
        //         B += Jacobian.at(i, j - 1); //b.at(i, j - 1);
        //         normJbar += (Jbar.at(i, j-1));
        //     }
        // }
     

        // sum up b to get B
       //cout << "nonlinearity index = " << nu << endl;
        