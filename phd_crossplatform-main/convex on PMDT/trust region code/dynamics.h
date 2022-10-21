#ifndef dynamics_H
#define dynamics_H


template <typename U>
U Smoothed_eclipse(U cs, U ct, DACE::AlgebraicVector<double> x_sun, DACE::AlgebraicVector<double> x_sat, U Rs , U Re)
{

  using DACE::cos;
  using DACE::pow;
  using DACE::sin;
  using DACE::sqrt;
  using std::cos;
  using std::pow;
  using std::sin;
  using std::sqrt;
  
  //cout << "Rs =" << Rs << endl;
 

  U a_SR = asin(Rs / x_sun.vnorm());
  U a_BR = asin(Re / x_sat.vnorm());

  U a_D = acos(x_sat.dot(x_sun) / x_sat.vnorm() / x_sun.vnorm());

  U gamma_L = 1 / (1 + exp(-cs * (a_D - ct * (a_SR + a_BR))));

  return gamma_L;
}

template <typename U>
DACE::AlgebraicVector<U> getSunPosVec(U JD, U LU)
{

  U pi = 4 * atan(1.0);

  U T_UT1 = (JD - 2451545.0) / 36525;

  U lambdaMSun = 280.460 + 36000.771 * T_UT1;

  U MSun = 357.5291092 + 35999.05034 * T_UT1;

  U lambdaEcliptic = lambdaMSun + 1.914666471 * sin(MSun * pi / 180.0) + 0.019994643 * sin(2 * MSun * pi / 180.0);

  U r_Sun = 1.000140612 - 0.016708617 * cos(MSun * pi / 180.0) - 0.000139589 * cos(2 * MSun * pi / 180.0);

  U epsilon = 23.439291 - 0.0130042 * T_UT1;

  DACE::AlgebraicVector<U> SunVec(3);
  SunVec[0] = r_Sun * cos(lambdaEcliptic * pi / 180.0);
  SunVec[1] = r_Sun * cos(epsilon * pi / 180.0) * sin(lambdaEcliptic * pi / 180.0);
  SunVec[2] = r_Sun * sin(epsilon * pi / 180.0) * sin(lambdaEcliptic * pi / 180.0);

  SunVec = 1.495978707e11 * SunVec / LU;
  
  //cout << SunVec << endl;
  return SunVec;
}

template <typename T, typename U>
double getEclipse(U t, DACE::AlgebraicVector<T> x, U cs , U ct,U LU, U Rs, U Re)
{

  DACE::AlgebraicVector<double> r_Earth2Sun(3);

  r_Earth2Sun = getSunPosVec(t + 2400000.5 + 60000, LU);

  // convert dace to double.
  DACE::AlgebraicVector<double> x_double(6);
  for (unsigned int i = 0; i < 6; i++)
  {
    x_double[i] = DACE::cons(x[i]);
  }

  DACE::AlgebraicVector<double> CART = x_double;

  DACE::AlgebraicVector<double> rSat2Earth(3);
  DACE::AlgebraicVector<double> rSat2Sun(3);

  for (unsigned int i = 0; i < 3; i++)
  {
    rSat2Earth[i] = -CART[i];
    rSat2Sun[i] = rSat2Earth[i] + r_Earth2Sun[i];
  }

  double numinus1 = Smoothed_eclipse(cs, ct, rSat2Sun, rSat2Earth, Rs , Re);

  return 1.0 - numinus1;
}

template <typename T, typename U>
DACE::AlgebraicVector<T> DynamicsCart(DACE::AlgebraicVector<T> x, U t, const double *cst)
{

  using DACE::cos;
  using DACE::pow;
  using DACE::sin;
  using DACE::sqrt;
  using std::cos;
  using std::pow;
  using std::sin;
  using std::sqrt;

  const double mu = cst[5];
  const double J2 = cst[6];
  const double Re = cst[7];

  DACE::AlgebraicVector<T> dxdt(6);

  T r2 = x[0] * x[0] + x[1] * x[1] + x[2] * x[2];
  // cout << DACE::cons(term) << endl;

  if (DACE::cons(r2) <= 0)
  {
    cout << "So the iter.tex has a mistake, its giving all states as zero, likely after infeasible run" << endl;
  }

  T r = sqrt(r2);

  T r3 = r * r * r;

  for (unsigned int i = 0; i < 3; i++)
  {
    dxdt[i] = x[3 + i];
    dxdt[i + 3] = -mu * x[i] / (r3);
  }

  //add J2
  //if (params.addJ2 == true){
  dxdt[3] = dxdt[3] - mu * x[0] / (r3)*1.5 * J2 * (Re / r) * (Re / r) * (1 - 5 * x[2] * x[2] / r2);
  dxdt[4] = dxdt[4] - mu * x[1] / (r3)*1.5 * J2 * (Re / r) * (Re / r) * (1 - 5 * x[2] * x[2] / r2);
  dxdt[5] = dxdt[5] - mu * x[2] / (r3)*1.5 * J2 * (Re / r) * (Re / r) * (3 - 5 * x[2] * x[2] / r2);

  //}

  if (cst[8] <2)
  {
    DACE::AlgebraicVector<T> u(3);
    for (unsigned int i = 0; i < 3; i++)
    {
      // cout << params.u0[i] << endl;
      u[i] = cst[i] + DACE::DA(i + 7); //
    }

    dxdt[3] += u[0];
    dxdt[4] += u[1];
    dxdt[5] += u[2];

    //cout << u[2] << endl;
  }

  return dxdt;
}

#endif