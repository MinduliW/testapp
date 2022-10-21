#ifndef osctomean_H
#define osctomean_H


#include <dace/dace.h>
#include <cmath>
#include <ctime>
#include <fstream>
using namespace std;
using namespace DACE;


template<typename T> T true2eccAnomaly(const T theta, const T e)
{
	return 2.0 * atan2(sqrt(1. - e)*sin(theta / 2.), sqrt(1. + e) * cos(theta / 2.));
}

template<typename T> T ecc2trueAnomaly(const T E, const T e)
{
	return 2.0 * atan2(sqrt(1. + e)*sin(E / 2.), sqrt(1. - e) * cos(E / 2.));
}

template<typename T> T mean2eccAnomaly(const T M, const T e)
{
	T E = M;

	for (int i = 0; i < 20; i++) {
		E = M + e*sin(E);
	}
	return E;
}

template<typename T> T mean2trueAnomaly(const T M, const T e)
{
	T E = mean2eccAnomaly(M, e);

	return ecc2trueAnomaly(E, e);
}


template<typename T> T true2meanAnomaly(const T theta, const T e)
{
	T E = true2eccAnomaly(theta, e);

	return E - e*sin(E);
}

template<typename T> AlgebraicVector<T> cart2kep(const AlgebraicVector<T>& rv, const double mu = 398600.4418)
{
	//const double mu = 398600.4415;
	AlgebraicVector<T> kep(6);

	AlgebraicVector<T> rr(3), vv(3);
	for (int i = 0; i < 3; i++)
	{
		rr[i] = rv[i];
		vv[i] = rv[i + 3];
	}

	T r = rr.vnorm();
	T v = vv.vnorm();
	AlgebraicVector<T> h = cross(rr, vv);
	//cout << h << endl;

	kep[0] = mu / (2.0 * (mu / r - pow(v, 2) / 2.0));

	T h1sqr = pow(h[0], 2);
	T h2sqr = pow(h[1], 2);

	T RAAN;
	if (cons(h1sqr + h2sqr) == 0.0)
	{
		RAAN = 0.0;
	}
	else
	{
		T sinOMEGA = h[0] / sqrt(h1sqr + h2sqr);
		T cosOMEGA = -1.0*h[1] / sqrt(h1sqr + h2sqr);
		if (cons(cosOMEGA) >= 0.0)
		{
			if (cons(sinOMEGA) >= 0.0)
			{
				RAAN = asin(h[0] / sqrt(h1sqr + h2sqr));
			}
			else
			{
				RAAN = 2.0 * M_PI + asin(h[0] / sqrt(h1sqr + h2sqr));
			}
		}
		else
		{
			if (cons(sinOMEGA) >= 0.0)
			{
				RAAN = acos(-1.0*h[1] / sqrt(h1sqr + h2sqr));
			}
			else
			{
				RAAN = 2.0 * M_PI - acos(-1.0*h[1] / sqrt(h1sqr + h2sqr));
			}
		}
	}

	//RAAN = real(RAAN);

	AlgebraicVector<T> ee = 1.0 / mu*(cross(vv, h)) - rr / vnorm(rr);
	T e = vnorm(ee);

	T i = acos(h[2] / vnorm(h));

	kep[1] = e;
	kep[2] = i;
	kep[3] = RAAN;

	T omega;
	T theta;
	if (cons(e) <= 1.0e-8 && cons(i) < 1.0e-8)
	{
		e = 0.0;
		omega = atan2(rr[1], rr[0]);
		theta = 0.0;
		kep[4] = omega;
		kep[5] = theta;
		return kep;
	}

	if (cons(e) <= 1.0e-8 && cons(i) >= 1.0e-8)
	{
		omega = 0;
		AlgebraicVector<T> P(3), Q(3), W(3);
		P[0] = cos(omega)*cos(RAAN) - sin(omega)*sin(i)*sin(RAAN);
		P[1] = -1.0*sin(omega)*cos(RAAN) - cos(omega)*cos(i)*sin(RAAN);
		P[2] = sin(RAAN)*sin(i);
		Q[0] = cos(omega)*sin(RAAN) + sin(omega)*cos(i)*cos(RAAN);
		Q[1] = -1.0*sin(omega)*sin(RAAN) + cos(omega)*cos(i)*cos(RAAN);
		Q[2] = -1.0*cos(RAAN)*sin(i);
		W[0] = sin(omega)*sin(i);
		W[1] = cos(omega)*sin(i);
		W[2] = cos(i);
		AlgebraicVector<T> rrt = P*rr[0] + Q*rr[1] + W*rr[2];
		theta = atan2(rrt[1], rrt[0]);
		kep[4] = omega;
		kep[5] = theta;
		return kep;
	}

	T dotRxE = dot(rr, ee);
	T RxE = vnorm(rr)*vnorm(ee);
	if (abs(cons(dotRxE)) > abs(cons(RxE)) && abs(cons(dotRxE)) - abs(cons(RxE)) < abs(numeric_limits<double>::epsilon()*cons(dotRxE)))
	{
		dotRxE -= numeric_limits<double>::epsilon()*dotRxE;
	}
	theta = acos(dotRxE / RxE);

	if (cons(dot(rr, vv)) < 0.0)
	{
		theta = 2.0 * M_PI - theta;
	}

	if (cons(i) <= 1.0e-8 && cons(e) >= 1.0e-8)
	{
		i = 0.0;
		omega = atan2(ee[1], ee[0]);
		kep[4] = omega;
		kep[5] = theta;
		return kep;
	}

	T sino = rr[2] / r / sin(i);
	T coso = (rr[0] * cos(RAAN) + rr[1] * sin(RAAN)) / r;
	T argLat;
	if (cons(coso) >= 0.0)
	{
		if (cons(sino) >= 0.0)
		{
			argLat = asin(rr[2] / r / sin(i));
		}
		else
		{
			argLat = 2.0 * M_PI + asin(rr[2] / r / sin(i));
		}
	}
	else
	{
		if (cons(sino) >= 0.0)
		{
			argLat = acos((rr[0] * cos(RAAN) + rr[1] * sin(RAAN)) / r);
		}
		else
		{
			argLat = 2.0 * M_PI - acos((rr[0] * cos(RAAN) + rr[1] * sin(RAAN)) / r);
		}
	}
	//argLat = real(argLat);
	omega = argLat - theta;

	if (cons(omega) < 0.0)
	{
		omega = omega + 2.0 * M_PI;
	}
	//omega = real(omega);

	kep[4] = omega;
	kep[5] = theta;

	return kep;
}

template<typename T> AlgebraicVector<T> kep2delaunay(const AlgebraicVector<T>& kep, const double mu)
{
	// Keplerian orbital elements
	T a = kep[0];
	T e = kep[1];
	T i = kep[2];
	T RAAN = kep[3];
	T omega = kep[4];
	T M = kep[5];

	// Delaunay elements
	AlgebraicVector<T> delaunay(6);
	delaunay[0] = M; // l
	delaunay[1] = omega; // g
	delaunay[2] = RAAN; // h
	delaunay[3] = sqrt(mu*a); // L
	delaunay[4] = sqrt(1 - pow(e, 2)) * delaunay[3]; // G
	delaunay[5] = cos(i)*delaunay[4]; // H

	return delaunay;
}



template<typename T> AlgebraicVector<T> delaunay2kep(const AlgebraicVector<T>& delaunay, const double mu)
{
	// Delaunay elements
	T l = delaunay[0]; // l
	T g = delaunay[1]; // g
	T h = delaunay[2]; // h
	T L = delaunay[3]; // L
	T G = delaunay[4]; // G
	T H = delaunay[5]; // H

	// Keplerian orbital elements
	AlgebraicVector<T> kep(6);
	kep[0] = pow(L, 2) / mu; // a
	kep[1] = sqrt(1 - pow((G / L), 2)); // e
	kep[2] = acos(H / G); // i
	kep[3] = h; // RAAN
	kep[4] = g; // omega
	kep[5] = l; // M

	return kep;
}


template<typename T> AlgebraicVector<T> kep2hill(const AlgebraicVector<T>& kep, const double mu)

{
    AlgebraicVector<T> hill(6);
    T u, f, p;
    
    p = kep[0]*(1.0 - kep[1]*kep[1]);
    f = kep[5];
    
    hill[4] = sqrt(mu*p);
    hill[0] = p/(1.0 + kep[1]*cos(f));
    hill[1] = f + kep[4];
    hill[2] = kep[3];
    hill[3] = (hill[4]/p)*kep[1]*sin(f);
    hill[5] = hill[4]*cos(kep[2]);
    
    return hill;
}

template<typename T> AlgebraicVector<T> hill2cart(const AlgebraicVector<T>& hill, const double mu)
//  hill[] = {r, th, nu, R, Th, N}
//  cart[] = {x, y, z, vx, vy, vz}
{
    AlgebraicVector<T> u(3);
    AlgebraicVector<T> cart(6);
    
    T r, th, nu, R, Th, ci, si;
    int i;
    
    r  = hill[0];
    th = hill[1];
    nu = hill[2];
    R  = hill[3];
    Th = hill[4];
    ci = hill[5]/hill[4];
    si = sqrt(1.0 - ci*ci);
    
    u[0] = cos(th)*cos(nu) - ci*sin(th)*sin(nu);
    u[1] = cos(th)*sin(nu) + ci*sin(th)*cos(nu);
    u[2] = si*sin(th);
    
    for(i = 0; i < 3; ++i)
        
        cart[i] = r*u[i];
    cart[3] = (R*cos(th) - Th*sin(th)/r)*cos(nu) - (R*sin(th) + Th*cos(th)/r)*sin(nu)*ci;
    cart[4] = (R*cos(th) - Th*sin(th)/r)*sin(nu) + (R*sin(th) + Th*cos(th)/r)*cos(nu)*ci;
    cart[5] = (R*sin(th) + Th*cos(th)/r)*si;
    
    return cart;
}


template<typename T> AlgebraicVector<T> hill2kep(const AlgebraicVector<T>& hill, const double mu)
//  hill[] = {r, th, nu, R, Th, Nu}
//  cart[] = {x, y, z, vx, vy, vz}
{
    AlgebraicVector<T> kep(6);
    
    T r  = hill[0];
    T th = hill[1];
    T nu = hill[2];
    T R  = hill[3];
    T Th = hill[4];
    T Nu = hill[5];
    
    T i = acos(Nu/Th);
    T cs  =   (-1.0 + pow(Th,2)/(mu*r))*cos(th) + (R*Th*sin(th))/mu;
    T ss  =  -((R*Th*cos(th))/mu) + (-1.0 + pow(Th,2)/(mu*r))*sin(th);
    T e = sqrt(cs*cs+ss*ss);
    T p = Th*Th/mu;
    T costrue = 1.0/e*(p/r-1.0);
    T f = acos(costrue);
    
    if (cons(R)<0.0) {
        f = 2.0*M_PI-f;
    }
    T a = p/(1-e*e);
    
    kep[0] = a;
    kep[1] = e;
    kep[2] = i;
    kep[3] = nu;
    kep[4] = th-f;
    kep[5] = f;
    
    return kep;
    
}


template<typename T> AlgebraicVector<T> osculating2meanHill(AlgebraicVector<T> hillOsc, double mu, double J2, double rE, T cont)
{
    
    // Mean Delaunay elements
    T r   = hillOsc[0]; // l
    T th  = hillOsc[1]; // g
    T nu  = hillOsc[2]; // h
    T R   = hillOsc[3]; // L
    T Th  = hillOsc[4]; // G
    T Nu  = hillOsc[5]; // H
    
    T ci = Nu/Th;
    T si = sqrt(1.0-ci*ci);
    T cs  =   (-1.0 + pow(Th,2)/(mu*r))*cos(th) + (R*Th*sin(th))/mu;
    T ss  =  -((R*Th*cos(th))/mu) + (-1.0 + pow(Th,2)/(mu*r))*sin(th);
    T e = sqrt(cs*cs+ss*ss);
    T eta  = sqrt(1.0-e*e);
    
    T beta = 1.0/(1.0+eta);
    
    T p = Th*Th/mu;
    T costrue = 1/e*(p/r-1);
    
    T f = acos(costrue);
    
    if (cons(R)<0.0) {
        f = 2.0*M_PI-f;
    }
    
    T M = true2meanAnomaly(f,e);
    
    T phi  = f - M;
    
    r = r +  (cont)*((pow(rE,2)*beta*J2)/(2.*r) - (3*pow(rE,2)*beta*J2*pow(si,2))/(4.*r) +
                     (pow(rE,2)*eta*J2*pow(mu,2)*r)/pow(Th,4) - (3*pow(rE,2)*eta*J2*pow(mu,2)*r*pow(si,2))/(2.*pow(Th,4)) +
                     (pow(rE,2)*J2*mu)/(2.*pow(Th,2)) - (pow(rE,2)*beta*J2*mu)/(2.*pow(Th,2)) -
                     (3.*pow(rE,2)*J2*mu*pow(si,2))/(4.*pow(Th,2)) + (3*pow(rE,2)*beta*J2*mu*pow(si,2))/(4.*pow(Th,2)) -
                     (pow(rE,2)*J2*mu*pow(si,2)*cos(2*th))/(4.*pow(Th,2)));
    
    
    th = th + (cont)*((-3.*pow(rE,2)*J2*pow(mu,2)*phi)/pow(Th,4) + (15.*pow(rE,2)*J2*pow(mu,2)*phi*pow(si,2))/(4.*pow(Th,4)) -
                      (5.*pow(rE,2)*J2*mu*R)/(2.*pow(Th,3)) - (pow(rE,2)*beta*J2*mu*R)/(2.*pow(Th,3)) +
                      (3.*pow(rE,2)*J2*mu*R*pow(si,2))/pow(Th,3) + (3.*pow(rE,2)*beta*J2*mu*R*pow(si,2))/(4.*pow(Th,3)) -
                      (pow(rE,2)*beta*J2*R)/(2.*r*Th) + (3.*pow(rE,2)*beta*J2*R*pow(si,2))/(4.*r*Th) +
                      (-(pow(rE,2)*J2*mu*R)/(2.*pow(Th,3)) + (pow(rE,2)*J2*mu*R*pow(si,2))/pow(Th,3))*cos(2.*th) +
                      (-(pow(rE,2)*J2*pow(mu,2))/(4.*pow(Th,4)) + (5.*pow(rE,2)*J2*pow(mu,2)*pow(si,2))/(8.*pow(Th,4)) +
                       (pow(rE,2)*J2*mu)/(r*pow(Th,2)) - (3.*pow(rE,2)*J2*mu*pow(si,2))/(2.*r*pow(Th,2)))*sin(2.*th));
    
    nu = nu + (cont)*((3.*pow(rE,2)*ci*J2*pow(mu,2)*phi)/(2.*pow(Th,4)) + (3.*pow(rE,2)*ci*J2*mu*R)/(2.*pow(Th,3)) +
                      (pow(rE,2)*ci*J2*mu*R*cos(2.*th))/(2.*pow(Th,3)) +
                      ((pow(rE,2)*ci*J2*pow(mu,2))/(4.*pow(Th,4)) - (pow(rE,2)*ci*J2*mu)/(r*pow(Th,2)))*sin(2.*th));
    
    
    R = R  + (cont)*(-(pow(rE,2)*beta*J2*R)/(2.*pow(r,2)) + (3.*pow(rE,2)*beta*J2*R*pow(si,2))/(4.*pow(r,2)) -
                     (pow(rE,2)*eta*J2*pow(mu,2)*R)/(2.*pow(Th,4)) + (3.*pow(rE,2)*eta*J2*pow(mu,2)*R*pow(si,2))/(4.*pow(Th,4)) +
                     (pow(rE,2)*J2*mu*pow(si,2)*sin(2.*th))/(2.*pow(r,2)*Th));
    
    
    Th = Th  + (cont)*(((pow(rE,2)*J2*pow(mu,2)*pow(si,2))/(4.*pow(Th,3)) - (pow(rE,2)*J2*mu*pow(si,2))/(r*Th))*cos(2.*th) -
                       (pow(rE,2)*J2*mu*R*pow(si,2)*sin(2.*th))/(2.*pow(Th,2)));
    
    Nu = Nu +  0;
    
    AlgebraicVector<T> hillMean(6);
    
    
    hillMean[0] = r;
    hillMean[1] = th;
    hillMean[2] = nu;
    hillMean[3] = R;
    hillMean[4] = Th;
    hillMean[5] = Nu;
    
    return hillMean;
}


template<typename T> AlgebraicVector<T> mean2osculatingHill(AlgebraicVector<T> hillMean, double mu, double J2, double rE, T cont)
{
    
    // Mean Delaunay elements
    T r  = hillMean[0]; // l
    T th = hillMean[1]; // g
    T nu = hillMean[2]; // h
    T R  = hillMean[3]; // L
    T Th = hillMean[4]; // G
    T Nu = hillMean[5]; // H
    
    T ci = Nu/Th;
    T si = sqrt(1.0-ci*ci);
    T cs  =   (-1.0 + pow(Th,2)/(mu*r))*cos(th) + (R*Th*sin(th))/mu;
    T ss  =  -((R*Th*cos(th))/mu) + (-1.0 + pow(Th,2)/(mu*r))*sin(th);
    T e = sqrt(cs*cs+ss*ss);
    T eta  = sqrt(1.0-e*e);
    T beta = 1.0/(1.0+eta);
    
    T p = Th*Th/mu;
    T costrue = 1/e*(p/r-1);
    
    T f = acos(costrue);
    
    if (cons(R)<0.0) {
        f = 2.0*M_PI-f;
    }
    
    T M = true2meanAnomaly(f,e);
    
    T phi  = f - M;
    
    r = r -  (cont)*((pow(rE,2)*beta*J2)/(2.*r) - (3*pow(rE,2)*beta*J2*pow(si,2))/(4.*r) +
                     (pow(rE,2)*eta*J2*pow(mu,2)*r)/pow(Th,4) - (3*pow(rE,2)*eta*J2*pow(mu,2)*r*pow(si,2))/(2.*pow(Th,4)) +
                     (pow(rE,2)*J2*mu)/(2.*pow(Th,2)) - (pow(rE,2)*beta*J2*mu)/(2.*pow(Th,2)) -
                     (3.*pow(rE,2)*J2*mu*pow(si,2))/(4.*pow(Th,2)) + (3*pow(rE,2)*beta*J2*mu*pow(si,2))/(4.*pow(Th,2)) -
                     (pow(rE,2)*J2*mu*pow(si,2)*cos(2*th))/(4.*pow(Th,2)));
    
    
    th = th - (cont)*((-3.*pow(rE,2)*J2*pow(mu,2)*phi)/pow(Th,4) + (15.*pow(rE,2)*J2*pow(mu,2)*phi*pow(si,2))/(4.*pow(Th,4)) -
                      (5.*pow(rE,2)*J2*mu*R)/(2.*pow(Th,3)) - (pow(rE,2)*beta*J2*mu*R)/(2.*pow(Th,3)) +
                      (3.*pow(rE,2)*J2*mu*R*pow(si,2))/pow(Th,3) + (3.*pow(rE,2)*beta*J2*mu*R*pow(si,2))/(4.*pow(Th,3)) -
                      (pow(rE,2)*beta*J2*R)/(2.*r*Th) + (3.*pow(rE,2)*beta*J2*R*pow(si,2))/(4.*r*Th) +
                      (-(pow(rE,2)*J2*mu*R)/(2.*pow(Th,3)) + (pow(rE,2)*J2*mu*R*pow(si,2))/pow(Th,3))*cos(2.*th) +
                      (-(pow(rE,2)*J2*pow(mu,2))/(4.*pow(Th,4)) + (5.*pow(rE,2)*J2*pow(mu,2)*pow(si,2))/(8.*pow(Th,4)) +
                       (pow(rE,2)*J2*mu)/(r*pow(Th,2)) - (3.*pow(rE,2)*J2*mu*pow(si,2))/(2.*r*pow(Th,2)))*sin(2.*th));
    
    nu = nu - (cont)*((3.*pow(rE,2)*ci*J2*pow(mu,2)*phi)/(2.*pow(Th,4)) + (3.*pow(rE,2)*ci*J2*mu*R)/(2.*pow(Th,3)) +
                      (pow(rE,2)*ci*J2*mu*R*cos(2.*th))/(2.*pow(Th,3)) +
                      ((pow(rE,2)*ci*J2*pow(mu,2))/(4.*pow(Th,4)) - (pow(rE,2)*ci*J2*mu)/(r*pow(Th,2)))*sin(2.*th));
    
    
    R = R  - (cont)*(-(pow(rE,2)*beta*J2*R)/(2.*pow(r,2)) + (3.*pow(rE,2)*beta*J2*R*pow(si,2))/(4.*pow(r,2)) -
                     (pow(rE,2)*eta*J2*pow(mu,2)*R)/(2.*pow(Th,4)) + (3.*pow(rE,2)*eta*J2*pow(mu,2)*R*pow(si,2))/(4.*pow(Th,4)) +
                     (pow(rE,2)*J2*mu*pow(si,2)*sin(2.*th))/(2.*pow(r,2)*Th));
    
    
    Th = Th  - (cont)*(((pow(rE,2)*J2*pow(mu,2)*pow(si,2))/(4.*pow(Th,3)) - (pow(rE,2)*J2*mu*pow(si,2))/(r*Th))*cos(2.*th) -
                       (pow(rE,2)*J2*mu*R*pow(si,2)*sin(2.*th))/(2.*pow(Th,2)));
    
    Nu = Nu +  0;
    
    AlgebraicVector<T> hillOsc(6);
    
    hillOsc[0] = r;
    hillOsc[1] = th;
    hillOsc[2] = nu;
    hillOsc[3] = R;
    hillOsc[4] = Th;
    hillOsc[5] = Nu;
    
    return hillOsc;
}

template<typename T> AlgebraicVector<T> osculating2mean(AlgebraicVector<T> delaunayOsc, double mu, double J2, double rE, T cont)
{
	// Osculating Delaunay elements
	T l = delaunayOsc[0]; // l
	T g = delaunayOsc[1]; // g
	T h = delaunayOsc[2]; // h
	T L = delaunayOsc[3]; // L
	T G = delaunayOsc[4]; // G
	T H = delaunayOsc[5]; // H

	// Osculating Keplerian elements
	AlgebraicVector<T> kep = delaunay2kep(delaunayOsc, mu);
	T a = kep[0];
	T e = kep[1];
	T incl = kep[2];

	T f = mean2trueAnomaly(l, e);	// true anomaly

	T eta = G / L;	// sqrt(1-e^2)
	T ci = cos(incl);	// cos(i)
	T si = sin(incl);	// sin(i)

	T lmod = modulusLocal(l, 2. * M_PI);
	T phi = modulusLocal(f - lmod, 2.0 * M_PI);	// true minus mean anomaly

	T p = a*(1.0 - pow(e, 2));

	T epsilon = 1.0;
	T r = p / (1.0 + e*cos(f));	// Radial distance

	T n = sqrt(mu / pow(a, 3));	// mean motion

	T R = e*G*sin(f) / (a*pow(eta, 2));


	T lmean = 
		l + (cont)*( (J2*R*pow(rE, 2)*epsilon) / (2.*pow(e, 2)*L*r) - (3.*J2*R*pow(si, 2)*pow(rE, 2)*epsilon) / (4.*pow(e, 2)*L*r) +
		(J2*R*pow(rE, 2)*epsilon*mu) / (2.*pow(e, 2)*pow(L, 3)*pow(eta, 2)) -
		(3.*J2*R*pow(si, 2)*pow(rE, 2)*epsilon*mu) / (4.*pow(e, 2)*pow(L, 3)*pow(eta, 2)) + (J2*R*pow(rE, 2)*epsilon*cos(f)) / (2.*e*L*r) -
		(3.*J2*R*pow(si, 2)*pow(rE, 2)*epsilon*cos(f)) / (4.*e*L*r) + (J2*R*pow(rE, 2)*epsilon*mu*cos(f)) / (2.*e*pow(L, 3)*pow(eta, 2)) -
		(3.*J2*R*pow(si, 2)*pow(rE, 2)*epsilon*mu*cos(f)) / (4.*e*pow(L, 3)*pow(eta, 2)) +
		(3.*J2*R*pow(si, 2)*pow(rE, 2)*epsilon*cos(f + 2. * g)) / (8.*e*L*r) +
		(3.*J2*R*pow(si, 2)*pow(rE, 2)*epsilon*mu*cos(f + 2. * g)) / (8.*e*pow(L, 3)*pow(eta, 2)) +
		(3.*J2*R*pow(si, 2)*pow(rE, 2)*epsilon*cos(2. * f + 2. * g)) / (4.*pow(e, 2)*L*r) +
		(3.*J2*R*pow(si, 2)*pow(rE, 2)*epsilon*mu*cos(2. * f + 2. * g)) / (4.*pow(e, 2)*pow(L, 3)*pow(eta, 2)) +
		(3.*J2*R*pow(si, 2)*pow(rE, 2)*epsilon*cos(3.*f + 2. * g)) / (8.*e*L*r) +
		(3.*J2*R*pow(si, 2)*pow(rE, 2)*epsilon*mu*cos(3.*f + 2. * g)) / (8.*e*pow(L, 3)*pow(eta, 2)) +
		(J2*pow(rE, 2)*epsilon*pow(mu, 2)*sin(f)) / (2.*e*pow(L, 4)*eta) -
		(3.*J2*pow(si, 2)*pow(rE, 2)*epsilon*pow(mu, 2)*sin(f)) / (4.*e*pow(L, 4)*eta) +
		(3.*J2*pow(si, 2)*pow(rE, 2)*epsilon*pow(mu, 2)*sin(f + 2. * g)) / (8.*e*pow(L, 4)*eta) +
		(J2*pow(si, 2)*pow(rE, 2)*epsilon*pow(mu, 2)*sin(3.*f + 2. * g)) / (8.*e*pow(L, 4)*eta));

	T gmean = 
		g + (cont)*(- (pow(G, 2)*J2*R*pow(rE, 2)*epsilon) / (2.*pow(e, 2)*pow(L, 3)*r*pow(eta, 3)) +
		(3.*pow(G, 2)*J2*R*pow(si, 2)*pow(rE, 2)*epsilon) / (4.*pow(e, 2)*pow(L, 3)*r*pow(eta, 3)) -
		(J2*R*pow(rE, 2)*epsilon*mu) / (2.*pow(e, 2)*pow(L, 3)*pow(eta, 3)) +
		(3.*J2*R*pow(si, 2)*pow(rE, 2)*epsilon*mu) / (4.*pow(e, 2)*pow(L, 3)*pow(eta, 3)) -
		(3.*J2*pow(rE, 2)*epsilon*pow(mu, 2)*phi) / (pow(L, 4)*pow(eta, 4)) +
		(15. * J2*pow(si, 2)*pow(rE, 2)*epsilon*pow(mu, 2)*phi) / (4.*pow(L, 4)*pow(eta, 4)) -
		(pow(G, 2)*J2*R*pow(rE, 2)*epsilon*cos(f)) / (2.*e*pow(L, 3)*r*pow(eta, 3)) +
		(3.*pow(G, 2)*J2*R*pow(si, 2)*pow(rE, 2)*epsilon*cos(f)) / (4.*e*pow(L, 3)*r*pow(eta, 3)) -
		(J2*R*pow(rE, 2)*epsilon*mu*cos(f)) / (2.*e*pow(L, 3)*pow(eta, 3)) +
		(3.*J2*R*pow(si, 2)*pow(rE, 2)*epsilon*mu*cos(f)) / (4.*e*pow(L, 3)*pow(eta, 3)) -
		(3.*pow(G, 2)*J2*R*pow(si, 2)*pow(rE, 2)*epsilon*cos(f + 2. * g)) / (8.*e*pow(L, 3)*r*pow(eta, 3)) -
		(3.*J2*R*pow(si, 2)*pow(rE, 2)*epsilon*mu*cos(f + 2. * g)) / (8.*e*pow(L, 3)*pow(eta, 3)) -
		(3.*pow(G, 2)*J2*R*pow(si, 2)*pow(rE, 2)*epsilon*cos(2. * f + 2. * g)) / (4.*pow(e, 2)*pow(L, 3)*r*pow(eta, 3)) -
		(3.*J2*R*pow(si, 2)*pow(rE, 2)*epsilon*mu*cos(2. * f + 2. * g)) / (4.*pow(e, 2)*pow(L, 3)*pow(eta, 3)) -
		(3.*pow(G, 2)*J2*R*pow(si, 2)*pow(rE, 2)*epsilon*cos(3.*f + 2. * g)) / (8.*e*pow(L, 3)*r*pow(eta, 3)) -
		(3.*J2*R*pow(si, 2)*pow(rE, 2)*epsilon*mu*cos(3.*f + 2. * g)) / (8.*e*pow(L, 3)*pow(eta, 3)) -
		(3.*e*J2*pow(rE, 2)*epsilon*pow(mu, 2)*sin(f)) / (pow(L, 4)*pow(eta, 4)) +
		(15. * e*J2*pow(si, 2)*pow(rE, 2)*epsilon*pow(mu, 2)*sin(f)) / (4.*pow(L, 4)*pow(eta, 4)) -
		(J2*pow(rE, 2)*epsilon*pow(mu, 2)*sin(f)) / (2.*e*pow(L, 4)*pow(eta, 2)) +
		(3.*J2*pow(si, 2)*pow(rE, 2)*epsilon*pow(mu, 2)*sin(f)) / (4.*e*pow(L, 4)*pow(eta, 2)) +
		(3.*e*J2*pow(rE, 2)*epsilon*pow(mu, 2)*sin(f + 2. * g)) / (4.*pow(L, 4)*pow(eta, 4)) -
		(15. * e*J2*pow(si, 2)*pow(rE, 2)*epsilon*pow(mu, 2)*sin(f + 2. * g)) / (8.*pow(L, 4)*pow(eta, 4)) -
		(3.*J2*pow(si, 2)*pow(rE, 2)*epsilon*pow(mu, 2)*sin(f + 2. * g)) / (8.*e*pow(L, 4)*pow(eta, 2)) +
		(3.*J2*pow(rE, 2)*epsilon*pow(mu, 2)*sin(2. * f + 2. * g)) / (4.*pow(L, 4)*pow(eta, 4)) -
		(15. * J2*pow(si, 2)*pow(rE, 2)*epsilon*pow(mu, 2)*sin(2. * f + 2. * g)) / (8.*pow(L, 4)*pow(eta, 4)) +
		(e*J2*pow(rE, 2)*epsilon*pow(mu, 2)*sin(3.*f + 2. * g)) / (4.*pow(L, 4)*pow(eta, 4)) -
		(5. * e*J2*pow(si, 2)*pow(rE, 2)*epsilon*pow(mu, 2)*sin(3.*f + 2. * g)) / (8.*pow(L, 4)*pow(eta, 4)) -
		(J2*pow(si, 2)*pow(rE, 2)*epsilon*pow(mu, 2)*sin(3.*f + 2. * g)) / (8.*e*pow(L, 4)*pow(eta, 2)));

	T hmean = 
		h + (cont)*( (3.*J2*pow(rE, 2)*epsilon*pow(mu, 2)*phi) / (2.*H*pow(L, 3)*pow(eta, 3)) -
		(3.*J2*pow(si, 2)*pow(rE, 2)*epsilon*pow(mu, 2)*phi) / (2.*H*pow(L, 3)*pow(eta, 3)) +
		(3.*e*J2*pow(rE, 2)*epsilon*pow(mu, 2)*sin(f)) / (2.*H*pow(L, 3)*pow(eta, 3)) -
		(3.*e*J2*pow(si, 2)*pow(rE, 2)*epsilon*pow(mu, 2)*sin(f)) / (2.*H*pow(L, 3)*pow(eta, 3)) -
		(3.*e*J2*pow(rE, 2)*epsilon*pow(mu, 2)*sin(f + 2. * g)) / (4.*H*pow(L, 3)*pow(eta, 3)) +
		(3.*e*J2*pow(si, 2)*pow(rE, 2)*epsilon*pow(mu, 2)*sin(f + 2. * g)) / (4.*H*pow(L, 3)*pow(eta, 3)) -
		(3.*J2*pow(rE, 2)*epsilon*pow(mu, 2)*sin(2. * f + 2. * g)) / (4.*H*pow(L, 3)*pow(eta, 3)) +
		(3.*J2*pow(si, 2)*pow(rE, 2)*epsilon*pow(mu, 2)*sin(2. * f + 2. * g)) / (4.*H*pow(L, 3)*pow(eta, 3)) -
		(e*J2*pow(rE, 2)*epsilon*pow(mu, 2)*sin(3.*f + 2. * g)) / (4.*H*pow(L, 3)*pow(eta, 3)) +
		(e*J2*pow(si, 2)*pow(rE, 2)*epsilon*pow(mu, 2)*sin(3.*f + 2. * g)) / (4.*H*pow(L, 3)*pow(eta, 3)));

	T Lmean = 
		L + (cont)*( (J2*pow(rE, 2)*epsilon*pow(mu, 2)) / (2.*pow(L, 3)*pow(eta, 3)) -
		(3.*J2*pow(si, 2)*pow(rE, 2)*epsilon*pow(mu, 2)) / (4.*pow(L, 3)*pow(eta, 3)) -
		(pow(a, 2)*J2*pow(rE, 2)*epsilon*pow(mu, 2)) / (2.*pow(L, 3)*pow(r, 2)*pow(eta, 2)) +
		(3.*pow(a, 2)*J2*pow(si, 2)*pow(rE, 2)*epsilon*pow(mu, 2)) / (4.*pow(L, 3)*pow(r, 2)*pow(eta, 2)) -
		(pow(a, 2)*e*J2*pow(rE, 2)*epsilon*pow(mu, 2)*cos(f)) / (2.*pow(L, 3)*pow(r, 2)*pow(eta, 2)) +
		(3.*pow(a, 2)*e*J2*pow(si, 2)*pow(rE, 2)*epsilon*pow(mu, 2)*cos(f)) / (4.*pow(L, 3)*pow(r, 2)*pow(eta, 2)) -
		(3.*pow(a, 2)*e*J2*pow(si, 2)*pow(rE, 2)*epsilon*pow(mu, 2)*cos(f + 2. * g)) / (8.*pow(L, 3)*pow(r, 2)*pow(eta, 2)) -
		(3.*pow(a, 2)*J2*pow(si, 2)*pow(rE, 2)*epsilon*pow(mu, 2)*cos(2. * f + 2. * g)) / (4.*pow(L, 3)*pow(r, 2)*pow(eta, 2)) -
		(3.*pow(a, 2)*e*J2*pow(si, 2)*pow(rE, 2)*epsilon*pow(mu, 2)*cos(3.*f + 2. * g)) / (8.*pow(L, 3)*pow(r, 2)*pow(eta, 2)));

	T Gmean = 
		G +(cont)*(- (3.*e*J2*pow(si, 2)*pow(rE, 2)*epsilon*pow(mu, 2)*cos(f + 2. * g)) / (4.*pow(L, 3)*pow(eta, 3)) -
		(3.*J2*pow(si, 2)*pow(rE, 2)*epsilon*pow(mu, 2)*cos(2. * f + 2. * g)) / (4.*pow(L, 3)*pow(eta, 3)) -
		(e*J2*pow(si, 2)*pow(rE, 2)*epsilon*pow(mu, 2)*cos(3.*f + 2. * g)) / (4.*pow(L, 3)*pow(eta, 3)));

	T Hmean = H;

	// Mean Delaunay elements
	AlgebraicVector<T> delaunayMean(6);
	delaunayMean[0] = lmean;
	delaunayMean[1] = gmean;
	delaunayMean[2] = hmean;
	delaunayMean[3] = Lmean;
	delaunayMean[4] = Gmean;
	delaunayMean[5] = Hmean;

	return delaunayMean;
}

// Averaged J2
template<typename T> AlgebraicVector<T> averagedJ2rhs(AlgebraicVector<T> x, const double mu, const double J2, const double rE, T cont)
{
	// Delaunay elements
	T l = x[0]; // l
	T g = x[1]; // g
	T h = x[2]; // h
	T L = x[3]; // L
	T G = x[4]; // G
	T H = x[5]; // H

	T eta = G / L; // sqrt(1.0 - pow(e, 2));
	T ci = H / G; // cos(i);
	T si = sin(acos(ci));

	T dldt = pow(mu, 2) / pow(L, 3) 
	       + (cont)*((3.0 * J2*pow(rE, 2)*pow(mu, 4)) / (2.0*pow(L, 7)*pow(eta, 3)) -
		  (9.0 * J2*pow(si, 2)*pow(rE, 2)*pow(mu, 4)) / (4.0*pow(L, 7)*pow(eta, 3)));

	T dgdt = (cont)*((3. * J2*pow(rE, 2)*pow(mu, 4)) / (2.*pow(L, 7)*pow(eta, 4)) -
		(9. * J2*pow(si, 2)*pow(rE, 2)*pow(mu, 4)) / (4.*pow(L, 7)*pow(eta, 4)) +
		(3. * pow(ci, 2)*J2*pow(rE, 2)*pow(mu, 4)) / (2.*G*pow(L, 6)*pow(eta, 3)));

	T dhdt = (cont)*(-(3. * pow(ci, 2)*J2*pow(rE, 2)*pow(mu, 4)) / (2.*H*pow(L, 6)*pow(eta, 3)));

	AlgebraicVector<T> ff(6);
	ff[0] = dldt;
	ff[1] = dgdt;
	ff[2] = dhdt;
	ff[3] = 0;
	ff[4] = 0;
	ff[5] = 0;

	return ff;

}

template<typename T> AlgebraicVector<T> mean2osculating(AlgebraicVector<T> delaunayMean, double mu, double J2, double rE, T cont)
{
	// Mean Delaunay elements
	T l = delaunayMean[0]; // l
	T g = delaunayMean[1]; // g
	T h = delaunayMean[2]; // h
	T L = delaunayMean[3]; // L
	T G = delaunayMean[4]; // G
	T H = delaunayMean[5]; // H

	// Osculating Keplerian elements
	AlgebraicVector<T> kep = delaunay2kep(delaunayMean, mu);
	T a = kep[0];
	T e = kep[1];
	T incl = kep[2];

	T f = mean2trueAnomaly(l, e);	// true anomaly

	T eta = G / L;	// sqrt(1-e^2)
	T ci = cos(incl);	// cos(i)
	T si = sin(incl);	// sin(i)

	T lmod = modulusLocal(l, 2.0 * M_PI);
	T phi = modulusLocal(f - lmod, 2.0 * M_PI);	// true minus mean anomaly

	T p = a*(1. - pow(e, 2));

	T epsilon = 1.0;
	T r = p / (1. + e*cos(f));	// Radial distance

	T n = sqrt(mu / pow(a, 3));	// mean motion

	T R = e*G*sin(f) / (a*pow(eta, 2));

	T losc =
		l+(cont)*(- (J2*R*pow(rE, 2)*epsilon) / (2.*pow(e, 2)*L*r) + (3.*J2*R*pow(si, 2)*pow(rE, 2)*epsilon) / (4.*pow(e, 2)*L*r) -
		(J2*R*pow(rE, 2)*epsilon*mu) / (2.*pow(e, 2)*pow(L, 3)*pow(eta, 2)) +
		(3.*J2*R*pow(si, 2)*pow(rE, 2)*epsilon*mu) / (4.*pow(e, 2)*pow(L, 3)*pow(eta, 2)) - (J2*R*pow(rE, 2)*epsilon*cos(f)) / (2.*e*L*r) +
		(3.*J2*R*pow(si, 2)*pow(rE, 2)*epsilon*cos(f)) / (4.*e*L*r) - (J2*R*pow(rE, 2)*epsilon*mu*cos(f)) / (2.*e*pow(L, 3)*pow(eta, 2)) +
		(3.*J2*R*pow(si, 2)*pow(rE, 2)*epsilon*mu*cos(f)) / (4.*e*pow(L, 3)*pow(eta, 2)) -
		(3.*J2*R*pow(si, 2)*pow(rE, 2)*epsilon*cos(f + 2. * g)) / (8.*e*L*r) -
		(3.*J2*R*pow(si, 2)*pow(rE, 2)*epsilon*mu*cos(f + 2. * g)) / (8.*e*pow(L, 3)*pow(eta, 2)) -
		(3.*J2*R*pow(si, 2)*pow(rE, 2)*epsilon*cos(2. * f + 2. * g)) / (4.*pow(e, 2)*L*r) -
		(3.*J2*R*pow(si, 2)*pow(rE, 2)*epsilon*mu*cos(2. * f + 2. * g)) / (4.*pow(e, 2)*pow(L, 3)*pow(eta, 2)) -
		(3.*J2*R*pow(si, 2)*pow(rE, 2)*epsilon*cos(3.*f + 2. * g)) / (8.*e*L*r) -
		(3.*J2*R*pow(si, 2)*pow(rE, 2)*epsilon*mu*cos(3.*f + 2. * g)) / (8.*e*pow(L, 3)*pow(eta, 2)) -
		(J2*pow(rE, 2)*epsilon*pow(mu, 2)*sin(f)) / (2.*e*pow(L, 4)*eta) +
		(3.*J2*pow(si, 2)*pow(rE, 2)*epsilon*pow(mu, 2)*sin(f)) / (4.*e*pow(L, 4)*eta) -
		(3.*J2*pow(si, 2)*pow(rE, 2)*epsilon*pow(mu, 2)*sin(f + 2. * g)) / (8.*e*pow(L, 4)*eta) -
		(J2*pow(si, 2)*pow(rE, 2)*epsilon*pow(mu, 2)*sin(3.*f + 2. * g)) / (8.*e*pow(L, 4)*eta));

	T gosc =
		g +(cont)*((pow(G, 2)*J2*R*pow(rE, 2)*epsilon) / (2.*pow(e, 2)*pow(L, 3)*r*pow(eta, 3)) -
		(3.*pow(G, 2)*J2*R*pow(si, 2)*pow(rE, 2)*epsilon) / (4.*pow(e, 2)*pow(L, 3)*r*pow(eta, 3)) +
		(J2*R*pow(rE, 2)*epsilon*mu) / (2.*pow(e, 2)*pow(L, 3)*pow(eta, 3)) -
		(3.*J2*R*pow(si, 2)*pow(rE, 2)*epsilon*mu) / (4.*pow(e, 2)*pow(L, 3)*pow(eta, 3)) +
		(3.*J2*pow(rE, 2)*epsilon*pow(mu, 2)*phi) / (pow(L, 4)*pow(eta, 4)) -
		(15. * J2*pow(si, 2)*pow(rE, 2)*epsilon*pow(mu, 2)*phi) / (4.*pow(L, 4)*pow(eta, 4)) +
		(pow(G, 2)*J2*R*pow(rE, 2)*epsilon*cos(f)) / (2.*e*pow(L, 3)*r*pow(eta, 3)) -
		(3.*pow(G, 2)*J2*R*pow(si, 2)*pow(rE, 2)*epsilon*cos(f)) / (4.*e*pow(L, 3)*r*pow(eta, 3)) +
		(J2*R*pow(rE, 2)*epsilon*mu*cos(f)) / (2.*e*pow(L, 3)*pow(eta, 3)) -
		(3.*J2*R*pow(si, 2)*pow(rE, 2)*epsilon*mu*cos(f)) / (4.*e*pow(L, 3)*pow(eta, 3)) +
		(3.*pow(G, 2)*J2*R*pow(si, 2)*pow(rE, 2)*epsilon*cos(f + 2. * g)) / (8.*e*pow(L, 3)*r*pow(eta, 3)) +
		(3.*J2*R*pow(si, 2)*pow(rE, 2)*epsilon*mu*cos(f + 2. * g)) / (8.*e*pow(L, 3)*pow(eta, 3)) +
		(3.*pow(G, 2)*J2*R*pow(si, 2)*pow(rE, 2)*epsilon*cos(2. * f + 2. * g)) / (4.*pow(e, 2)*pow(L, 3)*r*pow(eta, 3)) +
		(3.*J2*R*pow(si, 2)*pow(rE, 2)*epsilon*mu*cos(2. * f + 2. * g)) / (4.*pow(e, 2)*pow(L, 3)*pow(eta, 3)) +
		(3.*pow(G, 2)*J2*R*pow(si, 2)*pow(rE, 2)*epsilon*cos(3.*f + 2. * g)) / (8.*e*pow(L, 3)*r*pow(eta, 3)) +
		(3.*J2*R*pow(si, 2)*pow(rE, 2)*epsilon*mu*cos(3.*f + 2. * g)) / (8.*e*pow(L, 3)*pow(eta, 3)) +
		(3.*e*J2*pow(rE, 2)*epsilon*pow(mu, 2)*sin(f)) / (pow(L, 4)*pow(eta, 4)) -
		(15. * e*J2*pow(si, 2)*pow(rE, 2)*epsilon*pow(mu, 2)*sin(f)) / (4.*pow(L, 4)*pow(eta, 4)) +
		(J2*pow(rE, 2)*epsilon*pow(mu, 2)*sin(f)) / (2.*e*pow(L, 4)*pow(eta, 2)) -
		(3.*J2*pow(si, 2)*pow(rE, 2)*epsilon*pow(mu, 2)*sin(f)) / (4.*e*pow(L, 4)*pow(eta, 2)) -
		(3.*e*J2*pow(rE, 2)*epsilon*pow(mu, 2)*sin(f + 2. * g)) / (4.*pow(L, 4)*pow(eta, 4)) +
		(15. * e*J2*pow(si, 2)*pow(rE, 2)*epsilon*pow(mu, 2)*sin(f + 2. * g)) / (8.*pow(L, 4)*pow(eta, 4)) +
		(3.*J2*pow(si, 2)*pow(rE, 2)*epsilon*pow(mu, 2)*sin(f + 2. * g)) / (8.*e*pow(L, 4)*pow(eta, 2)) -
		(3.*J2*pow(rE, 2)*epsilon*pow(mu, 2)*sin(2. * f + 2. * g)) / (4.*pow(L, 4)*pow(eta, 4)) +
		(15. * J2*pow(si, 2)*pow(rE, 2)*epsilon*pow(mu, 2)*sin(2. * f + 2. * g)) / (8.*pow(L, 4)*pow(eta, 4)) -
		(e*J2*pow(rE, 2)*epsilon*pow(mu, 2)*sin(3.*f + 2. * g)) / (4.*pow(L, 4)*pow(eta, 4)) +
		(5. * e*J2*pow(si, 2)*pow(rE, 2)*epsilon*pow(mu, 2)*sin(3.*f + 2. * g)) / (8.*pow(L, 4)*pow(eta, 4)) +
		(J2*pow(si, 2)*pow(rE, 2)*epsilon*pow(mu, 2)*sin(3.*f + 2. * g)) / (8.*e*pow(L, 4)*pow(eta, 2)));

	T hosc =
		h + (cont)*( -(3.*J2*pow(rE, 2)*epsilon*pow(mu, 2)*phi) / (2.*H*pow(L, 3)*pow(eta, 3)) +
		(3.*J2*pow(si, 2)*pow(rE, 2)*epsilon*pow(mu, 2)*phi) / (2.*H*pow(L, 3)*pow(eta, 3)) -
		(3.*e*J2*pow(rE, 2)*epsilon*pow(mu, 2)*sin(f)) / (2.*H*pow(L, 3)*pow(eta, 3)) +
		(3.*e*J2*pow(si, 2)*pow(rE, 2)*epsilon*pow(mu, 2)*sin(f)) / (2.*H*pow(L, 3)*pow(eta, 3)) +
		(3.*e*J2*pow(rE, 2)*epsilon*pow(mu, 2)*sin(f + 2. * g)) / (4.*H*pow(L, 3)*pow(eta, 3)) -
		(3.*e*J2*pow(si, 2)*pow(rE, 2)*epsilon*pow(mu, 2)*sin(f + 2. * g)) / (4.*H*pow(L, 3)*pow(eta, 3)) +
		(3.*J2*pow(rE, 2)*epsilon*pow(mu, 2)*sin(2. * f + 2. * g)) / (4.*H*pow(L, 3)*pow(eta, 3)) -
		(3.*J2*pow(si, 2)*pow(rE, 2)*epsilon*pow(mu, 2)*sin(2. * f + 2. * g)) / (4.*H*pow(L, 3)*pow(eta, 3)) +
		(e*J2*pow(rE, 2)*epsilon*pow(mu, 2)*sin(3.*f + 2. * g)) / (4.*H*pow(L, 3)*pow(eta, 3)) -
		(e*J2*pow(si, 2)*pow(rE, 2)*epsilon*pow(mu, 2)*sin(3.*f + 2. * g)) / (4.*H*pow(L, 3)*pow(eta, 3)));

	T Losc =
		L + (cont)*(- (J2*pow(rE, 2)*epsilon*pow(mu, 2)) / (2.*pow(L, 3)*pow(eta, 3)) +
		(3.*J2*pow(si, 2)*pow(rE, 2)*epsilon*pow(mu, 2)) / (4.*pow(L, 3)*pow(eta, 3)) +
		(pow(a, 2)*J2*pow(rE, 2)*epsilon*pow(mu, 2)) / (2.*pow(L, 3)*pow(r, 2)*pow(eta, 2)) -
		(3.*pow(a, 2)*J2*pow(si, 2)*pow(rE, 2)*epsilon*pow(mu, 2)) / (4.*pow(L, 3)*pow(r, 2)*pow(eta, 2)) +
		(pow(a, 2)*e*J2*pow(rE, 2)*epsilon*pow(mu, 2)*cos(f)) / (2.*pow(L, 3)*pow(r, 2)*pow(eta, 2)) -
		(3.*pow(a, 2)*e*J2*pow(si, 2)*pow(rE, 2)*epsilon*pow(mu, 2)*cos(f)) / (4.*pow(L, 3)*pow(r, 2)*pow(eta, 2)) +
		(3.*pow(a, 2)*e*J2*pow(si, 2)*pow(rE, 2)*epsilon*pow(mu, 2)*cos(f + 2. * g)) / (8.*pow(L, 3)*pow(r, 2)*pow(eta, 2)) +
		(3.*pow(a, 2)*J2*pow(si, 2)*pow(rE, 2)*epsilon*pow(mu, 2)*cos(2. * f + 2. * g)) / (4.*pow(L, 3)*pow(r, 2)*pow(eta, 2)) +
		(3.*pow(a, 2)*e*J2*pow(si, 2)*pow(rE, 2)*epsilon*pow(mu, 2)*cos(3.*f + 2. * g)) / (8.*pow(L, 3)*pow(r, 2)*pow(eta, 2)));

	T Gosc =
		G + (cont)*( (3.*e*J2*pow(si, 2)*pow(rE, 2)*epsilon*pow(mu, 2)*cos(f + 2. * g)) / (4.*pow(L, 3)*pow(eta, 3)) +
		(3.*J2*pow(si, 2)*pow(rE, 2)*epsilon*pow(mu, 2)*cos(2. * f + 2. * g)) / (4.*pow(L, 3)*pow(eta, 3)) +
		(e*J2*pow(si, 2)*pow(rE, 2)*epsilon*pow(mu, 2)*cos(3.*f + 2. * g)) / (4.*pow(L, 3)*pow(eta, 3)));

	T Hosc = H;

	// Osculating Delaunay elements
	AlgebraicVector<T> delaunayOsc(6);
	delaunayOsc[0] = losc;
	delaunayOsc[1] = gosc;
	delaunayOsc[2] = hosc;
	delaunayOsc[3] = Losc;
	delaunayOsc[4] = Gosc;
	delaunayOsc[5] = Hosc;

	return delaunayOsc;
}

template <typename T> AlgebraicVector<T> kep2cart(const AlgebraicVector<T>& kep, const double mu = 398600.4415)
{
	/*member function to convert keplerian  classical element into Earth-Centred inertial reference frame element
	!< keplerian element kep = {a, e, i, RA, PA, TA}
	RA: rigth ascension of ascending node; PA: argument of periapsis; TA: true anomaly; i:orbital inclination
	!> return AlgebraicVector of ECI reference frame res = {x, y, z, dx, dy, dz}*/

	//const double mu = 398600.4415;

	T p = kep[0] * (1.0 - kep[1] * kep[1]);

	// position and velocity in perifocal refererence frame
	AlgebraicVector<T> rm(3), vm(3);
	rm[0] = p*cos(kep[5]) / (1.0 + kep[1] * cos(kep[5]));
	rm[1] = p*sin(kep[5]) / (1.0 + kep[1] * cos(kep[5]));
	rm[2] = 0.0;
	vm[0] = -1.0*sin(kep[5])*sqrt(mu / p);
	vm[1] = (kep[1] + cos(kep[5]))*sqrt(mu / p);
	vm[2] = 0.0;

	T cRA = cos(kep[3]);
	T sRA = sin(kep[3]);
	T cPA = cos(kep[4]);
	T sPA = sin(kep[4]);
	T ci = cos(kep[2]);
	T si = sin(kep[2]);

	T RR[3][3]; // rotational matrix from perifocal to eci reference frame
	RR[0][0] = cRA*cPA - sRA*ci*sPA;  RR[0][1] = -1.0*cRA*sPA - sRA*ci*cPA; RR[0][2] = sRA*si;
	RR[1][0] = sRA*cPA + cRA*ci*sPA;  RR[1][1] = -1.0*sRA*sPA + cRA*ci*cPA; RR[1][2] = -1.0*cRA*si;
	RR[2][0] = si*sPA;                RR[2][1] = si*cPA;					RR[2][2] = ci;

	AlgebraicVector<T> rr(3), vv(3);
	for (unsigned int i = 0; i<3; i++){
		rr[i] = 0.0;
		vv[i] = 0.0;
		for (unsigned int j = 0; j<3; j++){
			rr[i] = rr[i] + RR[i][j] * rm[j];
			vv[i] = vv[i] + RR[i][j] * vm[j];
		}
	}

	AlgebraicVector<T> res(6);
	res[0] = rr[0];
	res[1] = rr[1];
	res[2] = rr[2];
	res[3] = vv[0];
	res[4] = vv[1];
	res[5] = vv[2];

	return res;
}


template<typename T> AlgebraicVector<T> analyticJ2propHill(AlgebraicVector<T> xx0, double tof, T cont)
{
    
    const double mu = 398600.4418; //{[km^3/s^2]}
    const double rE = 6378.137; //{km}
    const double J2 = 1.08262668e-3; //{-}
    
    AlgebraicVector<T> kep0(6), kep0Mean(6), kepfMean(6), del0Mean(6), delfMean, hill0(6), hill0Mean(6), hillfMean(6), hillf(6), kepf(6), xxf(6);
    
    // transform xx0 to keplerian elements
    double Lsc = rE;
    double Vsc = sqrt(mu/rE);
    double Tsc = Lsc/Vsc;
    double muSc = mu/Lsc/Lsc/Lsc*Tsc*Tsc;
    
    xx0[0] = xx0[0]/Lsc;
    xx0[1] = xx0[1]/Lsc;
    xx0[2] = xx0[2]/Lsc;
    xx0[3] = xx0[3]/Vsc;
    xx0[4] = xx0[4]/Vsc;
    xx0[5] = xx0[5]/Vsc;
    
    kep0 = cart2kep(xx0, muSc); //-> convert true to mean anomaly!
    
    // trasnform keplerian elements to Hill
    
    hill0 = kep2hill(kep0, muSc);
    
    // from osculating to mean
    hill0Mean = osculating2meanHill(hill0, muSc, J2, rE/Lsc, cont);
    
    kep0Mean = hill2kep(hill0Mean, muSc);
    T meanAnomaly = true2meanAnomaly(kep0Mean[5], kep0Mean[1]);
    kep0Mean[5] = meanAnomaly;
    
    del0Mean = kep2delaunay(kep0Mean, muSc);
    
    delfMean = averagedJ2rhs(del0Mean, muSc, J2, rE/Lsc, cont);
    delfMean = delfMean*tof/Tsc+del0Mean;
    
    kepfMean = delaunay2kep(delfMean, muSc);
    
    T trueAnomaly = mean2trueAnomaly(kepfMean[5], kepfMean[1]);
    kepfMean[5] = trueAnomaly;
    
    hillfMean = kep2hill(kepfMean, muSc);
    
    // transform mean to osculating
    
    hillf = mean2osculatingHill(hillfMean, muSc, J2, rE/Lsc, cont);
    
    // transform keplerian elements to xxf
    
    xxf = hill2cart(hillf, muSc);
    
    xxf[0] = xxf[0]*Lsc;
    xxf[1] = xxf[1]*Lsc;
    xxf[2] = xxf[2]*Lsc;
    xxf[3] = xxf[3]*Vsc;
    xxf[4] = xxf[4]*Vsc;
    xxf[5] = xxf[5]*Vsc;
    
    // cout << "xxf" << xxf << endl;
    
    return xxf;
    
}


template<typename T> AlgebraicVector<T> analyticJ2prop(AlgebraicVector<T> xx0, double tof, T cont)
{

const double mu = 398600.4418; //{[km^3/s^2]}
const double rE = 6378.137; //{km}
const double J2 = 1.08262668e-3; //{-}

AlgebraicVector<T> kep0(6), del0(6), del0Mean(6), delfMean(6), delf(6), kepf(6), xxf(6);

// transform xx0 to keplerian elements 
double Lsc = rE;
double Vsc = sqrt(mu/rE);
double Tsc = Lsc/Vsc;
double muSc = mu/Lsc/Lsc/Lsc*Tsc*Tsc;

xx0[0] = xx0[0]/Lsc;
xx0[1] = xx0[1]/Lsc;
xx0[2] = xx0[2]/Lsc;
xx0[3] = xx0[3]/Vsc;
xx0[4] = xx0[4]/Vsc;
xx0[5] = xx0[5]/Vsc;


kep0 = cart2kep(xx0, muSc); //-> convert true to mean anomaly!

T meanAnomaly = true2meanAnomaly(kep0[5], kep0[1]);
kep0[5] = meanAnomaly; 
 
// trasnform keplerian elements to Delanauy 

del0 = kep2delaunay(kep0, muSc);

// from osculating to mean 

del0Mean = osculating2mean(del0, muSc, J2, rE/Lsc, cont);

// propagation of mean elements for tof

delfMean = averagedJ2rhs(del0Mean, muSc, J2, rE/Lsc, cont);

delfMean = delfMean*tof/Tsc+del0Mean;


// transform mean to osculating 

delf = mean2osculating(delfMean, muSc, J2, rE/Lsc, cont);

// transform osculating Delanay to keplerian elements 

kepf = delaunay2kep(delf, muSc);

T trueAnomaly = mean2trueAnomaly(kepf[5], kepf[1]);
kepf[5] = trueAnomaly; 

// transform keplerian elements to xxf

xxf = kep2cart(kepf, muSc);

xxf[0] = xxf[0]*Lsc;
xxf[1] = xxf[1]*Lsc;
xxf[2] = xxf[2]*Lsc;
xxf[3] = xxf[3]*Vsc;
xxf[4] = xxf[4]*Vsc;
xxf[5] = xxf[5]*Vsc;


// cout << "xxf" << xxf << endl;

return xxf;

}


// int main( void )
// {
    
//     int order, nvar;
//     double maxCont;
//     int i;

    
//     ifstream infile;
//     infile.open("input.dat");
//     infile >> order;
//     infile >> maxCont;
//     infile.close();
    
//     nvar = 4;
//     DA::init( order, nvar ); //order and nvar
//     DA::setEps(1e-300);
    
//     AlgebraicVector<DA> mapD(4), mapI(4), dx(4), x0plus(6), dxDA(4);
//     AlgebraicVector<DA> x0DA(6), xfDA(6);

//     AlgebraicVector<double> x0(6), xf(6), xfc(6), consMapD(4), cr(3);
      
// 	double tof;
	
// 	infile.open("initialStateKep.dat");
// 	if (infile.is_open())
// 	{
// 		infile >> x0[0];
// 		infile >> x0[1];
// 		infile >> x0[2];
// 		infile >> x0[3];
// 		infile >> x0[4];
// 		infile >> x0[5];
// 		infile >> tof;
// 	}
// 	else
// 	{
// 		cout << "Input file not found!" << endl;
// 		return 1;
// 	}
// 	infile.close();
		
// 	infile.open("desiredStateKep.dat");
// 	if (infile.is_open())
// 	{
// 		infile >> xf[0];
// 		infile >> xf[1];
// 		infile >> xf[2];
// 		infile >> xf[3];
// 		infile >> xf[4];
// 		infile >> xf[5];
// 	}
// 	else
// 	{
// 		cout << "Desired final state file not found!" << endl;
// 		return 1;
// 	}
// 	infile.close();
	
//     xfc = analyticJ2propHill(x0, tof, 1.0);
    
//     ofstream outresiduals;
//     outresiduals.open("residuals.dat");
//     outresiduals << setprecision(16);
//     for (int i = 0; i<3 ; i++) {
//         outresiduals << xfc[i]-xf[i]<<endl;};

//     ofstream velocity;
//     velocity.open("velocity.dat");
//     velocity << setprecision(16);
 
//     double dJ2 = 1.0/maxCont;
//     double cont = 0.0;
//     double scl = 0.01;

    
//     for (int j = 1; j <= maxCont; j++){	
        		    	
// 	for (int i = 0; i < 3; i++)
// 	{
// 		x0DA[i] = x0[i];
// 	}
	
	
// 	for (int i = 3; i < 6; i++)
// 	{
// 		x0DA[i] = x0[i] + scl*DA(i-2);
// 	}
	
	
// 	xfDA = analyticJ2propHill(x0DA, tof, (cont+DA(4)));
     
 
//     for (int i = 0; i<3 ; i++) {
//             velocity << cons(x0DA[i+3])<<endl;};
        
//     for (int i = 0; i<3 ; i++) {
//             velocity << cons(xfDA[i+3])<<endl;};


// 	mapD[0] = xfDA[0]-xf[0];
// 	mapD[1] = xfDA[1]-xf[1];
// 	mapD[2] = xfDA[2]-xf[2];
// 	mapD[3] = DA(4);
	
//     consMapD = mapD.cons();
        
//     //    cout << mapD << endl;
     
//     if (j>1) {
//         for (i = 0; i<3 ; i++) {
//         outresiduals << consMapD[i] <<endl;};
//         }
        
//     mapD = mapD-consMapD;
//     mapI = mapD.invert();
    
                
//     dx[0] = -consMapD[0] +0*DA(1);
//     dx[1] = -consMapD[1] +0*DA(2);
//     dx[2] = -consMapD[2] +0*DA(3);
//     dx[3] = dJ2;
    
//     dx = mapI.eval(dx);
    
//     for (i = 0; i<6 ; i++) {
//         x0plus[i] = x0[i];
//     }
     
//     for (i = 3; i<6 ; i++) {
//         x0plus[i] = x0plus[i]+scl*dx[i-3];
//     }
    
//    for (i = 0; i<6 ; i++) {
//         x0[i] = x0plus[i].cons();
//     }
    
//     cont = cont+dJ2;
        
//     cout << "current continuation" << endl << cont << endl;
        
//     // eliminate all the errors at that iteration
        
//     for (int i = 0; i < 3; i++)
//     {
//         x0DA[i] = x0[i];
//     }
        
//     for (int i = 3; i < 6; i++)
//         {
//         x0DA[i] = x0[i] + scl*DA(i-2);
//     }
        
//     xfDA = analyticJ2propHill(x0DA, tof, (cont +0*DA(4)));
        
//     mapD[0] = xfDA[0]-xf[0];
//     mapD[1] = xfDA[1]-xf[1];
//     mapD[2] = xfDA[2]-xf[2];
//     mapD[3] = DA(4);
        
//     consMapD = mapD.cons();
        
//     mapD = mapD-consMapD;
//     mapI = mapD.invert();
        
//     dx[0] = -consMapD[0] +0*DA(1);
//     dx[1] = -consMapD[1] +0*DA(2);
//     dx[2] = -consMapD[2] +0*DA(3);
//     dx[3] = 0*DA(4);
        
//     dx = mapI.eval(dx);
        
//     for (i = 0; i<6 ; i++) {
//           x0plus[i] = x0[i];
//     }
        
//     for (i = 3; i<6 ; i++) {
//         x0plus[i] = x0plus[i]+scl*dx[i-3];
//     }
        
//     for (i = 0; i<6 ; i++) {
//         x0[i] = x0plus[i].cons();
//     }
    
    
//     }
    
// 	// validation of the results 
// 	AlgebraicVector<double> x0plusCons(6), xfCons(6), dxCons(4);	
//     x0plusCons = x0plus.cons();
    
// 	xfCons = analyticJ2propHill(x0plusCons, tof, 1.0);

//     for (i = 0; i<3 ; i++) {
//         outresiduals << xfCons[i] - xf[i] <<endl;};
    
//     for (int i = 0; i<3 ; i++) {
//         velocity << x0plusCons[i+3]<<endl;};
    
//     for (int i = 0; i<3 ; i++) {
//         velocity << xfCons[i+3]<<endl;};
    
//     outresiduals.close();
//     velocity.close();

// 	// write post maneuver conditions
// 	ofstream outfile;
// 	outfile.open("postManeuverStateFromKep.dat");
// 	outfile << setprecision(16);
// 	for (i = 0; i<6 ; i++) {
// 	outfile << x0plusCons[i]<<endl;};
// 	outfile.close();
	
//     // final state conditions
// 	outfile.open("finalStateValidationFromKep.dat");
// 	outfile << setprecision(16);
// 	for (i = 0; i<6 ; i++) {
// 	outfile << xfCons[i]<<endl;};
// 	outfile.close();

// }

template<typename T> AlgebraicVector<T> oscCart2MeanCart(const AlgebraicVector<T>& rv, const double mu, const double J2, const double rE)
{
     	AlgebraicVector<T> xcart = cart2kep(rv, mu);

        AlgebraicVector<T> xhill  = kep2hill(xcart, mu);

        T cont = 1.0;

        AlgebraicVector<T> xhillMean  = osculating2meanHill(xhill,mu, J2, rE, cont);

        // hill to kepler 
         AlgebraicVector<T>  xcartmean = hill2cart(xhillMean, mu);

        return xcartmean;

}


#endif
