// cosmo.h 
// C++ cosmology library by Suman Bhattacharya 
// https://github.com/bsuman79/Cosmocal

#include <stdio.h>
#include <math.h>
#include <iostream>
#include <stdlib.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_integration.h>

using namespace std;

class cosmo {

    protected:
        float H0, Omega_M, Omega_b, rho_crit, a, Omega_k, Omega_L, wt;
        float PI, m_sun, G, mpc, gnewt;

    public:

    cosmo(float inp1, float inp2, float inp3, float inp4, float inp5){
        
        H0 = inp1;
        Omega_M = inp2;
        Omega_b = inp3;
        Omega_k = inp4;
        Omega_L = 1.0- Omega_M- Omega_k;
        wt= inp5;
        PI = 4.*atan(1.);
        m_sun = 1.98892e30;
        G = 6.67e-11/1.0e9; // in km^3 kg^-1 s^-2
        mpc = 3.0857e19; // in km
        gnewt = G*m_sun/mpc; //for cosm params (in km^2 Mpc msun^-1 s^-2)
    }

    float get_Omega_M(){
        return Omega_M;
    }

    float get_Omega_b(){
        return Omega_b;
    }

    float get_Omega_L(){
        return Omega_L;
    }  

    float get_Omega_k(){
        return Omega_k;
    }  

    float get_wt(){
        return wt;
    }

    float get_H0() {
        return H0;
    }


    float scale_fact(float redshift) {
        a = 1.0/(1.0+redshift);
        return a;
    }

    float hubblez(float redshift) {
        scale_fact(redshift);
        return H0 * Efact(redshift);
    }

    float calc_rho_crit(float redshift) { // this returns rho_crit in units Msun/Mpc^3
        float Hz = hubblez(redshift);
        rho_crit = pow(Hz,2)*3./(8.*PI*gnewt);
        return rho_crit;
    }

    float Omega_Mz(float redshift) {
        scale_fact(redshift);
        float Hz = hubblez(redshift);
        return (Omega_M/pow(a,3))*pow(H0/Hz,2);
    }

    float Delta_vir(float redshift) {
        float x;
        x  = Omega_Mz(redshift) - 1.0;
        return (18*pow(PI,2)) + (82*x) - (39*pow(x,2));
    }

    float Efact(float redshift) {
        scale_fact(redshift);
        return sqrt(Omega_M/pow(a,3) + Omega_L/pow(a, 3*wt+3) + (1.0-Omega_M-Omega_L)/pow(a,2));
    }

    float ang_diam(float redshift) {

        int dummy;
        float result;
        double x[1000], y[1000];
        double units = scale_fact(redshift)*3.0e5/H0;
   
        for(dummy=0;dummy<1000;dummy++){
            x[dummy]= 0.0+ redshift/999*dummy;
            y[dummy] = 1.0/(sqrt( Omega_M*pow(1+x[dummy],3) + Omega_L*pow(1+x[dummy], 3*wt+3) + (1.0-Omega_M-Omega_L)*pow(1+x[dummy],2)));
        }
        if (redshift==0.0) return 0.0;
        else {
            gsl_interp_accel *acc = gsl_interp_accel_alloc ();
            gsl_spline *spline = gsl_spline_alloc (gsl_interp_cspline, 1000);
            gsl_spline_init (spline, x, y, 1000);
            result = units*gsl_spline_eval_integ (spline, x[0], x[dummy-1], acc); 
            gsl_spline_free (spline);
            gsl_interp_accel_free (acc);
        }
        return result;
    }

    float comov_dist(float redshift){
        return (1.0+redshift)*ang_diam(redshift);
    }

};