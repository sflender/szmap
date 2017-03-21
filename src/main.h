#include <stdio.h>
#include <stdlib.h>
#include <math.h>
using namespace std;

float return_P200(float rho_crit, float f_b, float M200, float R200){
    double G = 4.3018e-9; //Mpc Mo^-1 (km/s)^2 
    double P200 = 200.0 * rho_crit * f_b * G * M200  / (2.0 * R200); //Msun km^2 / Mpc^3 / s^2
    return P200; //
}

double yproj_func_Battaglia(double l, void* p){
    float *params = (float *) p;
    float xi = params[0];
    float M200 = params[1];
    float z = params[2];
    float x = sqrt(xi*xi+l*l);
    float P0 = 18.1 *  pow(M200 / 1e14, 0.154)     * pow(1. + z , -0.758);
    float xc = 0.497 * pow(M200 / 1e14,-0.00865) * pow(1. + z, 0.731);
    float beta = 4.35 *pow(M200 / 1e14,0.0393)   * pow(1. + z, 0.415);
    float gamma = -0.3;
    return P0 * pow(x/xc,gamma) * pow( 1.0+(x/xc) , (-1.0)*beta );
}

float return_ysz_Battaglia(float xmax, float xi, float M200, float R200, float z, float rho_crit, float f_b){

    // xmax in units of R200 == maximum radius to integrate to
    // xi in units of R200 == projected distance to center

	double Msun_to_kg = 1.9889e30; //kg
	double Mpc_to_cm = 3.0857e24;
	double sigma_T = 6.652e-25; //cm^2
	double m_e = 9.11e-31; //kg
	double Joules_to_keV = 6.24e15;// joules in keV
	double clight = 2.9979e8; //m/s
	double mecsq = m_e * pow(clight,2) * Joules_to_keV; // keV

    gsl_integration_workspace * w = gsl_integration_workspace_alloc (10000); 

    double result, error;
    float params[3] = {xi, M200, z};
    gsl_function F;
    F.function = &yproj_func_Battaglia;
    F.params = &params;

    gsl_integration_qags (&F, 0, sqrt(xmax*xmax-xi*xi), 0, 1e-6, 10000, w, &result, &error); 

    double P200 = return_P200(rho_crit, f_b, M200, R200); // Msun km^2 / Mpc^3 / s^2
    P200 *= Msun_to_kg * 1e6 * Joules_to_keV / pow(Mpc_to_cm,3); // now in keV/cm^3

    float this_ysz = P200 * (sigma_T/mecsq) * 0.518 * (R200*Mpc_to_cm) * 2.0 * result;

    gsl_integration_workspace_free (w);
    return this_ysz;
}