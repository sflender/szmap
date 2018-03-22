#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_errno.h>
using namespace std;

float HuKravtsov(float z, float M200, float rhocrit, float Delta_vir, float h){
    //Eq. C10 in Hu&Kravstov to convert any mass to virial mass
    float a1 = 0.5116;
    float a2 = -0.4283;
    float a3 = -3.13e-3;
    float a4 = -3.52e-5;
    float Delta = Delta_vir / 200;
    float conc = (5.71 / pow(1 + z, 0.47)) * pow(M200 * h / 2e12, -0.084);//Duffy concetration from M200
    float B = -0.084;
    float R = pow(M200 / ((4 * M_PI / 3.) * 200 * rhocrit), 1/3.); //(Msun / Msun Mpc^(-3))1/3. -> Mpc    
    float Rs = R / conc;

    float A = log(1.+conc) - 1. + 1. / (1. + conc);
    float f = Delta * A / pow(conc, 3);
    float p = a2 + a3 * log(f) + a4 * pow(log(f), 2);
    float x = pow(a1 * pow(f, 2.*p) + pow(0.75, 2.), -0.5) + 2. * f;
    float Mvir = M200 * Delta * pow(1./(conc*x), 3);
    float Rvir = pow(3. * Mvir / Delta_vir / rhocrit / 4. / M_PI, 1./3.);
    return Mvir;
}

float return_P200(float rhocrit, float f_b, float M200, float R200){
    double G = 4.3018e-9; //Mpc Mo^-1 (km/s)^2 
    double P200 = 200.0 * rhocrit * f_b * G * M200  / (2.0 * R200); //Msun km^2 / Mpc^3 / s^2
    return P200; //
}

double yproj_func_Battaglia(double l, void* p){
    float *params = (float *) p;
    float xi = params[0];
    float M200 = params[1];
    float z = params[2];
    //printf("xi %f M200 %.2e z %f l %f\n", xi, M200, z, l);

    float x = sqrt(xi*xi+l*l);
    float P0 = 18.1 *  pow(M200 / 1e14, 0.154)     * pow(1. + z , -0.758);
    float xc = 0.497 * pow(M200 / 1e14,-0.00865) * pow(1. + z, 0.731);
    float beta = 4.35 *pow(M200 / 1e14,0.0393)   * pow(1. + z, 0.415);
    float gamma = -0.3;
    return P0 * pow(x/xc,gamma) * pow( 1.0+(x/xc) , (-1.0)*beta );
}

float return_ysz_Battaglia(float xmax, float xi, float M200, float R200, float z, float rhocrit, float f_b){

    // xmax in units of R200 == maximum radius to integrate to
    // xi in units of R200 == projected distance to center

	double Msun_to_kg = 1.9889e30; //kg
	double Mpc_to_cm = 3.0857e24;
	double sigma_T = 6.652e-25; //cm^2
	double m_e = 9.11e-31; //kg
	double Joules_to_keV = 6.24e15;// joules in keV
	double clight = 2.9979e8; //m/s
	double mecsq = m_e * pow(clight,2) * Joules_to_keV; // keV

    gsl_integration_workspace * w = gsl_integration_workspace_alloc (1e6); 

    double result, error;
    float params[3] = {xi, M200, z};
    gsl_function F;
    F.function = &yproj_func_Battaglia;
    F.params = &params;

    gsl_integration_qags (&F, 0, sqrt(xmax*xmax-xi*xi), 0, 1e-6, 1e6 , w, &result, &error); 

    double P200 = return_P200(rhocrit, f_b, M200, R200); // Msun km^2 / Mpc^3 / s^2
    P200 *= Msun_to_kg * 1e6 * Joules_to_keV / pow(Mpc_to_cm,3); // now in keV/cm^3

    float this_ysz = P200 * (sigma_T/mecsq) * 0.518 * (R200*Mpc_to_cm) * 2.0 * result;

    gsl_integration_workspace_free (w);
    return this_ysz;
}


float return_ysz_Beta_Birkinshaw(float beta, float ti, float M200, float z, float h, float Omega_M, float Omega_Mz, float Delta_vir, float ang_dia, float rhocrit){
    double sigma_T = 6.652e-25; //cm^2
    double Mpc_to_cm = 3.0857e24;
    double Mpc_to_cm_sq = 9.52154449e+48;
    float mecsq = 510.998; //keV


    float Mvir = HuKravtsov(z, M200, rhocrit, Delta_vir, h); 
    float Rvir = pow(3. * Mvir / Delta_vir / rhocrit / 4. / M_PI, 1./3.); //Mpc
    double Ne = 1.76e56 * Mvir; //Number of electrons

    float conc = (7.85 / pow(1. + z, 0.71)) * pow(Mvir * h / 2e12, -0.081); //Duffy concentration using Mvir
    float Rs = Rvir / conc; //scale radius of NFW profile
    

    gsl_sf_result gamma_numerator, gamma_denominator;

    gsl_sf_gamma_e(3. * beta - 0.5, &gamma_numerator);
    gsl_sf_gamma_e(3. * beta, &gamma_denominator);

    float Oz = Omega_M * Delta_vir / Omega_Mz;
    float Mnorm = Mvir * h / 1e15; 
    //float Rvir_WB = 9.5103 * pow(Oz, -1/3.) * pow(Mnorm, 1/3.) / h / (1 + z); //Based on Waizmann & Bartelmann Eq. 13. I found that Rvir and Rvir_WB in very good agreement

    float kbt = (1 + z)  * pow(Oz, 1/3.) * pow(Mnorm, 2/3.) / beta; //keV. Eq. 14 of WB
    float xiz = 0.14 * pow(1 + z, 1/5.); //WB Eq. 16
    float rc_1 = Rvir * xiz;
    float rc = 0.2 * Rs; //Sunil Golwala Sec. 2.2
    //printf("rc_1 %f rc %f\n", rc_1, rc);
    float thetac = rc / ang_dia;
    double ne0 = Ne / 4. / M_PI / pow(rc, 3.) / (1 / xiz + atan(xiz) - M_PI / 2.); //Number / Mpc^3. Eq. 17 of WB
    //printf("ti %f tc %f z %f Mvir %.2e Rvir %f rc %f ne0 %.2e kbt %f \n", ti, thetac, z, Mvir, Rvir, rc, ne0, kbt);
    double y0 = kbt * sigma_T * rc * ne0 * sqrt(M_PI) * gamma_numerator.val / (mecsq * Mpc_to_cm_sq * gamma_denominator.val); //Eqs. 68 & 69 of Birkinshaw 1999. rc in Mpc, ne in Mpc^3 amd I converted the unit of the equation from Mpc^2 to cm^2 
    
    double this_ysz = y0 * pow((1. + pow(ti / thetac, 2)), (0.5 - 3. * beta / 2.));// Birkinshaw Eq. 66
    //printf("y0 %.2e this_ysz %.2e\n", y0, this_ysz);
    return this_ysz;
}

double yproj_func_Beta_WB(double l, void* p){
    float *params = (float *) p;
    float ti = params[0];
    float xiz = params[1];
    float thetac = params[2];
    //printf("xiz %f, thetac %f l %f\n",xiz, thetac, l);
    float ixizsq = 1 / (xiz * xiz);
    float ttcsq = l * l / (thetac * thetac);
    float targ = sqrt((ixizsq - ttcsq) / (1 + ttcsq));
    return atan(targ) / sqrt(1 + ttcsq);
}


float return_ysz_Beta_WB(float tmax, float ti, float M200, float z, float h, float Omega_M, float Omega_Mz, float Delta_vir, float ang_dia, float rhocrit){
    //Using Eq. 14, 16, 17 and 19 of Waizmann & Bartelmann 
    double sigma_T = 6.652e-25; //cm^2
    float mecsq = 510.998; //keV
    double Mpc_to_cm_cube = 2.938e73;
    float beta = 0.75; //beta_T in Waizmann & Bartelmann 2009
    float Mvir = HuKravtsov(z, M200, rhocrit, Delta_vir, h); //Solar mass
    float Rvir = pow(3. * Mvir / Delta_vir / rhocrit / 4. / M_PI, 1./3.); //Mpc
    double Ne = 1.76e56 * Mvir; //Number of electrons
    //printf("Omega_Mz %f z %f Delta_vir %f Mvir %.2e h %f\n", Omega_Mz, z, Delta_vir, Mvir, h);
    float Oz = Omega_M * Delta_vir / Omega_Mz;
    float Mnorm = Mvir * h / 1e15; 
    //printf("Oz %f Mnorm %f 1/3 %f\n", pow(Oz, 1/3.), pow(Mnorm, (1/3.)), 1/3);
    //float Rvir_1 = 9.5103 * pow(Oz, -1/3.) * pow(Mnorm, 1/3.) / h / (1 + z); //Mpc
    //printf("Rvir %f, Rvir_1 %f \n", Rvir, Rvir_1);
    float kbt = (1 + z)  * pow(Oz, 1/3.) * pow(Mnorm, 2/3.) / beta; //keV
    float xiz = 0.14 * pow(1 + z, 1/5.); //Waizmann & Bartelmann Eq. 17
    float rc = Rvir * xiz; //Mpc
    float thetac = rc / ang_dia; //rad
    double ne0 = Ne / 4. / M_PI / pow(rc, 3.) / (1 / xiz + atan(xiz) - M_PI / 2.); // Number/Mpc^3
    //float ner = ne0  * (1 + (r / rc)**2.)**(-1.5 * 0.86)/3.086e24**3;
    ne0 = ne0 / Mpc_to_cm_cube; //Number/cm^3
    float y0 = 2 * kbt * sigma_T * rc * ne0 / mecsq; 

    //printf("ti %f tmax %f sqrt %f\n", ti, tmax, sqrt(tmax*tmax-ti*ti));
    gsl_integration_workspace * w = gsl_integration_workspace_alloc (1e6);

    double result, error;
    float params[3] = {ti, xiz, rc / ang_dia};
    gsl_function F;
    F.function = &yproj_func_Beta_WB;
    F.params = &params;
    double Rvir_rad = 1* Rvir / ang_dia;
    //printf("Rvir %f Rvir_rad %.2e, multi %.2e \n ", Rvir, Rvir_rad, sqrt(Rvir_rad * Rvir_rad - ti * ti));
    //printf("ti %f Rvir %f Rvir_rad %.2e, multi %.2e \n ", ti, Rvir, Rvir_rad, sqrt(Rvir_rad * Rvir_rad - ti * ti));
    if(Rvir_rad > ti){
        printf("ti %f Rvir_rad %f \n", ti, Rvir_rad);
        if (sqrt(Rvir_rad * Rvir_rad - ti * ti) <= Rvir_rad){
        printf("ti %f Rvir %f Rvir_rad %.2e, multi %.2e \n ", ti, Rvir, Rvir_rad, sqrt(Rvir_rad * Rvir_rad - ti * ti));
        gsl_integration_qags (&F, 0, sqrt(Rvir_rad*Rvir_rad-ti*ti), 0, 1e-6, 1e6 , w, &result, &error);
        printf("result %.2f \n", result);
    }}
    else {
        result = 0;
    }

    float this_ysz = y0 * 2.0 * result;
    //printf("kbt %f y0 %.2e result %.e this_y %.e \n", kbt, y0, result, this_ysz);
    //float p_e = kbt * ner; //I think 0.518 is included in the equation, therefore I don't need to mutiply this by 0.518
    //float this_ysz = p_e * 0.00402; //See Battaglia profile
    return this_ysz;
}



