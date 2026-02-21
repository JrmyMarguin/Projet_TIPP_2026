#pragma once
#include "Structure.h"

struct XSectionContrib {
    double gamma   = 0.0;
    double Z       = 0.0;
    double H       = 0.0;
    double gammaZ  = 0.0;
    double gammaH  = 0.0;
    double ZH      = 0.0;
};


class Physics {
/* Class containing all methods involving Physics calculations*/
private:
    LHE_Config config;

public:

    Physics(LHE_Config config);
    double ds(double theta); /* differential cross section*/
    XSectionContrib individual_ds(double theta); /* Individual contribution of a specific diagram(or interference term) to the differential cross-section */
    double total_x_section(double dtheta); /* total cross section */
    XSectionContrib total_x_section_contrib(double dtheta); /* Individual contribution to the total cross-section */
    double renormalisation(int k,double scale,double scale0,double alpha0,double* mf, double* Qf,double Nc); /* Computation of coupling constants after renormalization */
    double M_gamma_square(double h_e, double h_ebar,double h_t, double h_tbar,double theta); /* CMatrix Element squared / Interference Terms used directly in Monte Carlo sampling of helicities */
    double M_Z_square(double h_e, double h_ebar,double h_t, double h_tbar,double theta);
    double M_h_square(double h_e, double h_ebar,double h_t, double h_tbar);
    double M_gamma_H_interf(double h_e, double h_ebar,double h_t, double h_tbar,double theta);
    double M_Z_H_interf(double h_e, double h_ebar,double h_t, double h_tbar,double theta);
    double M_gamma_Z_interf(double h_e, double h_ebar,double h_t, double h_tbar,double theta);
};