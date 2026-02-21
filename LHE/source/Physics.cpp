#include <cmath>
#include "Physics.h"
#include "Structure.h"
#include <cstdio>


Physics::Physics(LHE_Config config) {

    this->config  = config;
}

/* Compute Differential Cross-Section / Cross-Section */
double Physics::ds(double theta)
{
    double s = config.sqrt_s*config.sqrt_s;
    double c = std::cos(theta);

    double Nc = config.Nc;
    double mZ = config.m_z;
    double mH = config.m_h;
    double rH = config.gamma_h;
    double rZ = config.gamma_z;
    double cep = config.cep;
    double cem = config.cem;
    double ctp = config.ctp;
    double ctm = config.ctm;
    double bt = config.bt;
    double be = config.be;
    double Kz = config.K_z;
    double Kh = config.K_h;
    double Kgamma = config.K_gamma;

    double term1 =
        config.k1 * Nc * Kgamma * Kgamma *
        (3.0 - bt*bt - be*be + bt*bt * be*be * c*c);

    double denomZ = (s - mZ*mZ)*(s - mZ*mZ) + (mZ*rZ)*(mZ*rZ);

    double term2 =
        config.k2 * (Nc * Kz * Kz / 4.0) *
        (s*s / denomZ) *
        (
            (cep*cep + cem*cem) * (ctp*ctp + ctm*ctm)
            * (1.0 + bt*bt * be*be * c*c)
          - 2.0 * (cep*cep - cem*cem)
                 * (ctp*ctp - ctm*ctm)
                 * bt * be * c

          + (cep*cep + cem*cem) * ctp * ctm
            * (1.0 - bt*bt) * (1.0 + be*be)

          + (ctp*ctp + ctm*ctm) * cep * cem
            * (1.0 + bt*bt) * (1.0 - be*be)

          + 4.0 * cep * cem * ctp * ctm
            * (1.0 - bt*bt) * (1.0 - be*be)
        );

    double denomH = (s - mH*mH)*(s - mH*mH) + (mH*rH)*(mH*rH);

    double term3 =
        config.k3 * (Nc * Kh * Kh / 16.0) *
        bt*bt * be*be * (1.0 - bt*bt) * (1.0 - be*be) *
        (s*s*s*s / denomH);

    double term4 =
        -config.k4 * (Nc * Kgamma * Kz / 2.0) *
        (s * (s - mZ*mZ) / denomZ) *
        (
            (cep + cem) * (ctp + ctm)
            * (3.0 - be*be - bt*bt + bt*bt * be*be * c*c)

          - 2.0 * (cep - cem) * (ctp - ctm) * bt * be * c
        );

    double term5 =
        config.k5 * (Nc * Kh * Kgamma / 2.0) *
        bt * be * (1.0 - bt*bt) * (1.0 - be*be) *
        s*s * (s - mH*mH) * c / denomH;

    double term6 =
        config.k6 * Nc * Kh * Kz *
        (cem + cep) * (ctm + ctp) * bt * be *
        (1.0 - bt*bt) * (1.0 - be*be) *
        s*s*s * ((s - mZ*mZ)*(s - mH*mH) + mH*mZ*rH*rZ) * c /
        (8.0 * denomH * denomZ);

    double pref =
        (2.0 / (config.sqrt_s * be))
        * (bt / (32.0 * M_PI * M_PI))
        / (4.0*config.sqrt_s);

    return pref * (term1 + term2 + term3 + term4 + term5 + term6);
}

XSectionContrib Physics::individual_ds(double theta)
{
    double s = config.sqrt_s*config.sqrt_s;
    double c = std::cos(theta);

    double Nc = config.Nc;
    double mZ = config.m_z;
    double mH = config.m_h;
    double rH = config.gamma_h;
    double rZ = config.gamma_z;
    double cep = config.cep;
    double cem = config.cem;
    double ctp = config.ctp;
    double ctm = config.ctm;
    double bt = config.bt;
    double be = config.be;
    double Kz = config.K_z;
    double Kh = config.K_h;
    double Kgamma = config.K_gamma;

    double term1 =
        config.k1 * Nc * Kgamma * Kgamma *
        (3.0 - bt*bt - be*be + bt*bt * be*be * c*c);

    double denomZ = (s - mZ*mZ)*(s - mZ*mZ) + (mZ*rZ)*(mZ*rZ);

    double term2 =
        config.k2 * (Nc * Kz * Kz / 4.0) *
        (s*s / denomZ) *
        (
            (cep*cep + cem*cem) * (ctp*ctp + ctm*ctm)
            * (1.0 + bt*bt * be*be * c*c)
          - 2.0 * (cep*cep - cem*cem)
                 * (ctp*ctp - ctm*ctm)
                 * bt * be * c

          + (cep*cep + cem*cem) * ctp * ctm
            * (1.0 - bt*bt) * (1.0 + be*be)

          + (ctp*ctp + ctm*ctm) * cep * cem
            * (1.0 + bt*bt) * (1.0 - be*be)

          + 4.0 * cep * cem * ctp * ctm
            * (1.0 - bt*bt) * (1.0 - be*be)
        );

    double denomH = (s - mH*mH)*(s - mH*mH) + (mH*rH)*(mH*rH);

    double term3 =
        config.k3 * (Nc * Kh * Kh / 16.0) *
        bt*bt * be*be * (1.0 - bt*bt) * (1.0 - be*be) *
        (s*s*s*s / denomH);

    double term4 =
        -config.k4 * (Nc * Kgamma * Kz / 2.0) *
        (s * (s - mZ*mZ) / denomZ) *
        (
            (cep + cem) * (ctp + ctm)
            * (3.0 - be*be - bt*bt + bt*bt * be*be * c*c)

          - 2.0 * (cep - cem) * (ctp - ctm) * bt * be * c
        );

    double term5 =
        config.k5 * (Nc * Kh * Kgamma / 2.0) *
        bt * be * (1.0 - bt*bt) * (1.0 - be*be) *
        s*s * (s - mH*mH) * c / denomH;

    double term6 =
        config.k6 * Nc * Kh * Kz *
        (cem + cep) * (ctm + ctp) * bt * be *
        (1.0 - bt*bt) * (1.0 - be*be) *
        s*s*s * ((s - mZ*mZ)*(s - mH*mH) + mH*mZ*rH*rZ) * c /
        (8.0 * denomH * denomZ);

    double pref =
        (2.0 / (config.sqrt_s * be))
        * (bt / (32.0 * M_PI * M_PI))
        / (4.0*config.sqrt_s);

    XSectionContrib result;

    result.gamma   = pref * term1;
    result.Z       = pref * term2;
    result.H       = pref * term3;
    result.gammaZ  = pref * term4;
    result.gammaH  = pref * term5;
    result.ZH      = pref * term6;

    return result;
}



double Physics::total_x_section(double dtheta) {
    double y = 0.0;
    double theta = 0.0;
    while(theta < M_PI) {
        y += 2*M_PI*ds(theta)*sin(theta)*dtheta;
        theta += dtheta;
    }
    return y;
}

XSectionContrib Physics::total_x_section_contrib(double dtheta)
{
    XSectionContrib total;

    double theta = 0.0;

    while(theta < M_PI)
    {
        XSectionContrib contrib = individual_ds(theta);

        double jacobian = 2.0 * M_PI * std::sin(theta) * dtheta;

        total.gamma   += contrib.gamma  * jacobian;
        total.Z       += contrib.Z      * jacobian;
        total.H       += contrib.H      * jacobian;
        total.gammaZ  += contrib.gammaZ * jacobian;
        total.gammaH  += contrib.gammaH * jacobian;
        total.ZH      += contrib.ZH     * jacobian;

        theta += dtheta;
    }

    return total;
}


/*1-loop renormalization using RGE equations */
double Physics::renormalisation(int k,double scale,double scale0,double alpha0,double* mf, double* Qf,double Nc) {

    double alphaRenorm;

    /*QED (k=0)*/
    if(k==0) {
        double sum = 0.0;

        /*Lepton*/
        for(int i=0; i < 3; i++) {
            if(mf[i]<scale) {
                sum += Qf[i]*Qf[i];
            }
        }
        /*Quark*/
        for(int i=3; i < 9; i++) {
            if(mf[i]<scale) {
                sum += Nc*Qf[i]*Qf[i];
            }
        }

        double b = 2*sum/(3*M_PI);

        alphaRenorm = alpha0/(1-b*alpha0*std::log(scale/scale0));

        return alphaRenorm;
    }

    /*QCD (k=1)*/
    if(k==1) {
        double Nf = 0;

        /*Quark*/
        for(int i=3; i<9; i++) {
            if(mf[i] < scale) Nf++;
        }

        double b = (11-(2.0*Nf/3.0))/(2*M_PI);

        alphaRenorm = alpha0/(1+b*alpha0*std::log(scale/scale0));

        return alphaRenorm;
    }

    return alpha0;
}


/* Compute Matrix Elements */
/* ========================================================= */
/* ===================  PHOTON  ============================ */
/* ========================================================= */

double Physics::M_gamma_square(double h_e,double h_ebar,
                               double h_t,double h_tbar,
                               double theta)
{
    double be=config.be;
    double bt=config.bt;
    double g=config.g_gamma;

    double c=std::cos(theta);
    double s2=std::sin(theta)*std::sin(theta);

    if(h_e==1 && h_ebar==1)
    {
        if(h_t==h_tbar)
            return g*(1-be*be)*(1-bt*bt)*c*c;

        if(h_t==-h_tbar)
            return g*(1-be*be)*s2;
    }

    if(h_e==-1 && h_ebar==-1)
    {
        if(h_t==h_tbar)
            return g*(1-be*be)*(1-bt*bt)*c*c;

        if(h_t==-h_tbar)
            return g*(1-be*be)*s2;
    }

    if(h_e==1 && h_ebar==-1)
    {
        if(h_t==h_tbar)
            return g*(1-bt*bt)*s2;

        if(h_t==1 && h_tbar==-1)
            return g*(1+c)*(1+c);

        if(h_t==-1 && h_tbar==1)
            return g*(1-c)*(1-c);
    }

    if(h_e==-1 && h_ebar==1)
    {
        if(h_t==h_tbar)
            return g*(1-bt*bt)*s2;

        if(h_t==-1 && h_tbar==1)
            return g*(1+c)*(1+c);

        if(h_t==1 && h_tbar==-1)
            return g*(1-c)*(1-c);
    }

    return 0;
}


/* ========================================================= */
/* ======================  Z  ============================== */
/* ========================================================= */

double Physics::M_Z_square(double h_e,double h_ebar,
                           double h_t,double h_tbar,
                           double theta)
{
    double be=config.be;
    double bt=config.bt;
    double cem=config.cem;
    double cep=config.cep;
    double ctm=config.ctm;
    double ctp=config.ctp;
    double gz=config.gz;

    double c=std::cos(theta);
    double s2=std::sin(theta)*std::sin(theta);

    if(h_e==1 && h_ebar==1 && h_t==1 && h_tbar==1)
        return gz*(1-be*be)*(1-bt*bt)*
               std::pow((c+1)*(cep*ctp+cem*ctm)+(c-1)*(cep*ctm+cem*ctp),2);

    if(h_e==-1&&h_ebar==-1&&h_t==-1&&h_tbar==-1)
        return gz*(1-be*be)*(1-bt*bt)*
               std::pow((c+1)*(cep*ctp+cem*ctm)+(c-1)*(cep*ctm+cem*ctp),2);

    if(h_e==1&&h_ebar==1&&h_t==-1&&h_tbar==-1)
        return gz*(1-be*be)*(1-bt*bt)*
               std::pow((c-1)*(cep*ctm+cem*ctp)+(c+1)*(cep*ctp+cem*ctm),2);

    if(h_e==-1&&h_ebar==-1&&h_t==1&&h_tbar==1)
        return gz*(1-be*be)*(1-bt*bt)*
               std::pow((c-1)*(cep*ctm+cem*ctp)+(c+1)*(cep*ctp+cem*ctm),2);

    if(h_e==1&&h_ebar==-1&&(h_t==h_tbar))
        return gz*s2*(1-bt*bt)*std::pow(ctp+ctm,2)*
               std::pow((1+be)*cep+(1-be)*cem,2);

    if(h_e==-1&&h_ebar==1&&(h_t==h_tbar))
        return gz*s2*(1-bt*bt)*std::pow(ctp+ctm,2)*
               std::pow((1+be)*cem+(1-be)*cep,2);

    if(h_e==1&&h_ebar==-1&&h_t==1&&h_tbar==-1)
        return gz*(1+c)*(1+c)*
               std::pow((1+be)*cep+(1-be)*cem,2)*
               std::pow((1+bt)*ctm+(1-bt)*ctp,2);

    if(h_e==-1&&h_ebar==1&&h_t==-1&&h_tbar==1)
        return gz*(1+c)*(1+c)*
               std::pow((1+be)*cem+(1-be)*cep,2)*
               std::pow((1+bt)*ctp+(1-bt)*ctm,2);

    if(h_e==-1&&h_ebar==1&&h_t==1&&h_tbar==-1)
        return gz*(1-c)*(1-c)*
               std::pow((1+be)*cem+(1-be)*cep,2)*
               std::pow((1+bt)*ctm+(1-bt)*ctp,2);

    if(h_e==1&&h_ebar==-1&&h_t==-1&&h_tbar==1)
        return gz*(1-c)*(1-c)*
               std::pow((1+be)*cep+(1-be)*cem,2)*
               std::pow((1+bt)*ctp+(1-bt)*ctm,2);

    if(h_e==1&&h_ebar==1&&h_t==1&&h_tbar==-1)
        return gz*s2*(1-be*be)*(cep+cem)*(cep+cem)*((1+bt)*ctm+(1-bt)*ctp)*((1+bt)*ctm+(1-bt)*ctp);
        
    if(h_e==-1&&h_ebar==-1&&h_t==1&&h_tbar==-1)
        return gz*s2*(1-be*be)*(cep+cem)*(cep+cem)*((1+bt)*ctm+(1-bt)*ctp)*((1+bt)*ctm+(1-bt)*ctp);

    return 0;
}


/* ========================================================= */
/* ====================  HIGGS  ============================ */
/* ========================================================= */

double Physics::M_h_square(double h_e,double h_ebar,
                           double h_t,double h_tbar)
{
    double be=config.be;
    double bt=config.bt;
    double gh=config.gh;

    if((h_e==h_ebar)&&(h_t==h_tbar))
        return gh*bt*bt*be*be*(1-bt*bt)*(1-be*be);

    return 0;
}


/* ========================================================= */
/* ===================  GAMMA-Z  =========================== */
/* ========================================================= */

double Physics::M_gamma_Z_interf(double h_e,double h_ebar,
                                 double h_t,double h_tbar,
                                 double theta)
{
    double be = config.be;
    double bt = config.bt;
    double cem = config.cem;
    double cep = config.cep;
    double ctm = config.ctm;
    double ctp = config.ctp;
    double g = config.gzgamma;

    double c = std::cos(theta);
    double s2 = 1.0 - c*c;

    // -------- h_e = - , h_ebar = + --------

    if(h_e==-1 && h_ebar==1 && h_t==-1 && h_tbar==1)
        return g*(1+c)*(1+c)*
               ((1+be)*cem + (1-be)*cep)*
               ((1+bt)*ctp + (1-bt)*ctm);

    if(h_e==-1 && h_ebar==1 && h_t==1 && h_tbar==-1)
        return g*(1-c)*(1-c)*
               ((1+be)*cem + (1-be)*cep)*
               ((1+bt)*ctm + (1-bt)*ctp);

    if(h_e==-1 && h_ebar==1 && h_t==h_tbar)
        return g*s2*(1-bt*bt)*(ctp+ctm)*
               ((1+be)*cem + (1-be)*cep);

    // -------- h_e = + , h_ebar = - --------

    if(h_e==1 && h_ebar==-1 && h_t==1 && h_tbar==-1)
        return g*(1+c)*(1+c)*
               ((1+be)*cep + (1-be)*cem)*
               ((1+bt)*ctm + (1-bt)*ctp);

    if(h_e==1 && h_ebar==-1 && h_t==-1 && h_tbar==1)
        return g*(1-c)*(1-c)*
               ((1+be)*cep + (1-be)*cem)*
               ((1+bt)*ctp + (1-bt)*ctm);

    if(h_e==1 && h_ebar==-1 && h_t==h_tbar)
        return g*s2*(1-bt*bt)*(ctp+ctm)*
               ((1+be)*cep + (1-be)*cem);

    // -------- h_e = h_ebar --------

    if(h_e==h_ebar && h_t==-1 && h_tbar==1)
        return g*s2*(1-be*be)*
               ( (cep+cem) )*
               ((1+bt)*ctp + (1-bt)*ctm);

    if(h_e==h_ebar && h_t==1 && h_tbar==-1)
        return g*s2*(1-be*be)*
               ( (cep+cem) )*
               ((1+bt)*ctm + (1-bt)*ctp);

    if(h_e==h_ebar && h_t==h_tbar)
        return g*c*(1-be*be)*(1-bt*bt)*
               (cep+cem)*(ctp+ctm);

    return 0.0;
}




/* ========================================================= */
/* ===================  GAMMA-H  =========================== */
/* ========================================================= */

double Physics::M_gamma_H_interf(double h_e,double h_ebar,
                                 double h_t,double h_tbar,
                                 double theta)
{
    double be=config.be;
    double bt=config.bt;
    double ggh=config.gh;

    double c=std::cos(theta);

    if((h_e==h_ebar)&&(h_t==h_tbar))
        return ggh*bt*bt*be*be*(1-bt*bt)*(1-be*be)*c;

    return 0;
}


/* ========================================================= */
/* ======================  Z-H  ============================ */
/* ========================================================= */

double Physics::M_Z_H_interf(double h_e,double h_ebar,
                             double h_t,double h_tbar,
                             double theta)
{
    double cem=config.cem;
    double cep=config.cep;
    double ctm=config.ctm;
    double ctp=config.ctp;
    double gzh=config.gzh;

    double c=std::cos(theta);

    if((h_e==1&&h_ebar==1&&h_t==1&&h_tbar==1)||
       (h_e==-1&&h_ebar==-1&&h_t==-1&&h_tbar==-1))
        return gzh*(c*(cep+cem)*(ctp+ctm)
                   +(cep-cem)*(ctp-ctm));

    if((h_e==1&&h_ebar==1&&h_t==-1&&h_tbar==-1)||
       (h_e==-1&&h_ebar==-1&&h_t==1&&h_tbar==1))
        return gzh*(c*(cep+cem)*(ctp+ctm)
                   -(cep-cem)*(ctp-ctm));

    return 0;
}

