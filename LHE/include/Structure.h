#pragma once

/* Data structures used in all classes, allowing fast and organized handling of all input and output parameters */

struct LHE_Particle {

    int beam_id;
    int part_status;

    int mother_1;
    int mother_2;

    int color_flow_1;
    int color_flow_2;

    double px;
    double py;
    double pz;
    double E;
    double mass;

    double mean_distance;
    int helicity;
    

};

struct LHE_Init {

    int PDG_ID_1;
    int PDG_ID_2;

    double beam_energy_1;
    double beam_energy_2;

    int Family_ID_1;
    int Family_ID_2;

    int Set_ID_1;
    int Set_ID_2;

    int weight_strategy; 
    int number_physical_processes;

    double Xsection_value;
    double Xsection_error;
    double maximum_weight_value;
    int Process_ID;

    
};

struct LHE_Event_Header {

    int number_generated_particle;
    int Process_ID;

    double event_weight;
    double renorm_scale;
    double alphaS;
    double alphaQED;

};

struct LHE_Config {

    /* Loaded from config.yaml*/
    int model;
    double sqrt_s;
    double m_e;
    double m_top;
    double m_z;
    double m_h;
    double m_w;
    double gamma_z;
    double gamma_h;
    double renorm_scale;
    double alphaS;
    double alphaQED;
    double Nc;
    double m_mu;
    double m_tau;
    double m_u;
    double m_d;
    double m_s;
    double m_c;
    double m_b;
    int n_events;
    int seed;
    int weight_strategy;
    int nb_threads;

    /* Activation of different Feynmann Diagramms + interference terms */
    int k1;
    int k2;
    int k3;
    int k4;
    int k5;
    int k6;
    
    /* Computed Parameters */
    double sin2_theta_W;
    double sinw;
    double cosw;
    double e;
    double gW;
    double K_gamma;
    double K_z;
    double K_h;
    double cem;
    double cep;
    double ctm;
    double ctp;
    double be;
    double bt;
    double g_gamma;
    double gh;
    double gz;
    double gzh;
    double gzgamma;

};





