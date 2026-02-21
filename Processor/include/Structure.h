#pragma once

#include <cmath>

struct LHE_Particle {

    int beam_id = 0;
    int part_status = 0;

    int mother_1 = 0;
    int mother_2 = 0;

    int color_flow_1 = 0;
    int color_flow_2 = 0;

    double px = 0;
    double py = 0;
    double pz = 0;
    double E = 0;
    double mass = 0;

    double mean_distance = 0;
    int helicity = 0;

    /* Observable Methods */

    static double inv_mass(const LHE_Particle& a, const LHE_Particle& b) {
        return std::sqrt((a.E + b.E)*(a.E+b.E) - (a.px + b.px)*(a.px + b.px) - (a.py + b.py)*(a.py + b.py) - (a.pz + b.pz)*(a.pz + b.pz));
    }
    
    double p_abs() const {
        return std::sqrt(px*px + py*py + pz*pz);
    }

    double pt() const {
        return std::sqrt(px*px + py*py);
    }

    double eta() const {
        double pabs = p_abs();
        return 0.5*std::log((pabs+pz)/(pabs-pz));
    }

    double costheta() const {
        double pabs = p_abs();
        return (pz/pabs);
    }
    
};

struct LHE_Init {

    int PDG_ID_1 = 0;
    int PDG_ID_2 = 0;

    double beam_energy_1 = 0;
    double beam_energy_2 = 0;

    int Family_ID_1 = 0;
    int Family_ID_2 = 0;

    int Set_ID_1 = 0;
    int Set_ID_2 = 0;

    int weight_strategy = 0; 
    int number_physical_processes = 0;

    double Xsection_value = 0;
    double Xsection_error = 0;
    double maximum_weight_value = 0;
    int Process_ID = 0;
};

struct LHE_Event_Header {

    int number_generated_particle = 0;
    int Process_ID = 0;

    double event_weight = 0;
    double renorm_scale = 0;
    double alphaS = 0;
    double alphaQED = 0;
};



