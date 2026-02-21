#pragma once
#include <string>
#include <cstdio>
#include "Structure.h"

class LHE {

    /* Class for writing/exporting generated events data to an LHE file format */
    private:
        LHE_Config config;
        LHE_Event_Header h;
        LHE_Init i;
    
    public:
        LHE(LHE_Config config,LHE_Event_Header h,LHE_Init i);
        void write_config(const std::string& name);
        void write_init(const std::string& name);
        void write_event_para(FILE* file,LHE_Particle e,LHE_Particle ebar,LHE_Particle t,LHE_Particle tbar);
        void write_end(const std::string& name);

};
