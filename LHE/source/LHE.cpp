#include "Structure.h"
#include <cstdio>
#include <cstring>
#include <string>
#include "LHE.h"
#include <sys/stat.h>
#include <sys/types.h>



LHE::LHE(LHE_Config config,LHE_Event_Header h,LHE_Init i) {
    this->config = config;
    this->h = h;
    this->i = i;
}

void LHE::write_config(const std::string& name){

    // 1) Create Results/ directory if it does not exist
    mkdir("Results", 0777);

    // 2) Create Results/<name>/ directory
    std::string dir = "Results/" + name;
    mkdir(dir.c_str(), 0777);

    // 3) Build full output file path
    std::string filename = dir + "/" + name + "_LHE.txt";

    FILE* file = std::fopen(filename.c_str(), "w");
    if (!file) {
        perror("fopen");
        return;
    }

    std::fprintf(file, "<LesHouchesEvents version=\"3.0\">\n");
	std::fprintf(file, "<header>\n");
	std::fprintf(file, "<configPhysics>\n");
	switch (config.model) {
        case 100:
            std::fprintf(file, "# Chosen model : Photon exchange\n");
            break;
        case 010:
            std::fprintf(file, "# Chosen model : Z-boson exchange\n");
            break;
        case 001:
            std::fprintf(file, "# Chosen model : Higgs boson exchange\n");
            break;
        case 110:
            std::fprintf(file, "# Chosen model : Photon + Z-boson exchange\n");
            break;
        case 101:
            std::fprintf(file, "# Chosen model : Photon + Higgs boson exchange\n");
            break;
        case 011:
            std::fprintf(file, "# Chosen model : Z-boson + Higgs boson exchange\n");
            break;
        case 111:
            std::fprintf(file, "# Chosen model : Photon + Z-boson + Higgs boson exchange\n");
            break;
    }
    std::fprintf(file, "# Collision energy : %lf\n", config.sqrt_s);
    std::fprintf(file, "############### MASSES and WIDTHS ##############\n");
    std::fprintf(file, "Electron mass     : %lf GeV\n", config.m_e);
    std::fprintf(file, "Top quark mass    : %lf GeV\n", config.m_top);
    std::fprintf(file, "Z-boson mass      : %lf GeV\n", config.m_z);
    std::fprintf(file, "Higgs boson mass  : %lf GeV\n", config.m_h);
    std::fprintf(file, "W-boson mass      : %lf GeV\n", config.m_w);
    std::fprintf(file, "Z-boson width     : %lf GeV\n", config.gamma_z);
    std::fprintf(file, "Higgs boson width : %lf GeV\n", config.gamma_h);
    std::fprintf(file, "################################################\n");
    std::fprintf(file, "# Renormalisation scale : %lf\n", config.sqrt_s); 
    std::fprintf(file, "# sin^2(theta_W) : %lf\n", config.sin2_theta_W);
    std::fprintf(file, "# QED coupling constant : %lf\n", config.alphaQED);
    std::fprintf(file, "# Strong coupling constant : %lf \n", config.alphaS);
    std::fprintf(file, "# Number of color : %lf\n", config.Nc);
	std::fprintf(file, "</configPhysics>\n");
	std::fprintf(file, "<configGenerator>\n");
	std::fprintf(file, "# Number of generated events : %lf\n", double(config.n_events));
	std::fprintf(file, "# Seed                       : %lf\n", double(config.seed));
	std::fprintf(file, "# Weight strategy            : %lf\n", double(config.weight_strategy));
    std::fprintf(file, "# Number of Threads          : %lf\n", double(config.nb_threads));
	std::fprintf(file, "</configGenerator>\n");
	std::fprintf(file, "</header>\n");
	std::fclose(file);
}

void LHE::write_init(const std::string& name) {

    std::string dir = "Results/" + name;

    // Build full output file path 
    std::string filename = dir + "/" + name + "_LHE.txt";

    FILE* file = std::fopen(filename.c_str(), "a"); // 'a' because the file was already created in write_config 
    if (!file) {
        perror("fopen");
        return;
    }

    std::fprintf(file, "<init>\n");
	std::fprintf(file, "%d %d %e %e %d %d %d %d %d %d\n",i.PDG_ID_1, i.PDG_ID_2, i.beam_energy_1, i.beam_energy_2,
		i.Family_ID_1, i.Family_ID_2, i.Set_ID_1, i.Set_ID_2, i.weight_strategy, i.number_physical_processes);
	std::fprintf(file, "%e %e %e %d\n", i.Xsection_value, i.Xsection_error, i.maximum_weight_value, i.Process_ID);
	std::fprintf(file, "<generator name='Madgraph_Eco+' version=1.0> </generator>\n");
	std::fprintf(file, "</init>\n");
    std::fclose(file);

}

void LHE::write_event_para(FILE* file,LHE_Particle e,LHE_Particle ebar,LHE_Particle t,LHE_Particle tbar) {

    std::fprintf(file, "<event>\n");
	std::fprintf(file, "%d %d %e %e %e %e\n", h.number_generated_particle, h.Process_ID, h.event_weight,
		h.renorm_scale, h.alphaS, h.alphaQED);
	std::fprintf(file, "%d %d    %d    %d    %d    %d %e %e %e %e %e %e %d\n",
		e.beam_id, e.part_status, e.mother_1, e.mother_2, e.color_flow_1, e.color_flow_2, e.px, e.py, e.pz,
		e.E, e.mass, e.mean_distance, e.helicity);
	std::fprintf(file, "%d %d    %d    %d    %d    %d %e %e %e %e %e %e %d\n",
		ebar.beam_id, ebar.part_status, ebar.mother_1, ebar.mother_2, ebar.color_flow_1, ebar.color_flow_2, ebar.px, ebar.py, ebar.pz,
		ebar.E, ebar.mass, ebar.mean_distance, ebar.helicity);
	std::fprintf(file, "%d %d    %d    %d    %d    %d %e %e %e %e %e %e %d\n",
		t.beam_id, t.part_status, t.mother_1, t.mother_2, t.color_flow_1, t.color_flow_2, t.px, t.py, t.pz,
		t.E, t.mass, t.mean_distance, t.helicity);
	std::fprintf(file, "%d %d    %d    %d    %d    %d %e %e %e %e %e %e %d\n",
		tbar.beam_id, tbar.part_status, tbar.mother_1, tbar.mother_2, tbar.color_flow_1, tbar.color_flow_2, tbar.px, tbar.py, tbar.pz,
		tbar.E, tbar.mass, tbar.mean_distance, tbar.helicity);
	std::fprintf(file, "</event>\n");


}

void LHE::write_end(const std::string& name) {
    std::string filename = "Results/" + name + "/" + name + "_LHE.txt";

    FILE* file = std::fopen(filename.c_str(), "a");
    if (!file) {
        perror("fopen");
        return;
    }

    std::fprintf(file, "</LesHouchesEvents>");
    std::fclose(file);  
}
