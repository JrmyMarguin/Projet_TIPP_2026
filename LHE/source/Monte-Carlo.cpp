#include <cstdio>
#include <random>
#include <cmath>
#include "Physics.h"
#include "Structure.h"
#include "Monte-Carlo.h"
#include "LHE.h"
#include <map>
#include <utility>
#include <array>
#include <ostream>
#include <iostream>
#include <omp.h>
#include <filesystem>



void MonteCarlo(char* name,LHE lhe,LHE_Config config,Physics P) {
    
    double s = config.sqrt_s*config.sqrt_s;

    /////////////////////////
    /*Determination of wmax*/
    /////////////////////////

    int Nscan = 50000;
    double wmax = 0;
    double wmin = 0;

    for(int i = 0; i < Nscan; i++) {
        double theta = M_PI * i / Nscan;  
        double w = P.ds(theta)*sin(theta);

        if(w > wmax)
            wmax = w;
        if(w < wmin) // in case of w < 0
            wmin = w;
    }

    



    ////////////////////////////////
    /*Monte Carlo event Generator*/
    //////////////////////////////

    
    int H[2] = {-1,1}; /* helicity values */
 
    
    double sum = 0.0;
    double sum2 = 0.0;
    int N = 0;

    omp_set_num_threads(config.nb_threads);

    printf("Calculation start...\n");
   
    #pragma omp parallel reduction(+:sum,sum2,N)
    { 
        int tid = omp_get_thread_num(); /* Parallelization */
        int nthreads = omp_get_num_threads();
        int local_events = config.n_events / nthreads;
        if(tid == nthreads-1)
            local_events += config.n_events % nthreads;


        std::mt19937 gen(config.seed + tid); /*seed generator for each thread*/

        char local_name[256];
        sprintf(local_name,"%s_thread%d_event.txt",name,tid); /*each thread writes to its own .txt file, all files are merged at the end*/
        int local_event = 0;

        std::uniform_real_distribution<double> dist_theta(0.0, M_PI);
        std::uniform_real_distribution<double> dist_phi(0.0, 2.0 * M_PI);
        std::uniform_real_distribution<double> dist_u(0.0, wmax-wmin);

        std::map<std::array<int,4>,double> helicity_MC; /* store M² values for each helicity configuration  : map {h_e,h_ebar,h_t,h_tbar} -> Value */


        std::string filename = std::string("Results/") + name + "/" + local_name;
        FILE* file = std::fopen(filename.c_str(), "a");
        if(!file) std::cout << "ERROR : Writing in LHE event file" << std::endl;

        while(local_event < local_events) {

            /* Generator angle theta */
            double theta = dist_theta(gen);
            theta = M_PI-theta; /* Madgraph convention*/
            double phi = dist_phi(gen);
            double w = P.ds(theta)*sin(theta);
            double w_2 = w - wmin;
            double u = dist_u(gen);

            sum += 2*M_PI*P.ds(theta)*sin(theta);
            sum2 += 2*M_PI*P.ds(theta)*sin(theta)*2*M_PI*P.ds(theta)*sin(theta);

            // store data in structure for each type of Particle
            LHE_Particle e{};
            LHE_Particle ebar{};
            LHE_Particle t{};
            LHE_Particle tbar{};

            /* Store initial parameters */
            e.beam_id=11, e.part_status=-1, e.mother_1=0, e.mother_2=0, e.color_flow_1=0, e.color_flow_2=0;
            e.px=+0.0000000000e+00, e.py=+0.0000000000e+00, e.pz=config.sqrt_s/2, e.E=config.sqrt_s/2, e.mass=config.m_e;
            e.mean_distance=0.0000e+00;

            ebar.beam_id=-11, ebar.part_status=-1, ebar.mother_1=0, ebar.mother_2=0, ebar.color_flow_1=0, ebar.color_flow_2=0;
            ebar.px=-0.0000000000e+00, ebar.py=-0.0000000000e+00, ebar.pz=-config.sqrt_s/2, ebar.E=config.sqrt_s/2, ebar.mass=config.m_e;
            ebar.mean_distance=0.0000e+00;

            t.beam_id=6, t.part_status=1, t.mother_1=1, t.mother_2=2, t.color_flow_1=501, t.color_flow_2=0;
            t.mean_distance=0.0000e+00;

            tbar.beam_id=-6, tbar.part_status=1, tbar.mother_1=1, tbar.mother_2=2, tbar.color_flow_1=0, tbar.color_flow_2=501;
            tbar.mean_distance=0.0000e+00;

            /* Acceptance - rejection Method for helicity*/
            if(u<w_2) {
                local_event++; 
                helicity_MC.clear();


                /* Helicity sampling at fixed theta */
                

                for (int h_e : H){
                    for (int h_ebar : H){
                        for (int h_t : H){
                            for (int h_tbar : H){

                                double w_helicity =
                                    config.k1*P.M_gamma_square(h_e, h_ebar, h_tbar, h_t, theta) +
                                    config.k2*P.M_Z_square(h_e, h_ebar, h_tbar, h_t, theta) +
                                    config.k3*P.M_h_square(h_e, h_ebar, h_tbar, h_t) -
                                    config.k4*P.M_gamma_Z_interf(h_e, h_ebar, h_tbar, h_t, theta) +
                                    config.k5*P.M_gamma_H_interf(h_e, h_ebar, h_tbar, h_t, theta) +
                                    config.k6*P.M_Z_H_interf(h_e, h_ebar, h_tbar, h_t, theta);


                                
                                
                                helicity_MC[{h_e,h_ebar,h_t,h_tbar}] = (1.0/4.0)*w_helicity; /*Non polarized beams (average over inital helicities)*/
                               
                                                                                                                       


                            }
                        }
                    }
                }

                double w_helicity_tot = 0; /* Compute total helicity weight function */
                for (const std::pair<const std::array<int,4>,double>& x : helicity_MC){ /*  map {h_e,h_ebar,h_t,h_tbar} -> Value */
                    w_helicity_tot += x.second;

                }



                /* Generate helicity according to cumulative weights function*/
                std::uniform_real_distribution<double> dist_helicity(0.0, w_helicity_tot); 
                double r = dist_helicity(gen);

                double cumul = 0;
                std::array<int,4> selected_helicity;

                for (const std::pair<const std::array<int,4>,double>& x : helicity_MC){ 
                    cumul += x.second;

                    if (r<cumul){
                        selected_helicity = x.first; /*selected helicity*/
                        break;
                    }
                }

                /* Save helicity configuration */
                e.helicity = selected_helicity[0];
                ebar.helicity = selected_helicity[1];
                t.helicity = selected_helicity[2];
                tbar.helicity = selected_helicity[3];
                
                /* Final state kinematics*/
                t.mass = config.m_top;
                t.E = sqrt(s)/2;
                t.px = sqrt(s)*config.bt*sin(theta)*cos(phi)/2;
                t.py = sqrt(s)*config.bt*sin(theta)*sin(phi)/2;
                t.pz = sqrt(s)*config.bt*cos(theta)/2;

                tbar.mass = config.m_top;
                tbar.E = sqrt(s)/2;
                tbar.px = -sqrt(s)*config.bt*sin(theta)*cos(phi)/2;
                tbar.py = -sqrt(s)*config.bt*sin(theta)*sin(phi)/2;
                tbar.pz = -sqrt(s)*config.bt*cos(theta)/2;

                lhe.write_event_para(file,e,ebar,t,tbar); 
            }

            N++;

            
        }
        std::fclose(file);
    }

    /* Monte Carlo integration of differential cross section with statistical error */
    const double GEV2_TO_BARN = 3.89379e-4;

    double error_f = sum/N;
    double error_f2 = sum2/N;

    double dsigma = 1e12 * GEV2_TO_BARN * M_PI * sqrt(error_f2-error_f*error_f)/sqrt(N);
    double sigma = 1e12 * GEV2_TO_BARN * M_PI * error_f;

    printf("Total cross section (MC integrator) : %.10e +- %.10e pb\n",sigma,dsigma);

    /*Merge all parralelized files into one ouput */
    printf("Merging files together...\n");

    std::string filename = std::string("Results/") + name + "/" + name + "_LHE.txt";

    FILE* final = fopen(filename.c_str(),"a");

    int nthreads = omp_get_max_threads();

    for(int i=0;i<nthreads;i++)
    {
        char tmp[512];
        sprintf(tmp,"Results/%s/%s_thread%d_event.txt",name,name,i);


        FILE* f = fopen(tmp,"r");
        if(!f) continue;

        char buffer[4096];
        while(fgets(buffer,4096,f))
            fputs(buffer,final);

        fclose(f);
        remove(tmp);
    }

    fclose(final);

    printf("LHE file write...\n");
    lhe.write_end(name);

    
    // Copy .txt file to .lhe extension
    std::string filetxt = std::string("Results/") + name + "/" + name + "_LHE.txt";
    std::string filelhe = std::string("Results/") + name + "/" + name + "_LHE.lhe";
    try {
        std::filesystem::copy_file(filetxt.c_str(), filelhe.c_str(), std::filesystem::copy_options::overwrite_existing);
    } catch (std::filesystem::filesystem_error& e) {
        std::cerr << "Erreur : " << e.what() << std::endl;
    }
}

