#include <yaml-cpp/yaml.h>
#include <cstdio> 
#include "Monte-Carlo.h"
#include "LHE.h"
#include "Structure.h"
#include "Physics.h"
#include <ctime>
#include <cstring>
#include <omp.h>
#include <iostream>

/* sudo apt install libyaml-cpp-dev */


int main(){

    printf("\n");
    printf("***************************************************************************************************\n");
    printf("*                                   e+ e- -> t t̄ Monte Carlo                                      *\n");
    printf("*                                        Madgraph Eco+                                            *\n");
    printf("*                                                                                                 *\n");
    printf("*   e+  \\                /  t       e+ \\                /  t     e+ \\                /  t         *\n");
    printf("*        \\              /               \\              /             \\              /             *\n");
    printf("*         \\     γ      /                 \\     Z      /               \\     H      /              *\n");
    printf("*          o~~~~~~~~~~o         +         o~~~~~~~~~~o        +        o----------o               *\n");
    printf("*         /            \\                 /            \\               /            \\              *\n");
    printf("*        /              \\               /              \\             /              \\             *\n");
    printf("*   e-  /                \\  t̄       e- /                \\  t̄     e- /                \\  t̄         *\n");
    printf("*                                                                                                 *\n");
    printf("*                           SAUER Julien / GOGEL Yann / MARGUIN Jérémy                            *\n");
    printf("***************************************************************************************************\n\n");

    while (true){

        char name[100];
        char name_config[100];
        char path[200];

        printf("Experiment name : ");
        scanf("%100s", name);

        /* clear remaining \n in stdin */
        int c;
        while ((c = getchar()) != '\n' && c != EOF);

        printf("Configuration file to load (default : 'Config' | ' ') : ");
        fgets(name_config, sizeof(name_config), stdin);


        
        name_config[strcspn(name_config, "\n")] = 0;

        
        if (name_config[0] == '\0') {
            strcpy(name_config, "Config");
        }


        snprintf(path, sizeof(path),
                "LHE/configuration/%s.yaml", name_config);

        printf("Modify configuration file ? (y|n) ");
        char rep;
        scanf(" %c", &rep);

        if (rep == 'y') {
            printf("Opening the configuration file...\n");

            char cmd[512];
            snprintf(cmd, sizeof(cmd), "nano %s", path);
            system(cmd);   // open nano ; program waits until nano exists

            printf("Configuration save. Reloading...\n");
        }
        else if (rep == 'n') {
            printf("No modification applied to configuration file\n");
        }
        else {
            printf("ERROR : choose y ou n\n");
            return 0;
        }


        YAML::Node config = YAML::LoadFile(path);

        LHE_Config config_value {};

        /* Load physics parameters */
        config_value.model = config["physics"]["model"].as<int>();
        config_value.sqrt_s = config["physics"]["sqrt_s"].as<double>();
        config_value.m_e = config["physics"]["m_e"].as<double>();
        config_value.m_top = config["physics"]["m_top"].as<double>();
        config_value.m_z = config["physics"]["m_z"].as<double>();
        config_value.m_h = config["physics"]["m_h"].as<double>();
        config_value.m_w = config["physics"]["m_w"].as<double>();
        config_value.gamma_z = config["physics"]["gamma_z"].as<double>();
        config_value.gamma_h = config["physics"]["gamma_h"].as<double>();
        config_value.renorm_scale = config["physics"]["renorm_scale"].as<double>();
        config_value.alphaS = config["physics"]["alphaS"].as<double>();
        config_value.alphaQED = config["physics"]["alphaQED"].as<double>();
        config_value.Nc = config["physics"]["Nc"].as<double>();
        config_value.m_mu  = config["physics"]["m_mu"].as<double>();
        config_value.m_tau = config["physics"]["m_tau"].as<double>();
        config_value.m_u = config["physics"]["m_u"].as<double>();
        config_value.m_d = config["physics"]["m_d"].as<double>();
        config_value.m_s = config["physics"]["m_s"].as<double>();
        config_value.m_c = config["physics"]["m_c"].as<double>();
        config_value.m_b = config["physics"]["m_b"].as<double>();

        /* Load Generator parameters*/
        config_value.n_events = config["generator"]["n_events"].as<int>();
        config_value.seed = config["generator"]["seed"].as<int>();
        config_value.weight_strategy = config["generator"]["weight_strategy"].as<int>();
        config_value.nb_threads = config["generator"]["nb_threads"].as<int>();


        /* Validation checks */

        if(config_value.nb_threads > omp_get_max_threads() || config_value.nb_threads <= 0){
            printf("\nERROR : Invalid number of threads for this machine\n");
            return 0;
        }

        if (config_value.m_top < 0 || config_value.m_e < 0 || config_value.gamma_z < 0 || config_value.gamma_h < 0 || config_value.renorm_scale < 0 || config_value.sqrt_s < 0 ) {
            printf("\nERROR : Negative value detected (mass/width/energy scale) \n");
            return 0; 
        }

        if (config_value.m_mu < 0 || config_value.m_tau < 0 || config_value.m_u < 0 || config_value.m_d < 0 || config_value.m_s < 0 || config_value.m_c < 0 || config_value.m_b < 0) {
            printf("\nERROR : Negative mass value detected (quark/mu/tau) \n");
            return 0; 
        }

        if (config_value.sqrt_s < 2 * config_value.m_top) {
            printf("\nERROR : Insufficient energy to produce a top pair \n");
            return 0;
        }



        /* Process Selection (QED + Z ; QED ; Z ; Higgs ; QED + Z + Higgs etc...)*/
        int k1 = 0, k2 = 0, k3 = 0, k4 = 0, k5 = 0, k6 = 0;

        switch (config_value.model) {
            case 100:
                k1 = 1;
                break;

            case 010:
                k2 = 1;
                break;

            case 001:
                k3 = 1;
                break;

            case 110:
                k1 = 1; k2 = 1; k4 = 1;
                break;

            case 101:
                k1 = 1; k3 = 1; k5 = 1;
                break;

            case 011:
                k2 = 1; k3 = 1; k6 = 1;
                break;

            case 111:
                k1 = 1; k2 = 1; k3 = 1;
                k4 = 1; k5 = 1; k6 = 1;
                break;

            default:
                printf("\nERROR : Invalid process selection\n");
                return 0;
        }

        config_value.k1 = k1;
        config_value.k2 = k2;
        config_value.k3 = k3;
        config_value.k4 = k4;
        config_value.k5 = k5;
        config_value.k6 = k6;

        printf("Configuration valid. Simulation can start.\n");

        

        /*Renormalization*/
        Physics R(config_value);
        double alphaRenormQED;
        double alphaRenormS;
        double mf[9] = {config_value.m_e,config_value.m_mu,config_value.m_tau,config_value.m_u,config_value.m_d,config_value.m_c,config_value.m_s,config_value.m_top,config_value.m_b};
        double Qf[9] = {-1.0, -1.0, -1.0, 2.0/3.0, -1.0/3.0,  2.0/3.0,-1.0/3.0,  2.0/3.0, -1.0/3.0};

        char renorm_char;

        printf("Apply running of the coupling constants? (y/n): ");
        scanf(" %c",&renorm_char);
        if (renorm_char == 'y' || renorm_char == 'Y')
        {
            alphaRenormQED = R.renormalisation(0,config_value.renorm_scale,config_value.m_z,config_value.alphaQED,mf,Qf,config_value.Nc);
            alphaRenormS = R.renormalisation(1,config_value.renorm_scale,config_value.m_z,config_value.alphaS,mf,Qf,config_value.Nc);
            config_value.alphaQED = alphaRenormQED;
            config_value.alphaS = alphaRenormS;
            printf("alphaRenormQED = %.10e\n",alphaRenormQED);
            printf("alphaRenormS = %.10e\n",alphaRenormS);
        }
        else if(renorm_char == 'n' || renorm_char == 'N')
        {
            alphaRenormQED = config_value.alphaQED;
            alphaRenormS = config_value.alphaS;
            printf("No renormalization\n");
        }


        ///////////////////////////////////////
        /*generate all derived constants*/
        /////////////////////////////////////
        config_value.sin2_theta_W = 1 - (config_value.m_w/config_value.m_z)*(config_value.m_w/config_value.m_z);
        config_value.sinw = sqrt(config_value.sin2_theta_W);
        config_value.cosw = sqrt(1-config_value.sin2_theta_W);
        config_value.e = sqrt(4*M_PI*alphaRenormQED);

        config_value.gW = config_value.e/config_value.sinw;

        config_value.K_gamma = 2*config_value.e*config_value.e/3;
        config_value.K_z = (config_value.gW*config_value.gW)/(config_value.cosw*config_value.cosw);
        config_value.K_h = config_value.gW*config_value.gW/(4*config_value.m_w*config_value.m_w);
        config_value.cem = -0.5 + config_value.sin2_theta_W;
        config_value.cep = config_value.sin2_theta_W;
        config_value.ctm = 0.5 - 2*config_value.sin2_theta_W/3;
        config_value.ctp = -2*config_value.sin2_theta_W/3;
        config_value.be = sqrt(1 - 4*config_value.m_e*config_value.m_e/(config_value.sqrt_s*config_value.sqrt_s));
        config_value.bt = sqrt(1 - 4*config_value.m_top*config_value.m_top/(config_value.sqrt_s*config_value.sqrt_s));
        
        config_value.g_gamma = (config_value.Nc*config_value.K_gamma*config_value.K_gamma)/4;
        config_value.gh = (config_value.Nc*(config_value.K_h*config_value.K_h)*(config_value.sqrt_s*config_value.sqrt_s*config_value.sqrt_s*config_value.sqrt_s))/(64*(config_value.sqrt_s*config_value.sqrt_s-config_value.m_h*config_value.m_h)*(config_value.sqrt_s*config_value.sqrt_s-config_value.m_h*config_value.m_h) + (config_value.m_h*config_value.gamma_h)*(config_value.m_h*config_value.gamma_h));
        config_value.gz = (config_value.Nc*(config_value.K_z*config_value.K_z)*(config_value.sqrt_s*config_value.sqrt_s*config_value.sqrt_s*config_value.sqrt_s))/(64*(config_value.sqrt_s*config_value.sqrt_s-config_value.m_z*config_value.m_z)*(config_value.sqrt_s*config_value.sqrt_s-config_value.m_z*config_value.m_z) + (config_value.m_z*config_value.gamma_z)*(config_value.m_z*config_value.gamma_z));
        config_value.gzh = (config_value.Nc * config_value.K_z * config_value.K_h / 32) * (config_value.bt*config_value.be*(1 - config_value.bt*config_value.bt) * (1 - config_value.be*config_value.be))*( pow(config_value.sqrt_s,6)*( (config_value.sqrt_s*config_value.sqrt_s - config_value.m_z*config_value.m_z)*(config_value.sqrt_s*config_value.sqrt_s - config_value.m_h*config_value.m_h) + config_value.m_h*config_value.m_z*config_value.gamma_h*config_value.gamma_z ) )/ ( ( pow((config_value.sqrt_s*config_value.sqrt_s - config_value.m_h*config_value.m_h),2) + pow((config_value.m_h*config_value.gamma_h),2) )*( pow((config_value.sqrt_s*config_value.sqrt_s - config_value.m_z*config_value.m_z),2) + pow((config_value.m_z*config_value.gamma_z),2) ) );
        double s = config_value.sqrt_s*config_value.sqrt_s;
        config_value.gzgamma = (s*(s-config_value.m_z*config_value.m_z)*config_value.Nc*config_value.K_z * config_value.K_gamma / 8) / ( pow((config_value.sqrt_s*config_value.sqrt_s - config_value.m_z*config_value.m_z),2) + pow((config_value.m_z*config_value.gamma_z),2) );

        Physics P(config_value);



        /*numerical Rieman intergration error calculation of the cross section*/

        const double GEV2_TO_BARN = 3.89379e-4;

        double dtheta_1 = 1e-3;
        double dtheta_2 = dtheta_1/2; 
        double sigma_1 = 1e12 * GEV2_TO_BARN * P.total_x_section(dtheta_1);
        double sigma_2 = 1e12 * GEV2_TO_BARN * P.total_x_section(dtheta_2);

        double sigma = sigma_2;
        double error = std::abs(sigma_2-sigma_1);

        printf("Total cross section (Riemann integrator) %.10e +- %.10e pb\n",sigma,error);

        /*Individual contribution for each diagramm/interference term*/

        XSectionContrib sigma_decomp = P.total_x_section_contrib(dtheta_2);
        double sigma_tot = 1e12 * GEV2_TO_BARN *  (sigma_decomp.gamma + sigma_decomp.Z + sigma_decomp.H + sigma_decomp.gammaZ + sigma_decomp.gammaH + sigma_decomp.ZH);

        double contr_y = 100*1e12 * GEV2_TO_BARN *  (sigma_decomp.gamma)/sigma_tot;
        double contr_z = 100*1e12 * GEV2_TO_BARN *  (sigma_decomp.Z)/sigma_tot;
        double contr_h = 100 * 1e12 * GEV2_TO_BARN *  (sigma_decomp.H)/sigma_tot;
        double contr_yz = 100 * 1e12 * GEV2_TO_BARN *  (sigma_decomp.gammaZ)/sigma_tot;
        double contr_yh = 100 * 1e12 * GEV2_TO_BARN *  (sigma_decomp.gammaH)/sigma_tot;
        double contr_zh = 100 * 1e12 * GEV2_TO_BARN *  (sigma_decomp.ZH)/sigma_tot;

        printf("Process Decomposition :\n%.10e %% My\n%.10e %% Mz\n%.10e Mh\n%.10e %% Myz\n%.10e %% Myh\n%.10e %% Mzh\n",contr_y,contr_z,contr_h,contr_yz,contr_yh,contr_zh);
        
        
        LHE_Init i{};
        i.PDG_ID_1=-11, i.PDG_ID_2=11;
        i.beam_energy_1=config_value.sqrt_s/2, i.beam_energy_2=config_value.sqrt_s/2;
        i.Family_ID_1=0, i.Family_ID_2=0;
        i.Set_ID_1=247000, i.Set_ID_2=247000;
        i.weight_strategy=config_value.weight_strategy, i.number_physical_processes=1;
        i.Xsection_value=sigma, i.Xsection_error=error, i.maximum_weight_value=sigma, i.Process_ID=1;

        LHE_Event_Header h{};
        h.number_generated_particle=4, h.Process_ID=1;
        h.event_weight=1, h.renorm_scale=config_value.renorm_scale;
        h.alphaS=alphaRenormS, h.alphaQED=alphaRenormQED;

        LHE lhe(config_value,h,i);

        lhe.write_config(name);
        lhe.write_init(name);

        time_t t0 = time(nullptr);
        MonteCarlo(name,lhe,config_value,P);
        time_t t1 = time(nullptr);
        printf("Elapsed time : %ld s\n", t1 - t0);


        /* loop while(true) test */
        char again; 
        std::cout << "Generate another file (y/n) ? : ";
        std::cin >> again;
        std::cout << "" << std::endl;

        if (again != 'Y' && again != 'y') {
            break;
        }

    }



    return 0;
}