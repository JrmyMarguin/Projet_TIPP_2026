#include <iostream>
#include <fstream>
#include <sstream>
#include <vector> 
#include <string>

#include <TH1D.h>
#include <TFile.h>
#include <TSystem.h>

#include "Structure.h"

/* install root  : 
$ sudo snap install root-framework
$ snap run root-framework */


int main() {

    std::cout << "***************************************************************************************************\n" << std::endl;
    std::cout << "*                                           e+ e- -> t t̄                                          *\n" << std::endl; 
    std::cout << "*                                           LHE Analyzer                                          *\n" << std::endl; 
    std::cout << "***************************************************************************************************\n" << std::endl; 

    while (true) {
        /* Reading directory filename */
        
        std::string name; 
        std::cout << "Directory name : ";
        std::cin >> name;
        
        std::string filename = std::string("./Results/") +  name + "/" +  name + "_LHE.txt"; /* Select _LHE.txt  file */

        std::ifstream file(filename); /* open file for reading */

        if (!file) {
            std::cout << "ERROR : unable to read _LHE.txt file \n" << std::endl;
            return 0;
        }

        /* Reading and store _LHE.txt (case : where all events corresponds to the same process i.e same number of inital particles and final particles */

        std::string line; /* Initalization of line reading object */
        int N_part=0; /* Number of Particles in event header */ 
        int in_event = false; /* Event reading flag*/
        int read_header = false; /* Event HEADER line flag*/

        std::vector<LHE_Particle> electron; /* storing data in vectors */
        std::vector<LHE_Particle> positron;
        std::vector<LHE_Particle> top;
        std::vector<LHE_Particle> anti_top;
        std::vector<LHE_Particle> other; /* if intermediate particles are added by MADGRAPH */
        

        while (std::getline(file,line)) {

            if(line.find("<event>") != std::string::npos){ /* enter <event> tag */
                in_event = true;
                read_header = true; /* not necessary but good for code readability */
                continue; /* read next line*/
            }

            if(line.find("</event>") != std::string::npos){ /* exiting </event> tag */
                in_event = false;
                N_part = 0;
                continue;
            }

            if (in_event) { /* inside <event> tag */

                if (read_header) {

                    read_header = false; 
                    std::stringstream sh(line); /* Read number of particles from event header*/
                    sh >> N_part;

                    for (int i=0; i<N_part;i++) { /* For each Particle */

                        std::getline(file,line); /* Read particle line */
                        std::stringstream sp(line);

                        LHE_Particle p; /* store read particle data into structure*/

                        sp >> p.beam_id >> p.part_status >> p.mother_1 >> p.mother_2 >> p.color_flow_1 >> p.color_flow_2 >> p.px >> p.py >> p.pz >> p.E >> p.mass >> p.mean_distance >> p.helicity;

                        if (p.beam_id == 11) {electron.push_back(p);} 
                        else if (p.beam_id == -11) {positron.push_back(p);}
                        else if (p.beam_id == 6) {top.push_back(p);}
                        else if (p.beam_id == -6) {anti_top.push_back(p);}
                        else {other.push_back(p);}
                    }

                    continue;
                }

            }

        }

        file.close();
        

    
        /* Writing output.root file */

        std::string output_name = std::string("./Results/") + name + "/" + name + ".root";

        TFile *output = new TFile(output_name.c_str(),"RECREATE"); /* RECREATE : Overwwrite existing .root file */

        TH1D *histo_pt = new TH1D("histo_pt",(name + ";pT (GeV);Events").c_str(),100,0,2*electron[0].E); /* Last three values : Number of Bins, Min Value, Max Value (in X-axis) */
        TH1D *histo_eta = new TH1D("histo_eta",(name + ";#eta;Events").c_str(),100,-3,3);
        TH1D *histo_costheta = new TH1D("histo_costheta",(name + ";cos#theta;Events").c_str(),25,-1,1);
        TH1D *histo_mtt = new TH1D("histo_mtt",(name + ";M (GeV);Events").c_str(),100,0,4*electron[0].E);
        TH1D *histo_helicity = new TH1D("histo_helicity",(name + ";helicity(e#bar{e}t#bar{t});Events").c_str(),16,0,16);

        /* Create bins for the 16 helicity combinations*/

        const char* helicity_label[16] = {
        "----", "---+", "--+-", "--++",
        "-+--", "-+-+", "-++-", "-+++",
        "+---", "+--+", "+-+-", "+-++",
        "++--", "++-+", "+++-", "++++"
        };

        for (int i = 0; i < 16; ++i) {
            histo_helicity->GetXaxis()->SetBinLabel(i+1, helicity_label[i]); 
        }   

    /* Filling histograms */

        for(const LHE_Particle& t : top){
            histo_pt->Fill(t.pt());
            histo_eta->Fill(t.eta());
            histo_costheta->Fill(t.costheta());
        }

        int N_top = top.size(); /* in this study case N_top = N_anti_top = N_electron = N_positron */

        for (int i =0; i<N_top; i++){
            
            histo_mtt->Fill(LHE_Particle::inv_mass(top[i],anti_top[i]));

            int he = electron[i].helicity;
            int he_bar = positron[i].helicity;
            int ht = top[i].helicity;
            int ht_bar = anti_top[i].helicity; 

            int index = (he == 1)*8 + (he_bar == 1)*4 + (ht == 1)*2 + (ht_bar == 1)*1 ; /* binary encoding index to avoid 16 if-else statement */ 
            histo_helicity->Fill(index); 

        }

        output->Write();
        output->Close();
        delete output;

        /* Call ROOT macro for plotting */
        std::string cmd = "root -l -b -q './Processor/source/plot.C(\"";
        cmd = cmd + output_name;
        cmd = cmd + "\")'";

        gSystem->Exec(cmd.c_str());

        /* Test of while(true) loop */
        char again; 
        std::cout << " Analyze another file (y/n) ? : ";
        std::cin >> again;

        if (again != 'Y' && again != 'y') {
            break;
        }

    }

    return 0;
}   