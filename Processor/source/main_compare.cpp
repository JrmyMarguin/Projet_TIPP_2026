#include <iostream>
#include <fstream>
#include <sstream>
#include <vector> 
#include <string>

#include <TH1D.h>
#include <TFile.h>
#include <TSystem.h>

#include "Structure.h"


int main() {

    std::cout << "***************************************************************************************************\n" << std::endl;
    std::cout << "*                                           e+ e- -> t t̄                                          *\n" << std::endl; 
    std::cout << "*                                           LHE Comparator                                        *\n" << std::endl; 
    std::cout << "***************************************************************************************************\n" << std::endl; 


    while(true) {
        // Read directory names and create filename parameters X and Y
        // passed to tne compare.C Function
        std::string X; 
        std::string Y; 
        std::cout << "MADGRAPH Directory: ";
        std::cin >> X;
        std::cout << "ECO+ Directory : ";
        std::cin >> Y;
    
        std::string filenameX = std::string("./Results/") +  X + "/" +  X + ".root"; // Select the .root files
        std::string filenameY = std::string("./Results/") +  Y + "/" +  Y + ".root"; // 

        // Check if file are readable
        std::ifstream fileX(filenameX); // open file for reading 
        std::ifstream fileY(filenameY);

        if (!fileX) {
            std::cout << "ERROR : unable to read MADGRAPH.root file \n" << std::endl;
            return 0;
        }
        if (!fileY) {
            std::cout << "ERROR : unable to read ECO+.root file \n" << std::endl;
            return 0;
        }


        // Call ROOT macro for plotting 
        std::string cmd2 = "root -l -b -q './Processor/source/compare.C(\"";
        cmd2 = cmd2 + X + "\", \"" + Y + "\")'";
        gSystem->Exec(cmd2.c_str());


        /* Test of while(true) loop*/
        char again; 
        std::cout << "Analyze another file (y/n) ? : ";
        std::cin >> again;

        if (again != 'Y' && again != 'y') {
            break;
        }

    }

    return 0;
    
}