#include <TFile.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TSystem.h>
#include <TStyle.h>


#include <iostream> 
#include <string>



void plot(const char* filename) {

    gStyle->SetOptStat(0); // Remove statistics box

    /* .root file reading */
    TFile *file = new TFile(filename,"READ"); 

    if(!file->IsOpen()) {
         std::cout << "ERROR : unable to open .root file\n" << std::endl;
         return;
    }

    TH1D *h_pt= (TH1D*) file->Get("histo_pt");
    TH1D *h_eta = (TH1D*) file->Get("histo_eta");
    TH1D *h_costheta = (TH1D*) file->Get("histo_costheta");
    TH1D *h_mtt = (TH1D*) file->Get("histo_mtt");
    TH1D *h_helicity = (TH1D*) file->Get("histo_helicity");

    if (!h_pt || !h_eta || !h_costheta || !h_mtt || !h_helicity){
        std::cout << "ERROR : unable to read all histograms\n" << std::endl;
        return;
    }

    /* Output directory creation*/
    std::string output_name = filename; /* Remove .root from the name */
    size_t pos = output_name.rfind(".root"); 

    if (pos != std::string::npos) {
        output_name = output_name.substr(0,pos);
    }

    output_name = output_name + "_output"; /* creation of a sub-directory <filename>_output  */

    gSystem->mkdir(output_name.c_str()); 


    /* Plotting / Save histogramms as .png */

    TCanvas *c = new TCanvas("c","plot",1500,1000);
    // X-axis adaptation 
    h_pt->GetXaxis()->SetRangeUser(0, (h_pt->GetXaxis()->GetBinUpEdge(h_pt->FindLastBinAbove(0)) + 1));
    h_pt->Draw("HIST");
    c->SaveAs(Form("%s/pt_top.png",output_name.c_str()));

    h_eta->Draw("HIST");
    c->SaveAs(Form("%s/eta_top.png",output_name.c_str()));

    h_costheta->Draw("HIST");
    h_costheta->SetMinimum(0.);
    h_costheta->SetMaximum(h_costheta->GetMaximum()*1.3);
    c->SaveAs(Form("%s/costheta_top.png",output_name.c_str()));

    h_mtt->Draw("HIST");
    c->SaveAs(Form("%s/mtt.png",output_name.c_str()));
    
    h_helicity->Draw("BAR");
    c->SaveAs(Form("%s/helicity.png",output_name.c_str()));
    
    file->Close();
    c->Close();

    delete file;
    delete c;

    return;
}