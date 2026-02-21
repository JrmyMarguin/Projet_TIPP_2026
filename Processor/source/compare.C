#include <TFile.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TSystem.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TMath.h>
#include <TLine.h>

#include <iostream> 
#include <string>
#include <fstream>


// X = MadGraph5 ; Y = Eco+


void compare(const char* X, const char* Y) {


    gStyle->SetTextFont(42);
    gStyle->SetOptStat(0); // Remove the stats box

    // .root file reading from filename
    std::string inputX = std::string("./Results/") + X + "/" + X + ".root";
    std::string inputY = std::string("./Results/") + Y + "/" + Y + ".root";
    const char* cinputX = inputX.c_str();
    const char* cinputY = inputY.c_str();

    TFile *fileX = new TFile(cinputX, "READ");
    TFile *fileY = new TFile(cinputY, "READ");

    // Check that files are properly opened
    if(!fileX->IsOpen() || !fileY->IsOpen()) {
         std::cout << "EERROR : Unable to open .root file\n" << std::endl;
         return;
    }
    

    // Comparison of histograms : cos(theta)

    TH1F *h_csX = (TH1F*) fileX->Get("histo_costheta");
    TH1F *h_csY = (TH1F*) fileY->Get("histo_costheta");
    h_csX->Sumw2(); h_csY->Sumw2();

    // Check statistics of the two histograms to compare
    if (h_csX->GetEntries() != h_csY->GetEntries()){
        std::cout << " Histogram from \"" << X << "\" has : " << h_csX->GetEntries() << " entries" << std::endl;
        std::cout << " While histogram from \"" << Y << "\" has : " << h_csY->GetEntries() << " entries\n" << std::endl;
        std::cout << "Please, select histograms with same statistics." << std::endl;
        return;
    }


    // Canvas creation and pad layout configuration
    TCanvas *c1 = new TCanvas("c1","c1",1000,800);
    
    // Pad position Parameter
    double x_left  = 0.0;
    double x_right = 1.0;
    double y_bottom = 0.0;
    double y_middle = 0.25;
    double y_top = 1.0;

    // Upper plot + relative error in the lower pad
    TPad *p1 = new TPad("p1","p1", x_left, y_middle, x_right, y_top);
    TPad *p2 = new TPad("p2","p2", x_left, y_bottom, x_right, y_middle);

    p1->SetTopMargin(0.1);
    p1->SetBottomMargin(0.0);
    p2->SetTopMargin(0.0);
    p2->SetBottomMargin(0.35);

    p1->Draw(); 
    p2->Draw(); 


    // Create X_Y director in ./Results_compare after comparison
    std::string output_compare = std::string("./Results_compare");
    std::string dir_name = output_compare + "/" + X + "_" + Y;
    gSystem->mkdir(dir_name.c_str());


    // Two cos(theta) histogramms on the same plot
    p1->cd();
    h_csX->SetLineColor(2);
    h_csY->SetLineColor(4);
    h_csX->SetTitle("Angular distribution of events in cos#theta");
    h_csX->SetMinimum(h_csX->GetMinimum()*0.7);
    h_csX->SetMaximum(h_csX->GetMaximum()*1.3);
    h_csX->GetXaxis()->SetLabelSize(0);
    h_csX->GetXaxis()->SetTitleSize(0);
    h_csX->Draw ("hist E1");
    h_csY->Draw ("hist E1 same");
    
    // Relative error between the two histograms (X-Y)/X
    p2->cd()->SetGrid();
    TH1F *h_rel = new TH1F("relative", "", 25, -1, 1);
    h_rel->Add(h_csX);
    h_rel->Add(h_csY, -1);
    h_rel->Divide(h_csX);

    h_rel->SetMarkerStyle(kMultiply);
    h_rel->SetMaximum(h_rel->GetMaximum()*2.0);
    h_rel->SetMinimum(abs(h_rel->GetMinimum())*(-2.0));
    h_rel->GetXaxis()->SetLabelSize(0.07);
    h_rel->GetYaxis()->SetTitle("Relative error");
    h_rel->GetYaxis()->SetNdivisions(505);
    h_rel->GetYaxis()->SetTitleSize(0.10);
    h_rel->GetYaxis()->SetLabelSize(0.08);
    h_rel->GetYaxis()->SetTitleOffset(0.4);
    h_rel->GetXaxis()->SetTitle("cos#theta");
    h_rel->GetXaxis()->SetTitleSize(0.12);
    h_rel->GetXaxis()->SetLabelSize(0.10);
    h_rel->GetXaxis()->SetTitleOffset(1.0);
    h_rel->Draw ("P");

    TLine *line_rel = new TLine(-1, 0, 1, 0);
    line_rel->SetLineStyle(1);
    line_rel->Draw();

    // Legend for cos(theta)
    auto *legend1 = new TLegend(0.4, 0.75, 0.60, 0.90); 
    p1->cd();
    legend1->AddEntry(h_csX, "MadGraph5", "l");
    legend1->AddEntry(h_csY, "Eco+", "l");
    legend1->SetTextSize(0.037);
    legend1->SetBorderSize(0);
    legend1->SetFillStyle(0);
    legend1->Draw ();

    std::string png_name1 = std::string("cos(theta)_") + X + "_" + Y;
    c1->SaveAs(Form("%s/%s.png", dir_name.c_str(), png_name1.c_str()));
    c1->Close();


    ///////////////////////////////////////
    // Comparison of helicity histogramms//
    ///////////////////////////////////////
    TCanvas *c2 = new TCanvas("c2","c2",1000,800);

    // Upper and Lower Pads
    TPad *p3 = new TPad("p3","p3", x_left, y_middle, x_right, y_top);
    TPad *p4 = new TPad("p4","p4", x_left, y_bottom, x_right, y_middle);

    p3->SetTopMargin(0.1);
    p3->SetBottomMargin(0.0);
    p4->SetTopMargin(0.0);
    p4->SetBottomMargin(0.35);

    p3->Draw();
    p4->Draw();


    TH1F *h_helX = (TH1F*) fileX->Get("histo_helicity");
    TH1F *h_helY = (TH1F*) fileY->Get("histo_helicity");

    h_helX->Sumw2(); h_helY->Sumw2();

    const char* helicity_label[16] = {
        "----", "---+", "--+-", "--++",
        "-+--", "-+-+", "-++-", "-+++",
        "+---", "+--+", "+-+-", "+-++",
        "++--", "++-+", "+++-", "++++"
    };

    // two helicity histogramms on the same plot
    p3->cd();
    h_helX->SetLineColor(2);
    h_helY->SetLineColor(4);
    h_helY->SetTitle("Helicities combinaisons of e-e+t#bar{t}");
    h_helY->SetMaximum(h_helX->GetMaximum()*1.3);
    h_helY->SetMinimum(1);
    h_helY->GetXaxis()->SetLabelSize(0);
    h_helY->GetXaxis()->SetTitleSize(0);
    h_helY->Draw ("hist E1");
    h_helX->Draw ("hist E1 same");
    
    // Relative Error between the two histogramms
    p4->cd()->SetGrid();
    TH1F *h_rel2 = new TH1F("relative2", "", 16, 0, 16);
    for (int i = 0; i < 16; ++i) {
        h_rel2->GetXaxis()->SetBinLabel(i+1, helicity_label[i]); 
    }  

    h_rel2->Add(h_helX);
    h_rel2->Add(h_helY, -1);
    h_rel2->Divide(h_helX);

    h_rel2->SetMaximum(h_rel2->GetMaximum()*3.3);
    h_rel2->SetMinimum(abs(h_rel2->GetMinimum())*(-3.));
    h_rel2->GetXaxis()->SetLabelSize(0.17);

    h_rel2->SetMarkerStyle(kMultiply);
    h_rel2->GetYaxis()->SetTitle("Relative error");
    h_rel2->GetYaxis()->SetNdivisions(505);
    h_rel2->GetYaxis()->SetTitleSize(0.10);
    h_rel2->GetYaxis()->SetLabelSize(0.08);
    h_rel2->GetYaxis()->SetTitleOffset(0.4);
    h_rel2->GetXaxis()->SetTitle("Helicities combinations");
    h_rel2->GetXaxis()->SetTitleSize(0.12);
    h_rel2->GetXaxis()->SetTitleOffset(1.0);
    h_rel2->Draw ("P");

    TLine *line_rel2 = new TLine(0, 0, 16, 0);
    line_rel2->SetLineStyle(1);
    line_rel2->Draw();

    // Legend for helicities
    auto *legend2 = new TLegend(0.7, 0.75, 0.9, 0.90);
    p3->cd();
    legend2->AddEntry(h_helX, "MadGraph5", "l");
    legend2->AddEntry(h_helY, "Eco+", "l");
    legend2->SetTextSize(0.037);
    legend2->SetBorderSize(0);
    legend2->SetFillStyle(0);
    legend2->Draw ();



    // Save result as .png
    std::string png_name2 = std::string("helicities_") + X + "_" + Y;
    c2->SaveAs(Form("%s/%s.png", dir_name.c_str(), png_name2.c_str()));



    // - Forward/Backward coefficient
    int bin0_X = h_csX->FindBin(0.0);
    double forward_X = h_csX->Integral(1, bin0_X);
    double backward_X = h_csX->Integral(bin0_X, h_csX->GetNbinsX());
    int bin0_Y = h_csY->FindBin(0.0);
    double forward_Y = h_csY->Integral(1, bin0_Y);
    double backward_Y = h_csY->Integral(bin0_Y, h_csY->GetNbinsY());

    double A_FB_X = (forward_X - backward_X)/(forward_X + backward_X);
    double A_FB_Y = (forward_Y - backward_Y)/(forward_Y + backward_Y);

    double sigma_A_X = sqrt((1 - pow(A_FB_X, 2))/(forward_X + backward_X));
    double sigma_A_Y = sqrt((1 - pow(A_FB_Y, 2))/(forward_Y + backward_Y));

    // - Global khi-square test (pvalue)
    double chi2 = h_csX->Chi2Test(h_csY, "WW CHI2");
    int dof = h_csX->GetNbinsX() - 1; // Dof Approximation
    double p_value = TMath::Prob(chi2, dof);



    // Bin by Bin statistical test
    double variance;
    double chi2_bin;
    TCanvas *c3 = new TCanvas("c3","c3",1000,800);
    TH1D* h_pvalue = new TH1D("h_pvalue", "", 25, -1, 1);
    h_pvalue->SetTitle("p-value bin-to-bin for cos#theta;cos#theta");
    h_pvalue->GetXaxis()->SetTitle("p-value");
    for (int i=0; i <= h_csX->GetNbinsX(); i++) {
        int NX = h_csX->GetBinContent(i); 
        int NY = h_csY->GetBinContent(i);
        variance = NX + NY;
        chi2_bin = pow(NX - NY, 2)/variance;

        h_pvalue->SetBinContent(i, TMath::Prob(chi2_bin, 1));
    }
    h_pvalue->SetMinimum(0);
    // h_pvalue->SetMaximum(1);
    h_pvalue->Draw("hist P*");

    TLine *p_crit = new TLine(-1, 0.05, 1, 0.05);
    p_crit->SetLineColor(kRed);
    p_crit->SetLineStyle(2);
    p_crit->Draw();


    // Save result as .png
    std::string png_name3 = std::string("p_value_") + X + "_" + Y;
    c3->SaveAs(Form("%s/%s.png", dir_name.c_str(), png_name3.c_str()));


    // Transverse momentum pT of top quark
    TH1F *h_ptX = (TH1F*) fileX->Get("histo_pt");
    TH1F *h_ptY = (TH1F*) fileY->Get("histo_pt");
    h_ptX->Sumw2(); h_ptY->Sumw2();


    // Canvas creation and pad layout configuration
    TCanvas *c4 = new TCanvas("c4","c4",1000,800);

    // Upper plot + relative error in the lower pad
    TPad *p5 = new TPad("p5","p5", x_left, y_middle, x_right, y_top);
    TPad *p6 = new TPad("p6","p6", x_left, y_bottom, x_right, y_middle);

    p5->SetTopMargin(0.1);
    p5->SetBottomMargin(0.0);
    p6->SetTopMargin(0.0);
    p6->SetBottomMargin(0.35);

    p5->Draw(); 
    p6->Draw(); 

    // Two pT histogramms on the same plot
    p5->cd();
    h_ptX->SetLineColor(2);
    h_ptY->SetLineColor(4);
    h_ptX->SetTitle("pT [GeV]");
    h_ptX->SetMinimum(1.);
    h_ptX->SetMaximum(h_ptX->GetMaximum()*1.3);
    int max_pt = h_ptX->GetXaxis()->GetBinUpEdge(h_ptX->FindLastBinAbove(0)) + 1;
    h_ptX->GetXaxis()->SetRangeUser(0, max_pt);
    h_ptY->GetXaxis()->SetRangeUser(0, max_pt);
    h_ptX->GetXaxis()->SetLabelSize(0);
    h_ptX->GetXaxis()->SetTitleSize(0);
    h_ptX->Draw ("hist E1");
    h_ptY->Draw ("hist E1 same");
    
    // Relative error between the two histograms (X-Y)/X
    p6->cd()->SetGrid();
    TH1F *h_rel3 = new TH1F("relative", "", 100, 0, 365);
    h_rel3->GetXaxis()->SetRangeUser(0, max_pt);
    h_rel3->Add(h_ptX);
    h_rel3->Add(h_ptY, -1);
    h_rel3->Divide(h_ptX);

    h_rel3->SetMarkerStyle(kMultiply);
    h_rel3->SetMaximum(h_rel3->GetMaximum()*3.5);
    h_rel3->SetMinimum(abs(h_rel3->GetMinimum())*(-3.5));
    h_rel3->GetXaxis()->SetLabelSize(0.07);
    h_rel3->GetYaxis()->SetTitle("Relative error");
    h_rel3->GetYaxis()->SetNdivisions(505);
    h_rel3->GetYaxis()->SetTitleSize(0.10);
    h_rel3->GetYaxis()->SetLabelSize(0.08);
    h_rel3->GetYaxis()->SetTitleOffset(0.4);
    h_rel3->GetXaxis()->SetTitle("pt [GeV]");
    h_rel3->GetXaxis()->SetTitleSize(0.12);
    h_rel3->GetXaxis()->SetLabelSize(0.10);
    h_rel3->GetXaxis()->SetTitleOffset(1.0);
    h_rel3->Draw ("P");

    TLine *line_rel3 = new TLine(0, 0, max_pt, 0);
    line_rel3->SetLineStyle(1);
    line_rel3->Draw();

    // Legend for pt
    auto *legend3 = new TLegend(0.4, 0.75, 0.60, 0.90); 
    p1->cd();
    legend3->AddEntry(h_csX, "MadGraph5", "l");
    legend3->AddEntry(h_csY, "Eco+", "l");
    legend3->SetTextSize(0.037);
    legend3->SetBorderSize(0);
    legend3->SetFillStyle(0);
    legend3->Draw ();

    std::string png_name4 = std::string("pT_") + X + "_" + Y;
    c4->SaveAs(Form("%s/%s.png", dir_name.c_str(), png_name4.c_str()));
    c4->Close();

    // - Mean of pT + uncertainty
    double mean_ptX = h_ptX->GetMean();
    double mean_ptY = h_ptY->GetMean();

    double sigma_mean_ptX = h_ptX->GetMeanError();
    double sigma_mean_ptY = h_ptY->GetMeanError();

    // - Invariant mass
    TH1F *h_mtt_X = (TH1F*) fileX->Get("histo_mtt");
    TH1F *h_mtt_Y = (TH1F*) fileY->Get("histo_mtt");

    double mtt_X = h_mtt_X->GetMean();
    double mtt_Y = h_mtt_Y->GetMean();

    double sigma_mtt_X = h_mtt_X->GetMeanError();
    double sigma_mtt_Y = h_mtt_Y->GetMeanError();

    // - Polarization of top quark
    double NtpY = h_helY->GetBinContent(7) + h_helY->GetBinContent(8) +
                 h_helY->GetBinContent(11) + h_helY->GetBinContent(12); 
    double NtmY = h_helY->GetBinContent(5) + h_helY->GetBinContent(6) +
                 h_helY->GetBinContent(9) + h_helY->GetBinContent(10);
    double polaY = abs(NtpY - NtmY)/(NtpY + NtmY);
    double sigma_polaY = sqrt((1 - pow(polaY,2))/(NtpY + NtmY));

    double NtpX = h_helX->GetBinContent(7) + h_helX->GetBinContent(8) +
                 h_helX->GetBinContent(11) + h_helX->GetBinContent(12); 
    double NtmX = h_helX->GetBinContent(5) + h_helX->GetBinContent(6) +
                 h_helX->GetBinContent(9) + h_helX->GetBinContent(10);
    double polaX = abs(NtpX - NtmX)/(NtpX + NtmX);
    double sigma_polaX = sqrt((1 - pow(polaX,2))/(NtpX + NtmX));



    // Output .txt file after comparison 
    std::string txt_name = dir_name + "/" + X + "_" + Y + ".txt";
    std::ofstream fichier(txt_name);


    if (fichier.is_open()) {
        fichier << "Asymmetry coefficients" << std::endl;
        fichier << " Forward-Backward asymmetry of MadGraph5 : A_FB = " << A_FB_X 
                << "  \xC2\xB1 " << sigma_A_X << std::endl;
        fichier << " Forward-Backward asymmetry of Eco+ : A_FB = " << A_FB_Y 
                << "  \xC2\xB1 " << sigma_A_Y << "\n"<< std::endl;

        fichier << "Global p-value" << std::endl;
        fichier << " Chi2 = " << chi2  << std::endl;
        fichier << " Degrees of freedom = " << dof  << std::endl;
        fichier << " p-value = " << p_value << "\n" << std::endl;

        fichier << "Transverse momentum pT" << std::endl;
        fichier << " MadGraph5 <pT> = " << mean_ptX << " \xC2\xB1 " << sigma_mean_ptX << " GeV" << std::endl;
        fichier << " Eco+ <pT> = " << mean_ptY << " \xC2\xB1 " << sigma_mean_ptY << " GeV" << "\n" << std::endl;

        fichier << "Invariant mass" << std::endl;
        fichier << " MadGraph5 mtt = " << mtt_X << " \xC2\xB1 " << sigma_mtt_X << " GeV" << std::endl;
        fichier << " Eco+ mtt = " << mtt_Y << " \xC2\xB1 " << sigma_mtt_Y << " GeV" << "\n" << std::endl;

        fichier << "Top quark polarization" << std::endl;
        fichier << " MadGraph5 P = " << polaX << " \xC2\xB1 " << sigma_polaX << std::endl;
        fichier << " Eco+ P = " << polaY << " \xC2\xB1 " << sigma_polaY << std::endl;


        fichier.close();  // Close the file
        // std::cout << "OK" << std::endl;
    } else {
        std::cout << "ERROR : Unable to write in .txt file output" << std::endl;
    }


    return;
}