#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <vector>
#include <utility>
#include <algorithm>
#include <stdint.h>
#include "TH1.h"
#include "TCanvas.h"
#include "TPad.h"
#include "stdlib.h" 
#include "TStyle.h"
#include "TH1D.h"
#include "TGraph.h"
#include "TAxis.h"
#include "TLegend.h"
#include "TFile.h"
#include "TDirectory.h"
#include "string.h"
#include "TChain.h"
#include "TTree.h"
#include "TROOT.h"
#include "TGaxis.h"
#include "TLatex.h"
#include "TBufferFile.h"
#include "TLorentzVector.h"
#include "TPaveStats.h"
#include <iomanip>      // std::setprecision
#include <sstream>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>

using namespace std;

void cms_data() {
    gROOT->Reset();
    gStyle->SetOptStat("nemruo");

    std::string inputFile = "Run2011A_DoubleElectron_merged.root";
    std::string outputDir = "output";
    std::string nameOutFile = "cmsdata_Run2011A_DoubleElectron";

    TChain* inputChain = new TChain("Events");
    inputChain->Add((inputFile).c_str());

    UInt_t nElectron;
    float Electron_pt[100];
    float Electron_eta[100];
    float Electron_deltaEtaSC[100];
    float Electron_phi[100];
    float Electron_mass[100];
    Int_t Electron_charge[100];
    float Electron_dxy[100];
    float Electron_dz[100];
    float Electron_sip3d[100];
    float Electron_pfRelIso03_all[100];
    bool Electron_isPFcand[100];
    UChar_t Electron_lostHits[100];
    bool Electron_isEE[100];
    bool Electron_isEB[100];

    vector< pair<int, float> > goodElectrons;

    TLorentzVector ele1, ele2;

    // Other
    double S;
    double S1;
    double S2;
    double SQM;

    double mZ = 91.1876;

    double s1, s2, s3, s4, s;
    double dx, dy, dz, rap, pz;

    UInt_t run;
    ULong64_t event;

      // No. of good electron
    TH1F *h_nge = new TH1F("NGoodElectron", "No. of good electron", 10, 0., 10.);
    h_nge->GetXaxis()->SetTitle("Number of Electrons");
    h_nge->GetYaxis()->SetTitle("Number of Events");


    // Z to 2e --------------------------------------------------------------------

    // ZTo2e mass spectrum with Gsf muon
    TH1F *h_mZ_2e = new TH1F("massZto2e", "mass of Z to 2e", 120, 40., 120.);
    h_mZ_2e->GetXaxis()->SetTitle("Invariant Mass for Nelectron=2 (in GeV/c^2)");
    h_mZ_2e->GetYaxis()->SetTitle("Number of Events");

    // Electron momentum
    TH1F *h_p_e = new TH1F("e_momentum", "Electron momentum", 200, 0., 200.);
    h_p_e->GetXaxis()->SetTitle("Momentum (GeV/c)");
    h_p_e->GetYaxis()->SetTitle("Number of Events");

    // Electron eT
    TH1F *h_et_e = new TH1F("e_eT", "Electron eT", 200, 0., 200.);
    h_et_e->GetXaxis()->SetTitle("eT (GeV/c)");
    h_et_e->GetYaxis()->SetTitle("Number of Events");

    // Electron pT before cuts
    TH1F *h_pt_e_b4 = new TH1F("b4_e_pT", "Electron pT", 200, 0, 200.);
    h_pt_e_b4->GetXaxis()->SetTitle("pT (GeV/c)");
    h_pt_e_b4->GetYaxis()->SetTitle("Number of Events");

    // Electron eta before cuts
    TH1F *h_eta_e_b4 = new TH1F("b4_e_eta", "Electron eta", 140, -3.5, 3.5);
    h_eta_e_b4->GetXaxis()->SetTitle("eta");
    h_eta_e_b4->GetYaxis()->SetTitle("Number of Events");

    // Electron phi
    TH1F *h_phi_e = new TH1F("e_phi", "Electron phi", 314, -3.17, 3.17);
    h_phi_e->GetXaxis()->SetTitle("Phi");
    h_phi_e->GetYaxis()->SetTitle("Number of Events");

    // Electron SuperCluster (SC) eta
    TH1F *h_sc_eta = new TH1F("e_SC_eta", "Electron SC eta", 140, -3.5, 3.5);
    h_sc_eta->GetXaxis()->SetTitle("Super Cluster Eta");
    h_sc_eta->GetYaxis()->SetTitle("Number of Events");

    // Relative isolation electron
    TH1F *h_relPFIso_e = new TH1F("e_RelPFIso", "R.PFIso", 100, 0., 5.);
    h_relPFIso_e->GetXaxis()->SetTitle("Relative Isolation e");
    h_relPFIso_e->GetYaxis()->SetTitle("Number of Events");

    // Relative isolation electron after cuts
    TH1F *h_relPFIso_e_after = new TH1F("after_e_RelPFIso", "R.PFIso", 100, 0., 5.);
    h_relPFIso_e_after->GetXaxis()->SetTitle("Relative Isolation e");
    h_relPFIso_e_after->GetYaxis()->SetTitle("Number of Events");

    // Transverse impact parameter with respect to primary vertex for Electron
    TH1F *h_dxy_e = new TH1F("e_dxy", "Electron dxy", 100, 0., 1.);
    h_dxy_e->GetXaxis()->SetTitle("Transverse impact parameter w.r.t. primvtx");
    h_dxy_e->GetYaxis()->SetTitle("Number of Events");

    // Electron pT after cuts for Zto2e
    TH1F *h_pt_e_after_Zto2e = new TH1F("after_e_pT_Zto2e", "e pT", 240, 0., 120.);
    h_pt_e_after_Zto2e->GetXaxis()->SetTitle("Electron pT (GeV/c)");
    h_pt_e_after_Zto2e->GetYaxis()->SetTitle("Number of Events");

    // Electron eta after cuts for Zto2e
    TH1F *h_eta_e_after_Zto2e = new TH1F("after_e_eta_Zto2e","e eta", 140, -3.5, 3.5);
    h_eta_e_after_Zto2e->GetXaxis()->SetTitle("Eta");
    h_eta_e_after_Zto2e->GetYaxis()->SetTitle("Number of Events");

    // Electron pT after cuts
    TH1F *h_pt_e_after = new TH1F("after_e_pT", "Electron pT", 240, 0., 120.);
    h_pt_e_after->GetXaxis()->SetTitle("Electron pT (GeV/c)");
    h_pt_e_after->GetYaxis()->SetTitle("Number of Events");

    // Electron eta after cuts
    TH1F *h_eta_e_after = new TH1F("after_e_eta", "Electron eta", 140, -3.5, 3.5);
    h_eta_e_after->GetXaxis()->SetTitle("eta");
    h_eta_e_after->GetYaxis()->SetTitle("Number of Events");

    // 2mu2e
    // Electron pT after cuts for 2mu2e
    TH1F *h_pt_e_after_2mu2e = new TH1F("after_e_pT_2mu2e", "e pT", 240, 0., 120.);
    h_pt_e_after_2mu2e->GetXaxis()->SetTitle("Electron pT (GeV/c)");
    h_pt_e_after_2mu2e->GetYaxis()->SetTitle("Number of Events");

    // Electron eta after cuts for 2mu2e
    TH1F *h_eta_e_after_2mu2e = new TH1F("after_e_eta_2mu2e","e eta", 140, -3.5, 3.5);
    h_eta_e_after_2mu2e->GetXaxis()->SetTitle("eta");
    h_eta_e_after_2mu2e->GetYaxis()->SetTitle("Number of Events");

    // SIP for muon
    TH1F *h_SIP3d_mu_b4 = new TH1F("SIP3d_mu", "SIP_3D for Muon", 100, 0., 10.);
    h_SIP3d_mu_b4->GetXaxis()->SetTitle("SIP_3D");
    h_SIP3d_mu_b4->GetYaxis()->SetTitle("Number of Events");

    // SIP for electron
    TH1F *h_SIP3d_e_b4 = new TH1F("SIP3d_e", "SIP_3D for Electron", 100, 0., 10.);
    h_SIP3d_e_b4->GetXaxis()->SetTitle("SIP_3D");
    h_SIP3d_e_b4->GetYaxis()->SetTitle("Number of Events");

    TH1F *h_misshite = new TH1F("e_misshit", "e track missing hits", 5, 0., 5.);
    h_misshite->GetXaxis()->SetTitle("gsfTrack Hit type");
    h_misshite->GetYaxis()->SetTitle("Number of Events");


    ////////////////////////
    // Data DoubleMu 2011 //
    ////////////////////////

    // Set all the status tree to 0  
    inputChain->SetBranchStatus("*", 0);
    inputChain->SetBranchStatus("nElectron", 1);  
    inputChain->SetBranchStatus("Electron_pt", 1); 
    inputChain->SetBranchStatus("Electron_eta", 1);
    inputChain->SetBranchStatus("Electron_deltaEtaSC", 1);
    inputChain->SetBranchStatus("Electron_phi", 1);
    inputChain->SetBranchStatus("Electron_mass", 1);
    inputChain->SetBranchStatus("Electron_charge", 1);
    inputChain->SetBranchStatus("Electron_dxy", 1);
    inputChain->SetBranchStatus("Electron_dz", 1);
    inputChain->SetBranchStatus("Electron_sip3d", 1);
    inputChain->SetBranchStatus("Electron_pfRelIso03_all", 1);
    inputChain->SetBranchStatus("Electron_isPFcand", 1);
    inputChain->SetBranchStatus("Electron_lostHits", 1);
    inputChain->SetBranchStatus("Electron_isEE", 1);
    inputChain->SetBranchStatus("Electron_isEB", 1);

    inputChain->SetBranchAddress("nElectron", &nElectron);  
    inputChain->SetBranchAddress("Electron_pt", Electron_pt); 
    inputChain->SetBranchAddress("Electron_eta", Electron_eta);
    inputChain->SetBranchAddress("Electron_deltaEtaSC", Electron_deltaEtaSC);
    inputChain->SetBranchAddress("Electron_phi", Electron_phi);
    inputChain->SetBranchAddress("Electron_mass", Electron_mass);
    inputChain->SetBranchAddress("Electron_charge", Electron_charge);
    inputChain->SetBranchAddress("Electron_dxy", Electron_dxy);
    inputChain->SetBranchAddress("Electron_dz", Electron_dz);
    inputChain->SetBranchAddress("Electron_sip3d", Electron_sip3d);
    inputChain->SetBranchAddress("Electron_pfRelIso03_all", Electron_pfRelIso03_all);
    inputChain->SetBranchAddress("Electron_isPFcand", Electron_isPFcand);
    inputChain->SetBranchAddress("Electron_lostHits", Electron_lostHits);
    inputChain->SetBranchAddress("Electron_isEE", Electron_isEE);
    inputChain->SetBranchAddress("Electron_isEB", Electron_isEB);

    //activating branches for debugging
    inputChain->SetBranchStatus("run", 1);
    inputChain->SetBranchStatus("event", 1);
    inputChain->SetBranchAddress("run", &run);  
    inputChain->SetBranchAddress("event", &event);  
  
    // ********************************* //  


    Int_t nevent = inputChain->GetEntries();

    for (int aa = 0; aa < 1000000; aa++) {
        std::cout << "Start" << std::endl;
    //  for (int aa=12000000;aa<13500000;aa++) {
    // if (aa > 100) {break;}
        if (aa % 100000 == 0) {cout << "loop " << aa << " / " << nevent << endl;}

        inputChain->GetEntry(aa);

        // Electrons ------------------------------------------------------------------------

        int nGoodElectrons = 0;
        goodElectrons.clear();

    // Loop over all electrons/event
        for (int ee = 0; ee < nElectron; ee++) {

            TLorentzVector ele;
            ele.SetPtEtaPhiM(Electron_pt[ee], Electron_eta[ee], Electron_phi[ee], Electron_mass[ee]);

            // Particle flow electrons
            if (Electron_isPFcand[ee]) {

                h_p_e->Fill(ele.P());
                h_et_e->Fill(ele.Et());
                h_pt_e_b4->Fill(Electron_pt[ee]);
                h_eta_e_b4->Fill(Electron_eta[ee]);
                h_sc_eta->Fill(Electron_deltaEtaSC[ee] + Electron_eta[ee]);
                h_phi_e->Fill(Electron_phi[ee]);
                h_misshite->Fill(Electron_lostHits[ee]);
                h_relPFIso_e->Fill(Electron_pfRelIso03_all[ee]);
                h_dxy_e->Fill(Electron_dxy[ee]);
                h_SIP3d_e_b4->Fill(Electron_sip3d[ee]);
    
                // Electron selection
                if (Electron_pt[ee] > 7. && abs(Electron_deltaEtaSC[ee] + Electron_eta[ee]) < 2.5) {

                    if (Electron_lostHits[ee] <= 1 && abs(Electron_sip3d[ee]) < 4.) {

                        if (abs(Electron_dxy[ee]) < 0.5 && abs(Electron_dz[ee]) < 1.) {

                            if (Electron_isEB[ee]) {

                                if (Electron_pfRelIso03_all[ee] < 0.4) {

                                    goodElectrons.push_back(make_pair(ee, Electron_pt[ee]));
                                    nGoodElectrons++;

                                }       
                            }

                            else if (Electron_isEE[ee]) {

                                if (Electron_pfRelIso03_all[ee] < 0.4) {

                                    goodElectrons.push_back(make_pair(ee, Electron_pt[ee]));
                                    nGoodElectrons++;

                                }
                            }  // EB/EE
                        }  // dxy, dz
                    }  // lost hits, SIP 3D
                }  // pT, SC eta
            }  // Particle flow electrons
        } // Loop over all electrons/event


        h_nge->Fill(nGoodElectrons);

        if (nGoodElectrons >= 2) {

            int e1 = goodElectrons.at(0).first;
            int e2 = goodElectrons.at(1).first;
     
            if (Electron_charge[e1] + Electron_charge[e2] == 0) {

                ele1.SetPtEtaPhiM(Electron_pt[e1], Electron_eta[e1], Electron_phi[e1], Electron_mass[e1]);
                ele2.SetPtEtaPhiM(Electron_pt[e2], Electron_eta[e2], Electron_phi[e2], Electron_mass[e2]);

                SQM = (0.5109989461e-3) * (0.5109989461e-3);
                S1 = sqrt( (ele1.P()*ele1.P() + SQM) * (ele2.P()*ele2.P() + SQM) );
                S2 = ele1.Px() * ele2.Px() + ele1.Py() * ele2.Py() + ele1.Pz() * ele2.Pz();
                S = sqrt(2. * (SQM + S1 - S2));

                h_mZ_2e->Fill(S);

                // pT and eta after the cuts
                for (int ee = 0; ee < nGoodElectrons; ee++) {

                    h_pt_e_after_Zto2e->Fill(Electron_pt[goodElectrons.at(ee).first]);
                    h_eta_e_after_Zto2e->Fill(Electron_eta[goodElectrons.at(ee).first]);

                }

            } // charge 1 + charge 2 = 0

        } // nGoodElectrons >= 2

    }

    TFile outputFile((nameOutFile + ".root").c_str(), "RECREATE");
    // TDirectory* dir = outputFile.mkdir("Events");
    // dir->cd();

    // Electron selection
    h_p_e->Write();
    h_et_e->Write();
    h_pt_e_b4->Write();
    h_eta_e_b4->Write();
    h_sc_eta->Write();
    h_phi_e->Write();
    h_misshite->Write();
    h_relPFIso_e->Write();
    h_dxy_e->Write();
    h_SIP3d_e_b4->Write();
    h_nge->Write();

    // Z to 2e
    h_mZ_2e->Write();
    h_pt_e_after_Zto2e->Write();
    h_eta_e_after_Zto2e->Write();

    outputFile.Close();
    gROOT->ProcessLine(".q");

}


int main() {

    cms_data();
    return 0;

}