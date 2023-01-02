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
#include <TApplication.h>
#include <TROOT.h>


 

using namespace std;

void cms_data() {
    gROOT->Reset();
    // gStyle->SetOptStat("nemruo");
    gStyle->SetOptStat("nemrous");

    std::string inputFile = "Run2011A_DoubleElectron_merged.root";
    std::string outputDir = "output";
    std::string nameOutFile = "cmsdata_Run2011A_DoubleElectron";

    TChain* inputChain = new TChain("Events");
    inputChain->Add((inputFile).c_str());

    UInt_t nElectron;
    // float Electron_pt[100];
    // float Electron_eta[100];
    // float Electron_deltaEtaSC[100];
    // float Electron_phi[100];
    // float Electron_mass[100];
    // Int_t Electron_charge[100];
    // float Electron_dxy[100];
    // float Electron_dz[100];
    // float Electron_sip3d[100];
    // float Electron_pfRelIso03_all[100];
    // bool Electron_isPFcand[100];
    // UChar_t Electron_lostHits[100];
    // bool Electron_isEE[100];
    // bool Electron_isEB[100];


	UInt_t nElec = 100, Electron_nNano[nElec];
	
	Int_t Electron_charge[nElec], Electron_cutBased[nElec], Electron_simId[nElec], Electron_tightCharge[nElec], Electron_vtxIdx[nElec], Electron_genPartIdx[nElec];
	UChar_t Electron_lostHits[nElec];
	Bool_t Electron_convVeto[nElec], Electron_convVetoOld[nElec], Electron_isEB[nElec], Electron_isEE[nElec], Electron_isNano[nElec], Electron_isPFcand[nElec];

    Float_t Electron_SCeta[nElec], Electron_convDcot[nElec], Electron_convDist[nElec],
        Electron_deltaEtaSC[nElec], Electron_deltaEtaSCtr[nElec],
        Electron_deltaPhiSC[nElec], Electron_deltaPhiSCtr[nElec], Electron_dr03EcalRecHitSumEt[nElec],
        Electron_dr03EcalRecHitSumEtOld[nElec], Electron_dr03HcalDepth1TowerSumEt[nElec], Electron_dr03HcalDepth1TowerSumEtOld[nElec],
        Electron_dr03HcalTowerSumEt[nElec], Electron_dr03TkSumPt[nElec], Electron_dr03TkSumPtOld[nElec], Electron_dxy[nElec],
        Electron_dxyErr[nElec], Electron_dz[nElec], Electron_dzErr[nElec], Electron_eInvMinusPInv[nElec],
        Electron_eInvMinusPInvOld[nElec], Electron_eta[nElec], Electron_hoe[nElec],
        Electron_ip3d[nElec], Electron_mass[nElec], Electron_pfRelIso03_all[nElec], Electron_pfRelIso03_chg[nElec],
        Electron_phi[nElec], Electron_pt[nElec], Electron_sieie[nElec],
        Electron_sieieR1[nElec], Electron_sip3d[nElec], Electron_x[nElec], Electron_y[nElec], Electron_z[nElec];

    Float_t PV_x, PV_y, PV_z, PV_chi2;

    vector< pair<int, float> > goodElectrons;

    TLorentzVector ele1, ele2;

    // Other
    double S;
    double S1;
    double S2;
    double SQM = (0.5109989461e-3) * (0.5109989461e-3);;

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
    TH1F* h_mZ_2e = new TH1F("massZto2e", "mass of Z to 2e", 120, 40., 120.);
    h_mZ_2e->GetXaxis()->SetTitle("Invariant Mass for Nelectron=2 (in GeV/c^2)");
    h_mZ_2e->GetYaxis()->SetTitle("Number of Events");


    TH1F* h_Electron_SCeta_after_Zto2e = new TH1F("Electron_SCeta", "Electron_SCeta", 100u, -3.14, 3.14);
    TH1F* h_Electron_charge_after_Zto2e = new TH1F("Electron_charge", "Electron_charge", 100u, -2, 2);
    TH1F* h_Electron_convDcot_after_Zto2e = new TH1F("Electron_convDcot", "Electron_convDcot", 100u, -2.0, 2.0);
    TH1F* h_Electron_convDist_after_Zto2e = new TH1F("Electron_convDist", "Electron_convDist", 100u, -6., 2.);
    TH1F* h_Electron_convVeto_after_Zto2e = new TH1F("Electron_convVeto", "Electron_convVeto", 100u, 1., 1.0006);
    TH1F* h_Electron_convVetoOld_after_Zto2e = new TH1F("Electron_convVetoOld", "Electron_convVetoOld", 100u, 1., 1.0006);
    TH1F* h_Electron_cutBased_after_Zto2e = new TH1F("Electron_cutBased", "Electron_cutBased", 100u, 0., 5.);
    TH1F* h_Electron_deltaEtaSC_after_Zto2e = new TH1F("Electron_deltaEtaSC", "Electron_deltaEtaSC", 100u, -0.10, 0.10);
    TH1F* h_Electron_deltaEtaSCtr_after_Zto2e = new TH1F("Electron_deltaEtaSCtr", "Electron_deltaEtaSCtr", 100u, -0.01, 0.01);
    TH1F* h_Electron_deltaPhiSC_after_Zto2e = new TH1F("Electron_deltaPhiSC", "Electron_deltaPhiSC", 100u, -0.06, 0.06);
    TH1F* h_Electron_deltaPhiSCtr_after_Zto2e = new TH1F("Electron_deltaPhiSCtr", "Electron_deltaPhiSCtr", 100u, -0.05, 0.05);
    TH1F* h_Electron_dr03EcalRecHitSumEt_after_Zto2e = new TH1F("Electron_dr03EcalRecHitSumEt", "Electron_dr03EcalRecHitSumEt", 100u, 0., 3.);
    TH1F* h_Electron_dr03EcalRecHitSumEtOld_after_Zto2e = new TH1F("Electron_dr03EcalRecHitSumEtOld", "Electron_dr03EcalRecHitSumEtOld", 100u, 0., 5.);
    TH1F* h_Electron_dr03HcalDepth1TowerSumEt_after_Zto2e = new TH1F("Electron_dr03HcalDepth1TowerSumEt", "Electron_dr03HcalDepth1TowerSumEt", 100u, 0., 1.);
    TH1F* h_Electron_dr03HcalDepth1TowerSumEtOld_after_Zto2e = new TH1F("Electron_dr03HcalDepth1TowerSumEtOld", "Electron_dr03HcalDepth1TowerSumEtOld", 100u, 0., 1.);
    TH1F* h_Electron_dr03HcalTowerSumEt_after_Zto2e = new TH1F("Electron_dr03HcalTowerSumEt", "Electron_dr03HcalTowerSumEt", 100u, 0., 1.5);
    TH1F* h_Electron_dr03TkSumPt_after_Zto2e = new TH1F("Electron_dr03TkSumPt", "Electron_dr03TkSumPt", 100u, 0., 1.5);
    TH1F* h_Electron_dr03TkSumPtOld_after_Zto2e = new TH1F("Electron_dr03TkSumPtOld", "Electron_dr03TkSumPtOld", 100u, 0., 10.);
    TH1F* h_Electron_dxy_after_Zto2e = new TH1F("Electron_dxy", "Electron_dxy", 100u, -0.03, 2.);
    TH1F* h_Electron_dxyErr_after_Zto2e = new TH1F("Electron_dxyErr", "Electron_dxyErr", 100u, 0., 0.01);
    TH1F* h_Electron_dz_after_Zto2e = new TH1F("Electron_dz", "Electron_dz", 100u, -0.04, 0.04);
    TH1F* h_Electron_dzErr_after_Zto2e = new TH1F("Electron_dzErr", "Electron_dzErr", 100u, 0., 0.01);
    TH1F* h_Electron_eInvMinusPInv_after_Zto2e = new TH1F("Electron_eInvMinusPInv", "Electron_eInvMinusPInv", 100u, -0.02, 0.01);
    TH1F* h_Electron_eInvMinusPInvOld_after_Zto2e = new TH1F("Electron_eInvMinusPInvOld", "Electron_eInvMinusPInvOld", 100u, -0.025, 0.025);
    TH1F* h_Electron_eta_after_Zto2e = new TH1F("Electron_eta", "Electron_eta", 100u, -3.14, 3.14);
    TH1F* h_Electron_genPartIdx_after_Zto2e = new TH1F("Electron_genPartIdx", "Electron_genPartIdx", 100u, -2., 2.);
    TH1F* h_Electron_hoe_after_Zto2e = new TH1F("Electron_hoe", "Electron_hoe", 100u, 0., 0.04);
    TH1F* h_Electron_ip3d_after_Zto2e = new TH1F("Electron_ip3d", "Electron_ip3d", 100u, 0., 0.08);
    TH1F* h_Electron_isEB_after_Zto2e = new TH1F("Electron_isEB", "Electron_isEB", 100u, 0., 1.04);
    TH1F* h_Electron_isEE_after_Zto2e = new TH1F("Electron_isEE", "Electron_isEE", 100u, 0., 1.04);
    TH1F* h_Electron_isNano_after_Zto2e = new TH1F("Electron_isNano", "Electron_isNano", 100u, 0., 1.04);
    TH1F* h_Electron_isPFcand_after_Zto2e = new TH1F("Electron_isPFcand", "Electron_isPFcand", 100u, 0., 1.04);
    TH1F* h_Electron_lostHits_after_Zto2e = new TH1F("Electron_lostHits", "Electron_lostHits", 100u, 0., 0.16);
    TH1F* h_Electron_mass_after_Zto2e = new TH1F("Electron_mass", "Electron_mass", 100u, -0.08, 0.08);
    TH1F* h_Electron_nNano_after_Zto2e = new TH1F("Electron_nNano", "Electron_nNano", 100u, 2., 3.);
    TH1F* h_Electron_pfRelIso03_all_after_Zto2e = new TH1F("Electron_pfRelIso03_all", "Electron_pfRelIso03_all", 100u, 0., 0.2);
    TH1F* h_Electron_pfRelIso03_chg_after_Zto2e = new TH1F("Electron_pfRelIso03_chg", "Electron_pfRelIso03_chg", 100u, 0., 0.1);
    TH1F* h_Electron_phi_after_Zto2e = new TH1F("Electron_phi", "Electron_phi", 100u, -3.14, 3.14);
    TH1F* h_Electron_pt_after_Zto2e = new TH1F("Electron_pt", "Electron_pt", 100u, 10., 70.);
    TH1F* h_Electron_sieie_after_Zto2e = new TH1F("Electron_sieie", "Electron_sieie", 100u, 0.005, 0.03);
    TH1F* h_Electron_sieieR1_after_Zto2e = new TH1F("Electron_sieieR1", "Electron_sieieR1", 100u, 0.005, 0.03);
    TH1F* h_Electron_simId_after_Zto2e = new TH1F("Electron_simId", "Electron_simId", 100u, -1.1, -0.9);
    TH1F* h_Electron_sip3d_after_Zto2e = new TH1F("Electron_sip3d", "Electron_sip3d", 100u, 0., 3.5);
    TH1F* h_Electron_tightCharge_after_Zto2e = new TH1F("Electron_tightCharge", "Electron_tightCharge", 100u, 0., 2.);
    TH1F* h_Electron_vtxIdx_after_Zto2e = new TH1F("Electron_vtxIdx", "Electron_vtxIdx", 100u, 0, 0);
    TH1F* h_Electron_x_after_Zto2e = new TH1F("Electron_x", "Electron_x", 100u, 0.055, 0.09);
    TH1F* h_Electron_y_after_Zto2e = new TH1F("Electron_y", "Electron_y", 100u, 0.025, 0.055);
    TH1F* h_Electron_z_after_Zto2e = new TH1F("Electron_z", "Electron_z", 100u, -15., 15.);
    TH1F* h_PV_x_after_Zto2e = new TH1F("PV_x", "PV_x", 100u, -1., 1.);
    TH1F* h_PV_y_after_Zto2e = new TH1F("PV_y", "PV_y", 100u, -1., 1.);
    TH1F* h_PV_z_after_Zto2e = new TH1F("PV_z", "PV_z", 100u, -1., 1.);
    TH1F* h_PV_chi2_after_Zto2e = new TH1F("PV_chi2", "PV_chi2", 100u, 0., 10.);



    h_Electron_SCeta_after_Zto2e->GetXaxis()->SetTitle("Electron_SCeta");
    h_Electron_charge_after_Zto2e->GetXaxis()->SetTitle("Electron_charge");
    h_Electron_convDcot_after_Zto2e->GetXaxis()->SetTitle("Electron_convDcot");
    h_Electron_convDist_after_Zto2e->GetXaxis()->SetTitle("Electron_convDist");
    h_Electron_convVeto_after_Zto2e->GetXaxis()->SetTitle("Electron_convVeto");
    h_Electron_convVetoOld_after_Zto2e->GetXaxis()->SetTitle("Electron_convVetoOld");
    h_Electron_cutBased_after_Zto2e->GetXaxis()->SetTitle("Electron_cutBased");
    h_Electron_deltaEtaSC_after_Zto2e->GetXaxis()->SetTitle("Electron_deltaEtaSC");
    h_Electron_deltaEtaSCtr_after_Zto2e->GetXaxis()->SetTitle("Electron_deltaEtaSCtr");
    h_Electron_deltaPhiSC_after_Zto2e->GetXaxis()->SetTitle("Electron_deltaPhiSC");
    h_Electron_deltaPhiSCtr_after_Zto2e->GetXaxis()->SetTitle("Electron_deltaPhiSCtr");
    h_Electron_dr03EcalRecHitSumEt_after_Zto2e->GetXaxis()->SetTitle("Electron_dr03EcalRecHitSumEt");
    h_Electron_dr03EcalRecHitSumEtOld_after_Zto2e->GetXaxis()->SetTitle("Electron_dr03EcalRecHitSumEtOld");
    h_Electron_dr03HcalDepth1TowerSumEt_after_Zto2e->GetXaxis()->SetTitle("Electron_dr03HcalDepth1TowerSumEt");
    h_Electron_dr03HcalDepth1TowerSumEtOld_after_Zto2e->GetXaxis()->SetTitle("Electron_dr03HcalDepth1TowerSumEtOld");
    h_Electron_dr03HcalTowerSumEt_after_Zto2e->GetXaxis()->SetTitle("Electron_dr03HcalTowerSumEt");
    h_Electron_dr03TkSumPt_after_Zto2e->GetXaxis()->SetTitle("Electron_dr03TkSumPt");
    h_Electron_dr03TkSumPtOld_after_Zto2e->GetXaxis()->SetTitle("Electron_dr03TkSumPtOld");
    h_Electron_dxy_after_Zto2e->GetXaxis()->SetTitle("Electron_dxy");
    h_Electron_dxyErr_after_Zto2e->GetXaxis()->SetTitle("Electron_dxyErr");
    h_Electron_dz_after_Zto2e->GetXaxis()->SetTitle("Electron_dz");
    h_Electron_dzErr_after_Zto2e->GetXaxis()->SetTitle("Electron_dzErr");
    h_Electron_eInvMinusPInv_after_Zto2e->GetXaxis()->SetTitle("Electron_eInvMinusPInv");
    h_Electron_eInvMinusPInvOld_after_Zto2e->GetXaxis()->SetTitle("Electron_eInvMinusPInvOld");
    h_Electron_eta_after_Zto2e->GetXaxis()->SetTitle("Electron_eta");
    h_Electron_genPartIdx_after_Zto2e->GetXaxis()->SetTitle("Electron_genPartIdx");
    h_Electron_hoe_after_Zto2e->GetXaxis()->SetTitle("Electron_hoe");
    h_Electron_ip3d_after_Zto2e->GetXaxis()->SetTitle("Electron_ip3d");
    h_Electron_isEB_after_Zto2e->GetXaxis()->SetTitle("Electron_isEB");
    h_Electron_isEE_after_Zto2e->GetXaxis()->SetTitle("Electron_isEE");
    h_Electron_isNano_after_Zto2e->GetXaxis()->SetTitle("Electron_isNano");
    h_Electron_isPFcand_after_Zto2e->GetXaxis()->SetTitle("Electron_isPFcand");
    h_Electron_lostHits_after_Zto2e->GetXaxis()->SetTitle("Electron_lostHits");
    h_Electron_mass_after_Zto2e->GetXaxis()->SetTitle("Electron_mass");
    h_Electron_nNano_after_Zto2e->GetXaxis()->SetTitle("Electron_nNano");
    h_Electron_pfRelIso03_all_after_Zto2e->GetXaxis()->SetTitle("Electron_pfRelIso03_all");
    h_Electron_pfRelIso03_chg_after_Zto2e->GetXaxis()->SetTitle("Electron_pfRelIso03_chg");
    h_Electron_phi_after_Zto2e->GetXaxis()->SetTitle("Electron_phi");
    h_Electron_pt_after_Zto2e->GetXaxis()->SetTitle("Electron_pt");
    h_Electron_sieie_after_Zto2e->GetXaxis()->SetTitle("Electron_sieie");
    h_Electron_sieieR1_after_Zto2e->GetXaxis()->SetTitle("Electron_sieieR1");
    h_Electron_simId_after_Zto2e->GetXaxis()->SetTitle("Electron_simId");
    h_Electron_sip3d_after_Zto2e->GetXaxis()->SetTitle("Electron_sip3d");
    h_Electron_tightCharge_after_Zto2e->GetXaxis()->SetTitle("Electron_tightCharge");
    h_Electron_vtxIdx_after_Zto2e->GetXaxis()->SetTitle("Electron_vtxIdx");
    h_Electron_x_after_Zto2e->GetXaxis()->SetTitle("Electron_x");
    h_Electron_y_after_Zto2e->GetXaxis()->SetTitle("Electron_y");
    h_Electron_z_after_Zto2e->GetXaxis()->SetTitle("Electron_z");



    h_Electron_SCeta_after_Zto2e->GetYaxis()->SetTitle("Number of Events");
    h_Electron_charge_after_Zto2e->GetYaxis()->SetTitle("Number of Events");
    h_Electron_convDcot_after_Zto2e->GetYaxis()->SetTitle("Number of Events");
    h_Electron_convDist_after_Zto2e->GetYaxis()->SetTitle("Number of Events");
    h_Electron_convVeto_after_Zto2e->GetYaxis()->SetTitle("Number of Events");
    h_Electron_convVetoOld_after_Zto2e->GetYaxis()->SetTitle("Number of Events");
    h_Electron_cutBased_after_Zto2e->GetYaxis()->SetTitle("Number of Events");
    h_Electron_deltaEtaSC_after_Zto2e->GetYaxis()->SetTitle("Number of Events");
    h_Electron_deltaEtaSCtr_after_Zto2e->GetYaxis()->SetTitle("Number of Events");
    h_Electron_deltaPhiSC_after_Zto2e->GetYaxis()->SetTitle("Number of Events");
    h_Electron_deltaPhiSCtr_after_Zto2e->GetYaxis()->SetTitle("Number of Events");
    h_Electron_dr03EcalRecHitSumEt_after_Zto2e->GetYaxis()->SetTitle("Number of Events");
    h_Electron_dr03EcalRecHitSumEtOld_after_Zto2e->GetYaxis()->SetTitle("Number of Events");
    h_Electron_dr03HcalDepth1TowerSumEt_after_Zto2e->GetYaxis()->SetTitle("Number of Events");
    h_Electron_dr03HcalDepth1TowerSumEtOld_after_Zto2e->GetYaxis()->SetTitle("Number of Events");
    h_Electron_dr03HcalTowerSumEt_after_Zto2e->GetYaxis()->SetTitle("Number of Events");
    h_Electron_dr03TkSumPt_after_Zto2e->GetYaxis()->SetTitle("Number of Events");
    h_Electron_dr03TkSumPtOld_after_Zto2e->GetYaxis()->SetTitle("Number of Events");
    h_Electron_dxy_after_Zto2e->GetYaxis()->SetTitle("Number of Events");
    h_Electron_dxyErr_after_Zto2e->GetYaxis()->SetTitle("Number of Events");
    h_Electron_dz_after_Zto2e->GetYaxis()->SetTitle("Number of Events");
    h_Electron_dzErr_after_Zto2e->GetYaxis()->SetTitle("Number of Events");
    h_Electron_eInvMinusPInv_after_Zto2e->GetYaxis()->SetTitle("Number of Events");
    h_Electron_eInvMinusPInvOld_after_Zto2e->GetYaxis()->SetTitle("Number of Events");
    h_Electron_eta_after_Zto2e->GetYaxis()->SetTitle("Number of Events");
    h_Electron_genPartIdx_after_Zto2e->GetYaxis()->SetTitle("Number of Events");
    h_Electron_hoe_after_Zto2e->GetYaxis()->SetTitle("Number of Events");
    h_Electron_ip3d_after_Zto2e->GetYaxis()->SetTitle("Number of Events");
    h_Electron_isEB_after_Zto2e->GetYaxis()->SetTitle("Number of Events");
    h_Electron_isEE_after_Zto2e->GetYaxis()->SetTitle("Number of Events");
    h_Electron_isNano_after_Zto2e->GetYaxis()->SetTitle("Number of Events");
    h_Electron_isPFcand_after_Zto2e->GetYaxis()->SetTitle("Number of Events");
    h_Electron_lostHits_after_Zto2e->GetYaxis()->SetTitle("Number of Events");
    h_Electron_mass_after_Zto2e->GetYaxis()->SetTitle("Number of Events");
    h_Electron_nNano_after_Zto2e->GetYaxis()->SetTitle("Number of Events");
    h_Electron_pfRelIso03_all_after_Zto2e->GetYaxis()->SetTitle("Number of Events");
    h_Electron_pfRelIso03_chg_after_Zto2e->GetYaxis()->SetTitle("Number of Events");
    h_Electron_phi_after_Zto2e->GetYaxis()->SetTitle("Number of Events");
    h_Electron_pt_after_Zto2e->GetYaxis()->SetTitle("Number of Events");
    h_Electron_sieie_after_Zto2e->GetYaxis()->SetTitle("Number of Events");
    h_Electron_sieieR1_after_Zto2e->GetYaxis()->SetTitle("Number of Events");
    h_Electron_simId_after_Zto2e->GetYaxis()->SetTitle("Number of Events");
    h_Electron_sip3d_after_Zto2e->GetYaxis()->SetTitle("Number of Events");
    h_Electron_tightCharge_after_Zto2e->GetYaxis()->SetTitle("Number of Events");
    h_Electron_vtxIdx_after_Zto2e->GetYaxis()->SetTitle("Number of Events");
    h_Electron_x_after_Zto2e->GetYaxis()->SetTitle("Number of Events");
    h_Electron_y_after_Zto2e->GetYaxis()->SetTitle("Number of Events");
    h_Electron_z_after_Zto2e->GetYaxis()->SetTitle("Number of Events");


    ////////////////////////
    // Data DoubleMu 2011 //
    ////////////////////////

    // Set all the status tree to 0  

    inputChain->SetBranchStatus("Electron_SCeta", 1);
    inputChain->SetBranchStatus("Electron_charge", 1);
    inputChain->SetBranchStatus("Electron_convDcot", 1);
    inputChain->SetBranchStatus("Electron_convDist", 1);
    inputChain->SetBranchStatus("Electron_convVeto", 1);
    inputChain->SetBranchStatus("Electron_convVetoOld", 1);
    inputChain->SetBranchStatus("Electron_cutBased", 1);
    inputChain->SetBranchStatus("Electron_deltaEtaSC", 1);
    inputChain->SetBranchStatus("Electron_deltaEtaSCtr", 1);
    inputChain->SetBranchStatus("Electron_deltaPhiSC", 1);
    inputChain->SetBranchStatus("Electron_deltaPhiSCtr", 1);
    inputChain->SetBranchStatus("Electron_dr03EcalRecHitSumEt", 1);
    inputChain->SetBranchStatus("Electron_dr03EcalRecHitSumEtOld", 1);
    inputChain->SetBranchStatus("Electron_dr03HcalDepth1TowerSumEt", 1);
    inputChain->SetBranchStatus("Electron_dr03HcalDepth1TowerSumEtOld", 1);
    inputChain->SetBranchStatus("Electron_dr03HcalTowerSumEt", 1);
    inputChain->SetBranchStatus("Electron_dr03TkSumPt", 1);
    inputChain->SetBranchStatus("Electron_dr03TkSumPtOld", 1);
    inputChain->SetBranchStatus("Electron_dxy", 1);
    inputChain->SetBranchStatus("Electron_dxyErr", 1);
    inputChain->SetBranchStatus("Electron_dz", 1);
    inputChain->SetBranchStatus("Electron_dzErr", 1);
    inputChain->SetBranchStatus("Electron_eInvMinusPInv", 1);
    inputChain->SetBranchStatus("Electron_eInvMinusPInvOld", 1);
    inputChain->SetBranchStatus("Electron_eta", 1);
    inputChain->SetBranchStatus("Electron_genPartIdx", 1);
    inputChain->SetBranchStatus("Electron_hoe", 1);
    inputChain->SetBranchStatus("Electron_ip3d", 1);
    inputChain->SetBranchStatus("Electron_isEB", 1);
    inputChain->SetBranchStatus("Electron_isEE", 1);
    inputChain->SetBranchStatus("Electron_isNano", 1);
    inputChain->SetBranchStatus("Electron_isPFcand", 1);
    inputChain->SetBranchStatus("Electron_lostHits", 1);
    inputChain->SetBranchStatus("Electron_mass", 1);
    inputChain->SetBranchStatus("Electron_nNano", 1);
    inputChain->SetBranchStatus("Electron_pfRelIso03_all", 1);
    inputChain->SetBranchStatus("Electron_pfRelIso03_chg", 1);
    inputChain->SetBranchStatus("Electron_phi", 1);
    inputChain->SetBranchStatus("Electron_pt", 1);
    inputChain->SetBranchStatus("Electron_sieie", 1);
    inputChain->SetBranchStatus("Electron_sieieR1", 1);
    inputChain->SetBranchStatus("Electron_simId", 1);
    inputChain->SetBranchStatus("Electron_sip3d", 1);
    inputChain->SetBranchStatus("Electron_tightCharge", 1);
    inputChain->SetBranchStatus("Electron_vtxIdx", 1);
    inputChain->SetBranchStatus("Electron_x", 1);
    inputChain->SetBranchStatus("Electron_y", 1);
    inputChain->SetBranchStatus("Electron_z", 1);
    inputChain->SetBranchStatus("nElectron", 1);
    inputChain->SetBranchStatus("PV_x", 1);
    inputChain->SetBranchStatus("PV_y", 1);
    inputChain->SetBranchStatus("PV_z", 1);
    inputChain->SetBranchStatus("PV_chi2", 1);


    inputChain->SetBranchAddress("Electron_SCeta", &Electron_SCeta);
    inputChain->SetBranchAddress("Electron_charge", &Electron_charge);
    inputChain->SetBranchAddress("Electron_convDcot", &Electron_convDcot);
    inputChain->SetBranchAddress("Electron_convDist", &Electron_convDist);
    inputChain->SetBranchAddress("Electron_convVeto", &Electron_convVeto);
    inputChain->SetBranchAddress("Electron_convVetoOld", &Electron_convVetoOld);
    inputChain->SetBranchAddress("Electron_cutBased", &Electron_cutBased);
    inputChain->SetBranchAddress("Electron_deltaEtaSC", &Electron_deltaEtaSC);
    inputChain->SetBranchAddress("Electron_deltaEtaSCtr", &Electron_deltaEtaSCtr);
    inputChain->SetBranchAddress("Electron_deltaPhiSC", &Electron_deltaPhiSC);
    inputChain->SetBranchAddress("Electron_deltaPhiSCtr", &Electron_deltaPhiSCtr);
    inputChain->SetBranchAddress("Electron_dr03EcalRecHitSumEt", &Electron_dr03EcalRecHitSumEt);
    inputChain->SetBranchAddress("Electron_dr03EcalRecHitSumEtOld", &Electron_dr03EcalRecHitSumEtOld);
    inputChain->SetBranchAddress("Electron_dr03HcalDepth1TowerSumEt", &Electron_dr03HcalDepth1TowerSumEt);
    inputChain->SetBranchAddress("Electron_dr03HcalDepth1TowerSumEtOld", &Electron_dr03HcalDepth1TowerSumEtOld);
    inputChain->SetBranchAddress("Electron_dr03HcalTowerSumEt", &Electron_dr03HcalTowerSumEt);
    inputChain->SetBranchAddress("Electron_dr03TkSumPt", &Electron_dr03TkSumPt);
    inputChain->SetBranchAddress("Electron_dr03TkSumPtOld", &Electron_dr03TkSumPtOld);
    inputChain->SetBranchAddress("Electron_dxy", &Electron_dxy);
    inputChain->SetBranchAddress("Electron_dxyErr", &Electron_dxyErr);
    inputChain->SetBranchAddress("Electron_dz", &Electron_dz);
    inputChain->SetBranchAddress("Electron_dzErr", &Electron_dzErr);
    inputChain->SetBranchAddress("Electron_eInvMinusPInv", &Electron_eInvMinusPInv);
    inputChain->SetBranchAddress("Electron_eInvMinusPInvOld", &Electron_eInvMinusPInvOld);
    inputChain->SetBranchAddress("Electron_eta", &Electron_eta);
    inputChain->SetBranchAddress("Electron_genPartIdx", &Electron_genPartIdx);
    inputChain->SetBranchAddress("Electron_hoe", &Electron_hoe);
    inputChain->SetBranchAddress("Electron_ip3d", &Electron_ip3d);
    inputChain->SetBranchAddress("Electron_isEB", &Electron_isEB);
    inputChain->SetBranchAddress("Electron_isEE", &Electron_isEE);
    inputChain->SetBranchAddress("Electron_isNano", &Electron_isNano);
    inputChain->SetBranchAddress("Electron_isPFcand", &Electron_isPFcand);
    inputChain->SetBranchAddress("Electron_lostHits", &Electron_lostHits);
    inputChain->SetBranchAddress("Electron_mass", &Electron_mass);
    inputChain->SetBranchAddress("Electron_nNano", &Electron_nNano);
    inputChain->SetBranchAddress("Electron_pfRelIso03_all", &Electron_pfRelIso03_all);
    inputChain->SetBranchAddress("Electron_pfRelIso03_chg", &Electron_pfRelIso03_chg);
    inputChain->SetBranchAddress("Electron_phi", &Electron_phi);
    inputChain->SetBranchAddress("Electron_pt", &Electron_pt);
    inputChain->SetBranchAddress("Electron_sieie", &Electron_sieie);
    inputChain->SetBranchAddress("Electron_sieieR1", &Electron_sieieR1);
    inputChain->SetBranchAddress("Electron_simId", &Electron_simId);
    inputChain->SetBranchAddress("Electron_sip3d", &Electron_sip3d);
    inputChain->SetBranchAddress("Electron_tightCharge", &Electron_tightCharge);
    inputChain->SetBranchAddress("Electron_vtxIdx", &Electron_vtxIdx);
    inputChain->SetBranchAddress("Electron_x", &Electron_x);
    inputChain->SetBranchAddress("Electron_y", &Electron_y);
    inputChain->SetBranchAddress("Electron_z", &Electron_z);
    inputChain->SetBranchAddress("nElectron", &nElectron);
    inputChain->SetBranchAddress("PV_x", &PV_x);
    inputChain->SetBranchAddress("PV_y", &PV_y);
    inputChain->SetBranchAddress("PV_z", &PV_z);
    inputChain->SetBranchAddress("PV_chi2", &PV_chi2);


    //activating branches for debugging
    inputChain->SetBranchStatus("run", 1);
    inputChain->SetBranchStatus("event", 1);
    inputChain->SetBranchAddress("run", &run);  
    inputChain->SetBranchAddress("event", &event);
  
    // ********************************* //  


    Int_t nevent = inputChain->GetEntries();

    for (int aa = 0; aa < nevent; aa++) {
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

                // h_p_e->Fill(ele.P());
                // h_et_e->Fill(ele.Et());
                // h_pt_e_b4->Fill(Electron_pt[ee]);
                // h_eta_e_b4->Fill(Electron_eta[ee]);
                // h_sc_eta->Fill(Electron_deltaEtaSC[ee] + Electron_eta[ee]);
                // h_phi_e->Fill(Electron_phi[ee]);
                // h_misshite->Fill(Electron_lostHits[ee]);
                // h_relPFIso_e->Fill(Electron_pfRelIso03_all[ee]);
                // h_dxy_e->Fill(Electron_dxy[ee]);
                // h_SIP3d_e_b4->Fill(Electron_sip3d[ee]);
    
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
                if (S >= 80. && S <= 100.) {

                        h_PV_x_after_Zto2e->Fill(PV_x);
                        h_PV_y_after_Zto2e->Fill(PV_y);
                        h_PV_z_after_Zto2e->Fill(PV_z);
                        h_PV_chi2_after_Zto2e->Fill(PV_chi2);

                    for (int ee = 0; ee < nGoodElectrons; ee++) {

                        // h_pt_e_after_Zto2e->Fill(Electron_pt[goodElectrons.at(ee).first]);
                        // h_eta_e_after_Zto2e->Fill(Electron_eta[goodElectrons.at(ee).first]);
                        h_Electron_SCeta_after_Zto2e->Fill(Electron_SCeta[goodElectrons.at(ee).first]);
                        h_Electron_charge_after_Zto2e->Fill(Electron_charge[goodElectrons.at(ee).first]);
                        h_Electron_convDcot_after_Zto2e->Fill(Electron_convDcot[goodElectrons.at(ee).first]);
                        h_Electron_convDist_after_Zto2e->Fill(Electron_convDist[goodElectrons.at(ee).first]);
                        h_Electron_convVeto_after_Zto2e->Fill(Electron_convVeto[goodElectrons.at(ee).first]);
                        h_Electron_convVetoOld_after_Zto2e->Fill(Electron_convVetoOld[goodElectrons.at(ee).first]);
                        h_Electron_cutBased_after_Zto2e->Fill(Electron_cutBased[goodElectrons.at(ee).first]);
                        h_Electron_deltaEtaSC_after_Zto2e->Fill(Electron_deltaEtaSC[goodElectrons.at(ee).first]);
                        h_Electron_deltaEtaSCtr_after_Zto2e->Fill(Electron_deltaEtaSCtr[goodElectrons.at(ee).first]);
                        h_Electron_deltaPhiSC_after_Zto2e->Fill(Electron_deltaPhiSC[goodElectrons.at(ee).first]);
                        h_Electron_deltaPhiSCtr_after_Zto2e->Fill(Electron_deltaPhiSCtr[goodElectrons.at(ee).first]);
                        h_Electron_dr03EcalRecHitSumEt_after_Zto2e->Fill(Electron_dr03EcalRecHitSumEt[goodElectrons.at(ee).first]);
                        h_Electron_dr03EcalRecHitSumEtOld_after_Zto2e->Fill(Electron_dr03EcalRecHitSumEtOld[goodElectrons.at(ee).first]);
                        h_Electron_dr03HcalDepth1TowerSumEt_after_Zto2e->Fill(Electron_dr03HcalDepth1TowerSumEt[goodElectrons.at(ee).first]);
                        h_Electron_dr03HcalDepth1TowerSumEtOld_after_Zto2e->Fill(Electron_dr03HcalDepth1TowerSumEtOld[goodElectrons.at(ee).first]);
                        h_Electron_dr03HcalTowerSumEt_after_Zto2e->Fill(Electron_dr03HcalTowerSumEt[goodElectrons.at(ee).first]);
                        h_Electron_dr03TkSumPt_after_Zto2e->Fill(Electron_dr03TkSumPt[goodElectrons.at(ee).first]);
                        h_Electron_dr03TkSumPtOld_after_Zto2e->Fill(Electron_dr03TkSumPtOld[goodElectrons.at(ee).first]);
                        h_Electron_dxy_after_Zto2e->Fill(Electron_dxy[goodElectrons.at(ee).first]);
                        h_Electron_dxyErr_after_Zto2e->Fill(Electron_dxyErr[goodElectrons.at(ee).first]);
                        h_Electron_dz_after_Zto2e->Fill(Electron_dz[goodElectrons.at(ee).first]);
                        h_Electron_dzErr_after_Zto2e->Fill(Electron_dzErr[goodElectrons.at(ee).first]);
                        h_Electron_eInvMinusPInv_after_Zto2e->Fill(Electron_eInvMinusPInv[goodElectrons.at(ee).first]);
                        h_Electron_eInvMinusPInvOld_after_Zto2e->Fill(Electron_eInvMinusPInvOld[goodElectrons.at(ee).first]);
                        h_Electron_eta_after_Zto2e->Fill(Electron_eta[goodElectrons.at(ee).first]);
                        h_Electron_genPartIdx_after_Zto2e->Fill(Electron_genPartIdx[goodElectrons.at(ee).first]);
                        h_Electron_hoe_after_Zto2e->Fill(Electron_hoe[goodElectrons.at(ee).first]);
                        h_Electron_ip3d_after_Zto2e->Fill(Electron_ip3d[goodElectrons.at(ee).first]);
                        h_Electron_isEB_after_Zto2e->Fill(Electron_isEB[goodElectrons.at(ee).first]);
                        h_Electron_isEE_after_Zto2e->Fill(Electron_isEE[goodElectrons.at(ee).first]);
                        h_Electron_isNano_after_Zto2e->Fill(Electron_isNano[goodElectrons.at(ee).first]);
                        h_Electron_isPFcand_after_Zto2e->Fill(Electron_isPFcand[goodElectrons.at(ee).first]);
                        h_Electron_lostHits_after_Zto2e->Fill(Electron_lostHits[goodElectrons.at(ee).first]);
                        h_Electron_mass_after_Zto2e->Fill(Electron_mass[goodElectrons.at(ee).first]);
                        h_Electron_nNano_after_Zto2e->Fill(Electron_nNano[goodElectrons.at(ee).first]);
                        h_Electron_pfRelIso03_all_after_Zto2e->Fill(Electron_pfRelIso03_all[goodElectrons.at(ee).first]);
                        h_Electron_pfRelIso03_chg_after_Zto2e->Fill(Electron_pfRelIso03_chg[goodElectrons.at(ee).first]);
                        h_Electron_phi_after_Zto2e->Fill(Electron_phi[goodElectrons.at(ee).first]);
                        h_Electron_pt_after_Zto2e->Fill(Electron_pt[goodElectrons.at(ee).first]);
                        h_Electron_sieie_after_Zto2e->Fill(Electron_sieie[goodElectrons.at(ee).first]);
                        h_Electron_sieieR1_after_Zto2e->Fill(Electron_sieieR1[goodElectrons.at(ee).first]);
                        h_Electron_simId_after_Zto2e->Fill(Electron_simId[goodElectrons.at(ee).first]);
                        h_Electron_sip3d_after_Zto2e->Fill(Electron_sip3d[goodElectrons.at(ee).first]);
                        h_Electron_tightCharge_after_Zto2e->Fill(Electron_tightCharge[goodElectrons.at(ee).first]);
                        h_Electron_vtxIdx_after_Zto2e->Fill(Electron_vtxIdx[goodElectrons.at(ee).first]);
                        h_Electron_x_after_Zto2e->Fill(Electron_x[goodElectrons.at(ee).first]);
                        h_Electron_y_after_Zto2e->Fill(Electron_y[goodElectrons.at(ee).first]);
                        h_Electron_z_after_Zto2e->Fill(Electron_z[goodElectrons.at(ee).first]);
                        // h_nElectron_after_Zto2e->Fill(nElectron[goodElectrons.at(ee).first]);

                    }
                }

            } // charge 1 + charge 2 = 0

        } // nGoodElectrons >= 2

    }

    TFile outputFile((nameOutFile + ".root").c_str(), "RECREATE");
    // TDirectory* dir = outputFile.mkdir("Events");
    // dir->cd();

    // Electron selection
    // h_p_e->Write();
    // h_et_e->Write();
    // h_pt_e_b4->Write();
    // h_eta_e_b4->Write();
    // h_sc_eta->Write();
    // h_phi_e->Write();
    // h_misshite->Write();
    // h_relPFIso_e->Write();
    // h_dxy_e->Write();
    // h_SIP3d_e_b4->Write();
    h_nge->Write();

    // Z to 2e
    h_mZ_2e->Write();
    h_Electron_SCeta_after_Zto2e->Write();
    h_Electron_charge_after_Zto2e->Write();
    h_Electron_convDcot_after_Zto2e->Write();
    h_Electron_convDist_after_Zto2e->Write();
    h_Electron_convVeto_after_Zto2e->Write();
    h_Electron_convVetoOld_after_Zto2e->Write();
    h_Electron_cutBased_after_Zto2e->Write();
    h_Electron_deltaEtaSC_after_Zto2e->Write();
    h_Electron_deltaEtaSCtr_after_Zto2e->Write();
    h_Electron_deltaPhiSC_after_Zto2e->Write();
    h_Electron_deltaPhiSCtr_after_Zto2e->Write();
    h_Electron_dr03EcalRecHitSumEt_after_Zto2e->Write();
    h_Electron_dr03EcalRecHitSumEtOld_after_Zto2e->Write();
    h_Electron_dr03HcalDepth1TowerSumEt_after_Zto2e->Write();
    h_Electron_dr03HcalDepth1TowerSumEtOld_after_Zto2e->Write();
    h_Electron_dr03HcalTowerSumEt_after_Zto2e->Write();
    h_Electron_dr03TkSumPt_after_Zto2e->Write();
    h_Electron_dr03TkSumPtOld_after_Zto2e->Write();
    h_Electron_dxy_after_Zto2e->Write();
    h_Electron_dxyErr_after_Zto2e->Write();
    h_Electron_dz_after_Zto2e->Write();
    h_Electron_dzErr_after_Zto2e->Write();
    h_Electron_eInvMinusPInv_after_Zto2e->Write();
    h_Electron_eInvMinusPInvOld_after_Zto2e->Write();
    h_Electron_eta_after_Zto2e->Write();
    h_Electron_genPartIdx_after_Zto2e->Write();
    h_Electron_hoe_after_Zto2e->Write();
    h_Electron_ip3d_after_Zto2e->Write();
    h_Electron_isEB_after_Zto2e->Write();
    h_Electron_isEE_after_Zto2e->Write();
    h_Electron_isNano_after_Zto2e->Write();
    h_Electron_isPFcand_after_Zto2e->Write();
    h_Electron_lostHits_after_Zto2e->Write();
    h_Electron_mass_after_Zto2e->Write();
    h_Electron_nNano_after_Zto2e->Write();
    h_Electron_pfRelIso03_all_after_Zto2e->Write();
    h_Electron_pfRelIso03_chg_after_Zto2e->Write();
    h_Electron_phi_after_Zto2e->Write();
    h_Electron_pt_after_Zto2e->Write();
    h_Electron_sieie_after_Zto2e->Write();
    h_Electron_sieieR1_after_Zto2e->Write();
    h_Electron_simId_after_Zto2e->Write();
    h_Electron_sip3d_after_Zto2e->Write();
    h_Electron_tightCharge_after_Zto2e->Write();
    h_Electron_vtxIdx_after_Zto2e->Write();
    h_Electron_x_after_Zto2e->Write();
    h_Electron_y_after_Zto2e->Write();
    h_Electron_z_after_Zto2e->Write();

    outputFile.Close();
    gROOT->ProcessLine(".q");

}


int main() {

    cms_data();
    return 0;

}