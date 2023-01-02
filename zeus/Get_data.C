// Electrons and muons are always sorted from higher to lower pT
/*
2021 DESY Summer Student Program - fixed bugs from 2019 project. 
This code is for the 2011 Run A nanoAODplus ntuples. 
Same code can be run on any NanoAODplus ntuple . 
*/


#include <stdio.h>
#ifdef WINDOWS
   #include <direct.h>
   #define Define_CurrentDir _getcwd
#else
   #include <unistd.h>
   #define Define_CurrentDir getcwd
#endif
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

#include <iostream>
#include <iomanip>      // std::setprecision
#include <sstream>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>

using namespace std;

void Get_data(string in_file = "", string dir_out= "", string name_out_file = "") {

  gROOT->Reset();
  gStyle->SetOptStat("nemruo"); 

/*
  // Path to nanoAOD files
  string inDir = "/nfs/dust/cms/user/zulaiha/PhD/all_output/";

  // nanoAOD files
  // Data 2011
  string DataMu2011 = "2011/Data/Data11_DoubleMuRunA_ZeroBias2e.root";
  string DataE2011 = "2011/Data/Data11_DoubleERunA_ZeroBias2e.root";

  TChain * tDataMu2011 = new TChain("Events");
  tDataMu2011->Add((inDir + DataMu2011).c_str());

  TChain * tDataE2011 = new TChain("Events");
  tDataE2011->Add((inDir + DataE2011).c_str());

  // MC 2011
  string MCHZZ2011 = "2011/MC/MC11_HZZ_ZeroBias2e.root";
  string MCZZto4mu2011 = "2011/MC/MC11_ZZto4mu_ZeroBias2e.root";
  string MCZZto4e2011 = "2011/MC/MC11_ZZto4e_ZeroBias2e.root";
  string MCZZto2mu2e2011 = "2011/MC/MC11_ZZto2mu2e_ZeroBias2e.root";

  TChain * tMCHZZ2011 = new TChain("Events");
  tMCHZZ2011->Add((inDir + MCHZZ2011).c_str());

  TChain * tMCZZto4mu2011 = new TChain("Events");
  tMCZZto4mu2011->Add((inDir + MCZZto4mu2011).c_str());

  TChain * tMCZZto4e2011 = new TChain("Events");
  tMCZZto4e2011->Add((inDir + MCZZto4e2011).c_str());

  TChain * tMCZZto2mu2e2011 = new TChain("Events");
  tMCZZto2mu2e2011->Add((inDir + MCZZto2mu2e2011).c_str());

  // Data 2012
  string DataMu2012_1 = "2012/Data/Data12_DoubleMuParkedRunB_ZeroBias2e_0to4.root";
  string DataMu2012_2 = "2012/Data/Data12_DoubleMuParkedRunB_ZeroBias2e_5to9.root";
  string DataMu2012_3 = "2012/Data/Data12_DoubleMuParkedRunC_ZeroBias2e_0to2.root";
  string DataMu2012_4 = "2012/Data/Data12_DoubleMuParkedRunC_ZeroBias2e_3to9.root";
  string DataE2012_1 = "2012/Data/Data12_DoubleEParkedRunB_ZeroBias2e_0to4.root";
  string DataE2012_2 = "2012/Data/Data12_DoubleEParkedRunB_ZeroBias2e_5to9.root";
  string DataE2012_3 = "2012/Data/Data12_DoubleEParkedRunC_ZeroBias2e_0to4.root";
  string DataE2012_4 = "2012/Data/Data12_DoubleEParkedRunC_ZeroBias2e_5to9.root";

  TChain * tDataMu2012 = new TChain("Events");
  tDataMu2012->Add((inDir + DataMu2012_1).c_str());
  tDataMu2012->Add((inDir + DataMu2012_2).c_str());
  tDataMu2012->Add((inDir + DataMu2012_3).c_str());
  tDataMu2012->Add((inDir + DataMu2012_4).c_str());

  TChain * tDataE2012 = new TChain("Events");
  tDataE2012->Add((inDir + DataE2012_1).c_str());
  tDataE2012->Add((inDir + DataE2012_2).c_str());
  tDataE2012->Add((inDir + DataE2012_3).c_str());
  tDataE2012->Add((inDir + DataE2012_4).c_str());

  // MC 2012
  string MCHZZ2012 = "2012/MC/MC12_HZZ_ZeroBias2e.root";
  string MCZZto4mu2012 = "2012/MC/MC12_ZZto4mu_ZeroBias2e.root";
  string MCZZto4e2012 = "2012/MC/MC12_ZZto4e_ZeroBias2e.root";
  string MCZZto2mu2e2012 = "2012/MC/MC12_ZZto2mu2e_ZeroBias2e.root";

  TChain * tMCHZZ2012 = new TChain("Events");
  tMCHZZ2012->Add((inDir + MCHZZ2012).c_str());

  TChain * tMCZZto4mu2012 = new TChain("Events");
  tMCZZto4mu2012->Add((inDir + MCZZto4mu2012).c_str());

  TChain * tMCZZto4e2012 = new TChain("Events");
  tMCZZto4e2012->Add((inDir + MCZZto4e2012).c_str());

  TChain * tMCZZto2mu2e2012 = new TChain("Events");
  tMCZZto2mu2e2012->Add((inDir + MCZZto2mu2e2012).c_str());

  // TChains
  TChain * nanoAODs[12] = {tDataMu2011, tDataE2011, tMCHZZ2011, tMCZZto4mu2011, tMCZZto4e2011, tMCZZto2mu2e2011,
                       tDataMu2012, tDataE2012, tMCHZZ2012, tMCZZto4mu2012, tMCZZto4e2012, tMCZZto2mu2e2012};

  // Directory names in the .root out file
  string DirNames[12] = {"DataDoubleMu2011", "DataDooubleE2011", "MCHZZ2011", "MCZZto4mu2011", "MCZZto4e2011", "MCZZto2mu2e2011",
                         "DataDoubleMu2012", "DataDoubleE2012", "MCHZZ2012", "MCZZto4mu2012", "MCZZto4e2012", "MCZZto2mu2e2012"};

*/
  
   ////////////////////////// Write histograms ////////////////////////////
    char LOCAL_DIR[FILENAME_MAX];
    std::cout << "Directory where the code are being running: " << LOCAL_DIR<<endl;
    string topDir = LOCAL_DIR + dir_out;

  
  // ** C. Check subdirectory structure for requested options and create directories if necessary
	// * i. Check if topdir exists
	struct stat sb;
	if (!(stat(topDir.c_str(), &sb) == 0 && S_ISDIR(sb.st_mode))){
		cout << "top-level director, " << topDir << " , DNE. Creating now" << endl;
		mkdir(topDir.c_str(), S_IRWXU | S_IRWXG | S_IROTH);
	}

	// * ii. Create subdir
	std::string sampleDir = "/plots";

	topDir = (topDir + sampleDir + "/").c_str();
	if (!(stat(topDir.c_str(), &sb) == 0 && S_ISDIR(sb.st_mode))){
		cout << "sample subdirectory , " << topDir << " , DNE. Creating now" << endl;
		mkdir(topDir.c_str(), S_IRWXU | S_IRWXG | S_IROTH);
	}

  
  // TChain
  TChain * t = new TChain("Events");
  t->Add((in_file).c_str());  
  //t->Add((inDir + DataMu2011B).c_str());
  
  ///////////////////////
  // Declare variables //
  ///////////////////////
  int aflag=0,bflag=0;
  // Muon
  UInt_t nMuon;
  float Muon_pt[128];
  float Muon_eta[128];
  float Muon_phi[128];
  float Muon_mass[128];
  int Muon_charge[128];
  float Muon_dxy[128]; 
  float Muon_dz[128];
  float Muon_sip3d[128]; 
  float Muon_pfRelIso04_all[128];
  bool Muon_isPFcand[128];
  bool Muon_isGlobal[128];

  vector< pair<int, float> > goodMuons;

  TLorentzVector muon1, muon2, muon3, muon4;

  // Electron
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

  TLorentzVector ele1, ele2, ele3, ele4;

  // Other
  double S;
  double S1;
  double S2;
  double SQM;

  double mZ = 91.1876;
  double sqm1 = (0.105658) * (0.105658);
  double sqme = (0.0005109989) * (0.0005109989);

  double s1, s2, s3, s4, s;
  double dx,dy,dz, rap, pz;

  double mZ12, mZ34, mZ13, mZ24, mZ14, mZ23;
  double mass4mu, pt_4mu, eta_4mu, phi_4mu;
  double px4mu, py4mu, pz4mu, E4mu;
  double pt_mu1, pt_mu2, pt_mu3, pt_mu4;
  double eta_mu1, eta_mu2, eta_mu3, eta_mu4;
  double phi_mu1, phi_mu2, phi_mu3, phi_mu4;
  int cas_mu1, cas_mu2, cas_mu3, cas_mu4;

  double px_mu1, px_mu2, px_mu3, px_mu4;
  double py_mu1, py_mu2, py_mu3, py_mu4;
  double pz_mu1, pz_mu2, pz_mu3, pz_mu4;

  double E_mu1, E_mu2, E_mu3, E_mu4;

  double mZa, mZb;

  double eZ12, eZ34, eZ13, eZ24, eZ14, eZ23; 
  double pxZ12, pxZ34, pxZ13, pxZ24,  pxZ14, pxZ23;
  double pyZ12, pyZ34, pyZ13, pyZ24, pyZ14, pyZ23;
  double pzZ12, pzZ34, pzZ13, pzZ24, pzZ14, pzZ23;
  double pZ12, pZ34, pZ13, pZ24, pZ14, pZ23;
  double pTZ12, pTZ34, pTZ13, pTZ24, pTZ14, pTZ23;

  double dZ12, dZ34, dZ13, dZ24, dZ14, dZ23;
  double dZc1, dZc2, dZc3;
  double eZa, pxZa, pyZa, pzZa, pTZa;
  double eZb, pxZb, pyZb, pzZb, pTZb;

  TLorentzVector p4Za, p4Zb, p4H;

  double mass4e, pt_4e, eta_4e, phi_4e;
  double px4e, py4e, pz4e, E4e;
  double pt_e1, pt_e2, pt_e3, pt_e4;
  double eta_e1, eta_e2, eta_e3, eta_e4;
  double phi_e1, phi_e2, phi_e3, phi_e4;
  int cas_e1, cas_e2, cas_e3, cas_e4;

  double px_e1, px_e2, px_e3, px_e4;
  double py_e1, py_e2, py_e3, py_e4;
  double pz_e1, pz_e2, pz_e3, pz_e4;

  double E_e1, E_e2, E_e3, E_e4;

  double mass2mu2e, pt_2mu2e, eta_2mu2e, phi_2mu2e;
  double px2mu2e, py2mu2e, pz2mu2e, E2mu2e;
  double pt_2mu1, pt_2mu2, pt_2e1, pt_2e2;
  double eta_2mu1, eta_2mu2, eta_2e1, eta_2e2;
  double phi_2mu1, phi_2mu2, phi_2e1, phi_2e2;
  int cas_2mu1, cas_2mu2, cas_2e1, cas_2e2;

  double px_2mu1, px_2mu2, px_2e1, px_2e2;
  double py_2mu1, py_2mu2, py_2e1, py_2e2;
  double pz_2mu1, pz_2mu2, pz_2e1, pz_2e2;

  double E_2mu1, E_2mu2, E_2e1, E_2e2;

  // Declare histograms

  // No. of good reco muon 
  TH1F *h_ngmu = new TH1F("NGoodRecMuons", "No. of good reco muon", 10, 0., 10.);
  h_ngmu->GetXaxis()->SetTitle("Number of RMuons");
  h_ngmu->GetYaxis()->SetTitle("Number of Events");

  // No. of good electron
  TH1F *h_nge = new TH1F("NGoodElectron", "No. of good electron", 10, 0., 10.);
  h_nge->GetXaxis()->SetTitle("Number of Electrons");
  h_nge->GetYaxis()->SetTitle("Number of Events");

  // Z to 2mu -------------------------------------------------------------------

  // ZTo2mu mass spectrum with reco muon
  TH1F *h_mZ_2mu = new TH1F("massZto2muon", "mass of Z to 2 muon",120, 40., 120.);
  h_mZ_2mu->GetXaxis()->SetTitle("Invariant mass of dimuon (in GeV/c^2)");
  h_mZ_2mu->GetYaxis()->SetTitle("Number of Events");

  // Momentum of Reco Muon 
  TH1F *h_p_reco = new TH1F("RM_momentum", "RM Momentum", 200, 0., 200.);
  h_p_reco->GetXaxis()->SetTitle("Momentum (GeV/c)");
  h_p_reco->GetYaxis()->SetTitle("Number of Events");

  // pT of Reco Muon
  TH1F *h_pt_reco_b4 = new TH1F("b4_RM_pt", "RM pT", 200, 0., 200.);
  h_pt_reco_b4->GetXaxis()->SetTitle("pT (GeV/c)");
  h_pt_reco_b4->GetYaxis()->SetTitle("Number of Events");

  // Eta of Reco Muon
  TH1F *h_eta_reco_b4 = new TH1F("b4_RM_eta", "RM Eta", 140, -3.5, 3.5);
  h_eta_reco_b4->GetXaxis()->SetTitle("eta");
  h_eta_reco_b4->GetYaxis()->SetTitle("Number of Events");

  // Phi of Reco Muon
  TH1F *h_phi_reco = new TH1F("RM_phi", "RM Phi", 314, -3.15, 3.15);
  h_phi_reco->GetXaxis()->SetTitle("Phi");
  h_phi_reco->GetYaxis()->SetTitle("Number of Events");

  // Transverse impact parameter with respect to primary vertex of Reco Muon
  TH1F *h_dxy_mu = new TH1F("RM_dxy", "RM dxy", 100, 0., 1.);
  h_dxy_mu->GetXaxis()->SetTitle("Transverse impact parameter w.r.t. BS");
  h_dxy_mu->GetYaxis()->SetTitle("Number of Events");

  // Relative isolation of Reco Muon
  TH1F *h_relPFIso_mu = new TH1F("RM_RelPFIso", "R.PFIso", 100, 0., 5.);
  h_relPFIso_mu->GetXaxis()->SetTitle("Relative Isolation mu");
  h_relPFIso_mu->GetYaxis()->SetTitle("Number of Events");

  // Relative isolation of Reco Muon after cuts
  TH1F *h_relPFIso_mu_after = new TH1F("after_RM_RelPFIso", "R.PFIso", 100, 0., 5.);
  h_relPFIso_mu_after->GetXaxis()->SetTitle("Relative Isolation mu");
  h_relPFIso_mu_after->GetYaxis()->SetTitle("Number of Events");

  // pT reco muon after cuts for Z to 2mu
  TH1F *h_pt_after_Zto2mu = new TH1F("after_RM_pt_Z2mu", "Mu pT", 200, 0., 200.);
  h_pt_after_Zto2mu->GetXaxis()->SetTitle("pT (GeV/c)");
  h_pt_after_Zto2mu->GetYaxis()->SetTitle("Number of Events");

  // Eta reco muon after cuts for Z to 2mu
  TH1F *h_eta_after_Zto2mu = new TH1F("after_RM_eta_Z2mu","Mu eta", 140, -3.5, 3.5);
  h_eta_after_Zto2mu->GetXaxis()->SetTitle("Eta");
  h_eta_after_Zto2mu->GetYaxis()->SetTitle("Number of Events");

  // pT reco muon after cuts
  TH1F *h_pt_after = new TH1F("after_RM_pt", "Muon pT", 200, 0., 200.); 
  h_pt_after->GetXaxis()->SetTitle("pT (GeV/c)");
  h_pt_after->GetYaxis()->SetTitle("Number of Events");

  // Eta reco muon after cuts
  TH1F *h_eta_after = new TH1F("after_RM_eta", "Muon eta", 140, -3.5, 3.5);
  h_eta_after->GetXaxis()->SetTitle("eta");
  h_eta_after->GetYaxis()->SetTitle("Number of Events");

  // 2mu2e case
  // Relative isolation of 2mu for 2mu2e after cuts
  TH1F *h_relPFIso_2mu_after = new TH1F("after_relPFIso_2mu","R.PFIso", 50, 0., 5.);
  h_relPFIso_2mu_after->GetXaxis()->SetTitle("Relative Isolation mu");
  h_relPFIso_2mu_after->GetYaxis()->SetTitle("Number of Events");

  // Relative isolation of 2e for 2mu2e after cuts
  TH1F *h_relPFIso_2e_after = new TH1F("after_relPFIso_2e","R.PFIso", 50, 0., 5.);
  h_relPFIso_2e_after->GetXaxis()->SetTitle("Relative Isolation e");
  h_relPFIso_2e_after->GetYaxis()->SetTitle("Number of Events");

  // pT muon after cuts for 2mu2e
  TH1F *h_pt_after_2mu2e = new TH1F("after_pt_2mu2e", "Muon pT", 200, 0., 200.); 
  h_pt_after_2mu2e->GetXaxis()->SetTitle("pT (GeV/c)");
  h_pt_after_2mu2e->GetYaxis()->SetTitle("Number of Events");

  // Eta muon after cuts for 2mu2e
  TH1F *h_eta_after_2mu2e = new TH1F("after_eta_2mu2e", "Muon eta", 140, -3.5, 3.5);
  h_eta_after_2mu2e->GetXaxis()->SetTitle("eta");
  h_eta_after_2mu2e->GetYaxis()->SetTitle("Number of Events");


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

  // ZZ to 4mu ------------------------------------------------------------------

  // These histograms are for 4muon reconstruction with different combinations
  
  // First combination: 1234
  // Mass Z12
  TH1F *h_mZ12_4mu = new TH1F("mZ12_4mu", "mass of Z12", 75, 0., 150.);
  h_mZ12_4mu->GetXaxis()->SetTitle("Invariant mass of dimuon (in GeV/c^2)");
  h_mZ12_4mu->GetYaxis()->SetTitle("Number of Events");

  // Mass Z34
  TH1F *h_mZ34_4mu = new TH1F("mZ34_4mu", "mass of Z34", 75, 0., 150.);
  h_mZ34_4mu->GetXaxis()->SetTitle("Invariant mass of dimuon (in GeV/c^2)");
  h_mZ34_4mu->GetYaxis()->SetTitle("Number of Events");

  // Second combination: 1324
  // Mass Z13
  TH1F *h_mZ13_4mu = new TH1F("mZ13_4mu", "mass of Z13", 75, 0., 150.);
  h_mZ13_4mu->GetXaxis()->SetTitle("Invariant mass of dimuon (in GeV/c^2)");
  h_mZ13_4mu->GetYaxis()->SetTitle("Number of Events");

  // Mass Z24
  TH1F *h_mZ24_4mu = new TH1F("mZ24_4mu", "mass of Z24", 75, 0., 150.);
  h_mZ24_4mu->GetXaxis()->SetTitle("Invariant mass of dimuon");
  h_mZ24_4mu->GetYaxis()->SetTitle("Number of Events");

  // Third combination: 1423
  // Mass Z14
  TH1F *h_mZ14_4mu = new TH1F("mZ14_4mu", "mass of Z14", 75, 0., 150.);
  h_mZ14_4mu->GetXaxis()->SetTitle("Invariant mass of dimuon (in GeV/c^2)");
  h_mZ14_4mu->GetYaxis()->SetTitle("Number of Events");

  // Mass Z23
  TH1F *h_mZ23_4mu = new TH1F("mZ23_4mu", "mass of Z23", 75, 0., 150.);
  h_mZ23_4mu->GetXaxis()->SetTitle("Invariant mass of dimuon (in GeV/c^2)");
  h_mZ23_4mu->GetYaxis()->SetTitle("Number of Events");

  // Mass Za: mass of ZTo2mu closest to Z mass
  TH1F *h_mZa_4mu = new TH1F("mZa_4mu", "mass Za", 120, 0., 120.);
  h_mZa_4mu->GetXaxis()->SetTitle("Invariant mass of dimuon (in GeV/c^2)");
  h_mZa_4mu->GetYaxis()->SetTitle("Number of Events");

  // Mass Zb: mass of ZTo2mu not closest to Z mass
  TH1F *h_mZb_4mu = new TH1F("mZb_4mu", "mass Zb", 120, 0., 120.);
  h_mZb_4mu->GetXaxis()->SetTitle("Invariant mass of dimuon (in GeV/c^2)");
  h_mZb_4mu->GetYaxis()->SetTitle("Number of Events");

  // 4muon mass spectrum with reco muon (paper 7 TeV)
  TH1F *h_m1_m4mu = new TH1F("mass4mu_7TeV", "mass of 4 muon", 51, 98., 608.);
  h_m1_m4mu->GetXaxis()->SetTitle("Invariant mass of 4muons (in GeV/c^2)");
  h_m1_m4mu->GetYaxis()->SetTitle("Number of Events");

  // 4muon mass spectrum with reco muon (paper 8TeV)
  TH1F *h_m2_m4mu = new TH1F("mass4mu_8TeV", "mass of 4 muon", 74, 70., 810.);
  h_m2_m4mu->GetXaxis()->SetTitle("Invariant mass of 4muons (in GeV/c^2)");
  h_m2_m4mu->GetYaxis()->SetTitle("Number of Events");

  // 4muon mass spectrum with reco muon (paper 8TeV lower range)
  TH1F *h_m3_m4mu = new TH1F("mass4mu_8TeV_low", "mass of 4 muon", 37, 70., 181.);
  h_m3_m4mu->GetXaxis()->SetTitle("Invariant mass of 4muons (in GeV/c^2)");
  h_m3_m4mu->GetYaxis()->SetTitle("Number of Events");

  // 4muon mass spectrum with reco muon (full mass range)
  TH1F *h_m4_m4mu = new TH1F("mass4mu_full", "mass of 4 muon", 300, 0., 900.);
  h_m4_m4mu->GetXaxis()->SetTitle("Invariant mass of 4muons (in GeV/c^2)");
  h_m4_m4mu->GetYaxis()->SetTitle("Number of Events");


  // ZZ to 4e --------------------------------------------------------------------

  // These histograms are for 4electron reconstruction with different combinations
  
  // First combination: 1234
  // Mass Z12
  TH1F *h_mZ12_4e = new TH1F("mZ12_4e", "mass of Z12", 75, 0., 150.);
  h_mZ12_4e->GetXaxis()->SetTitle("Invariant mass of dielectron (in GeV/c^2)");
  h_mZ12_4e->GetYaxis()->SetTitle("Number of Events");

  // Mass Z34
  TH1F *h_mZ34_4e = new TH1F("mZ34_4e", "mass of Z34", 75, 0., 150.);
  h_mZ34_4e->GetXaxis()->SetTitle("Invariant mass of dielectron (in GeV/c^2)");
  h_mZ34_4e->GetYaxis()->SetTitle("Number of Events");

  // Second combination: 1324
  // Mass Z13
  TH1F *h_mZ13_4e = new TH1F("mZ13_4e", "mass of Z13", 75, 0., 150.);
  h_mZ13_4e->GetXaxis()->SetTitle("Invariant mass of dielectron (in GeV/c^2)");
  h_mZ13_4e->GetYaxis()->SetTitle("Number of Events");

  // Mass Z24
  TH1F *h_mZ24_4e = new TH1F("mZ24_4e", "mass of Z24", 75, 0., 150.);
  h_mZ24_4e->GetXaxis()->SetTitle("Invariant mass of dielectron (in GeV/c^2)");
  h_mZ24_4e->GetYaxis()->SetTitle("Number of Events");

  // Third combination: 1423
  // Mass Z14
  TH1F *h_mZ14_4e = new TH1F("mZ14_4e", "mass of Z14", 75, 0., 150.);
  h_mZ14_4e->GetXaxis()->SetTitle("Invariant mass of dielectron (in GeV/c^2)");
  h_mZ14_4e->GetYaxis()->SetTitle("Number of Events");

  // Mass Z23
  TH1F *h_mZ23_4e = new TH1F("mZ23_4e", "mass of Z23", 75, 0., 150.);
  h_mZ23_4e->GetXaxis()->SetTitle("Invariant mass of dielectron (in GeV/c^2)");
  h_mZ23_4e->GetYaxis()->SetTitle("Number of Events");

  // Mass Za: mass of Z closest to Z mass
  TH1F *h_mZa_4e = new TH1F("mZa_4e", "mass Za", 120, 0., 120.);
  h_mZa_4e->GetXaxis()->SetTitle("Invariant mass of dielectron (in GeV/c^2)");
  h_mZa_4e->GetYaxis()->SetTitle("Number of Events");

  // Mass Zb: mass of Z not closest to Z mass
  TH1F *h_mZb_4e = new TH1F("mZb_4e", "mass Zb", 120, 0., 120.);
  h_mZb_4e->GetXaxis()->SetTitle("Invariant mass of dielectron (in GeV/c^2)");
  h_mZb_4e->GetYaxis()->SetTitle("Number of Events");

  // 4electron mass spectrum with Gsf Electron (paper 7TeV)
  TH1F *h_m1_m4e = new TH1F("mass4e_7TeV", "mass of 4 electron", 51, 98., 608.);
  h_m1_m4e->GetXaxis()->SetTitle("Invariant mass of 4e (in GeV/c^2)");
  h_m1_m4e->GetYaxis()->SetTitle("Number of Events");

  // 4electron mass spectrum with Gsf Electron (paper 8TeV)
  TH1F *h_m2_m4e = new TH1F("mass4e_8TeV", "mass of 4 electron", 74, 70., 810.);
  h_m2_m4e->GetXaxis()->SetTitle("Invariant mass of 4e (in GeV/c^2)");
  h_m2_m4e->GetYaxis()->SetTitle("Number of Events");

  // 4electron mass spectrum with Gsf Electron (paper 8TeV lower range)
  TH1F *h_m3_m4e = new TH1F("mass4e_8TeV_low", "mass 4 electron", 37, 70., 181.);
  h_m3_m4e->GetXaxis()->SetTitle("Invariant mass of 4e (in GeV/c^2)");
  h_m3_m4e->GetYaxis()->SetTitle("Number of Events");

  // 4electron mass spectrum with Gsf Electron (full mass range)
  TH1F *h_m4_m4e = new TH1F("mass4e_full", "mass of 4 electron", 300, 0., 900.);
  h_m4_m4e->GetXaxis()->SetTitle("Invariant mass of 4e (in GeV/c^2)");
  h_m4_m4e->GetYaxis()->SetTitle("Number of Events");

  // ZZ to 2mu 2e ----------------------------------------------------------------

  // These histograms are for 2mu2e reconstruction with different combinations
  
  // Mass of Z to 2mu from 2mu2e
  TH1F *h_mZmu_2mu2e = new TH1F("massZmu_2mu2e", "mass Z2mu:2mu2e", 75, 0., 150.);
  h_mZmu_2mu2e->GetXaxis()->SetTitle("Invariant mass of Z1 (in GeV/c^2)");
  h_mZmu_2mu2e->GetYaxis()->SetTitle("Number of Events");

  // Mass of Z to 2e from 2mu2e
  TH1F *h_mZe_2mu2e = new TH1F("massZe_2mu2e", "mass Z2e:2mu2e", 75, 0., 150.);
  h_mZe_2mu2e->GetXaxis()->SetTitle("Invariant Mass of Z2 (in GeV/c^2)");
  h_mZe_2mu2e->GetYaxis()->SetTitle("Number of Events");

  // Mass Za: mass of Z1 closest to Z mass
  TH1F *h_mZa_2mu2e = new TH1F("mZa_2mu2e", "mass Z higher", 120, 0., 120.);
  h_mZa_2mu2e->GetXaxis()->SetTitle("Invariant mass of dimuon (in GeV/c^2)");
  h_mZa_2mu2e->GetYaxis()->SetTitle("Number of Events");

  // Mass Zb: mass of Z2 not closest to Z mass
  TH1F *h_mZb_2mu2e = new TH1F("mZb_2mu2e", "mass Z lower", 120, 0., 120.);
  h_mZb_2mu2e->GetXaxis()->SetTitle("Invariant mass of dimuon (in GeV/c^2)");
  h_mZb_2mu2e->GetYaxis()->SetTitle("Number of Events");

  // 2muon 2electron mass spectrum (paper 7TeV)
  TH1F *h_m1_m2mu2e = new TH1F("mass2mu2e_7TeV", "mass of 2mu2e", 51, 98., 608.);
  h_m1_m2mu2e->GetXaxis()->SetTitle("Invariant mass of 2mu2e (in GeV/c^2)");
  h_m1_m2mu2e->GetYaxis()->SetTitle("Number of Events");
   
  // These histograms are all for 7 TeV case. 

  // 2muon 2electron mass spectrum (paper 8TeV)
  TH1F *h_m2_m2mu2e = new TH1F("mass2mu2e_8TeV", "mass of 2mu2e", 74, 70., 810.);
  h_m2_m2mu2e->GetXaxis()->SetTitle("Invariant mass of 2mu2e (in GeV/c^2)");
  h_m2_m2mu2e->GetYaxis()->SetTitle("Number of Events");

  // 2muon 2electron mass spectrum (paper 8TeV lower range)
  TH1F *h_m3_m2mu2e = new TH1F("mass2mu2e_8TeV_low", "mass 2mu2e", 37, 70., 181.);
  h_m3_m2mu2e->GetXaxis()->SetTitle("Invariant mass of 2mu2e (in GeV/c^2)");
  h_m3_m2mu2e->GetYaxis()->SetTitle("Number of Events");

  // 2muons 2electrons mass spectrum full mass range
  TH1F *h_m4_m2mu2e = new TH1F("mass2mu2e_full", "mass of 2mu2e", 300, 0., 900.);
  h_m4_m2mu2e->GetXaxis()->SetTitle("Invariant mass of 2mu2e (in GeV/c^2)");
  h_m4_m2mu2e->GetYaxis()->SetTitle("Number of Events");
  UInt_t run;
  ULong64_t event;

  ////////////////////////
  // Data DoubleMu 2011 //
  ////////////////////////

  // Set all the status tree to 0  
  t->SetBranchStatus("*", 0);

  // Activate the variables we want to use
  t->SetBranchStatus("nMuon", 1);  
  t->SetBranchStatus("Muon_pt", 1); 
  t->SetBranchStatus("Muon_eta", 1);
  t->SetBranchStatus("Muon_phi", 1);
  t->SetBranchStatus("Muon_mass", 1);
  t->SetBranchStatus("Muon_charge", 1);
  t->SetBranchStatus("Muon_dxy", 1);
  t->SetBranchStatus("Muon_dz", 1);
  t->SetBranchStatus("Muon_sip3d", 1);
  t->SetBranchStatus("Muon_pfRelIso04_all", 1);
  t->SetBranchStatus("Muon_isPFcand", 1);
  t->SetBranchStatus("Muon_isGlobal", 1);

  t->SetBranchStatus("nElectron", 1);  
  t->SetBranchStatus("Electron_pt", 1); 
  t->SetBranchStatus("Electron_eta", 1);
  t->SetBranchStatus("Electron_deltaEtaSC", 1);
  t->SetBranchStatus("Electron_phi", 1);
  t->SetBranchStatus("Electron_mass", 1);
  t->SetBranchStatus("Electron_charge", 1);
  t->SetBranchStatus("Electron_dxy", 1);
  t->SetBranchStatus("Electron_dz", 1);
  t->SetBranchStatus("Electron_sip3d", 1);
  t->SetBranchStatus("Electron_pfRelIso03_all", 1);
  t->SetBranchStatus("Electron_isPFcand", 1);
  t->SetBranchStatus("Electron_lostHits", 1);
  t->SetBranchStatus("Electron_isEE", 1);
  t->SetBranchStatus("Electron_isEB", 1);

  t->SetBranchAddress("nMuon", &nMuon);  
  t->SetBranchAddress("Muon_pt", Muon_pt); 
  t->SetBranchAddress("Muon_eta", Muon_eta);
  t->SetBranchAddress("Muon_phi", Muon_phi);
  t->SetBranchAddress("Muon_mass", Muon_mass);
  t->SetBranchAddress("Muon_charge", Muon_charge);
  t->SetBranchAddress("Muon_dxy", Muon_dxy);
  t->SetBranchAddress("Muon_dz", Muon_dz);
  t->SetBranchAddress("Muon_sip3d", Muon_sip3d);
  t->SetBranchAddress("Muon_pfRelIso04_all", Muon_pfRelIso04_all);
  t->SetBranchAddress("Muon_isPFcand", Muon_isPFcand);
  t->SetBranchAddress("Muon_isGlobal", Muon_isGlobal);

  t->SetBranchAddress("nElectron", &nElectron);  
  t->SetBranchAddress("Electron_pt", Electron_pt); 
  t->SetBranchAddress("Electron_eta", Electron_eta);
  t->SetBranchAddress("Electron_deltaEtaSC", Electron_deltaEtaSC);
  t->SetBranchAddress("Electron_phi", Electron_phi);
  t->SetBranchAddress("Electron_mass", Electron_mass);
  t->SetBranchAddress("Electron_charge", Electron_charge);
  t->SetBranchAddress("Electron_dxy", Electron_dxy);
  t->SetBranchAddress("Electron_dz", Electron_dz);
  t->SetBranchAddress("Electron_sip3d", Electron_sip3d);
  t->SetBranchAddress("Electron_pfRelIso03_all", Electron_pfRelIso03_all);
  t->SetBranchAddress("Electron_isPFcand", Electron_isPFcand);
  t->SetBranchAddress("Electron_lostHits", Electron_lostHits);
  t->SetBranchAddress("Electron_isEE", Electron_isEE);
  t->SetBranchAddress("Electron_isEB", Electron_isEB);
  

  //activating branches for debugging
  t->SetBranchStatus("run", 1);
  t->SetBranchStatus("event", 1);
  t->SetBranchAddress("run", &run);  
  t->SetBranchAddress("event", &event);  
  
  // ********************************* //  

  Int_t nevent = t->GetEntries();

  // Muon and electron selection (because of ZZ to 2mu 2e)

  for (int aa = 0; aa < nevent; aa++) {
  //  for (int aa=12000000;aa<13500000;aa++) {
    // if (aa > 100) {break;}
    if (aa % 100000 == 0) {cout << "loop " << aa << " / " << nevent << endl;}

    t->GetEntry(aa);

    // Muons ----------------------------------------------------------------------------

    int nGoodMuons = 0;
    goodMuons.clear();
    int good_flag=0;
    // Loop over all muons/event
    for (int mm = 0; mm < nMuon; mm++) {

      TLorentzVector muon;
      muon.SetPtEtaPhiM(Muon_pt[mm], Muon_eta[mm], Muon_phi[mm], Muon_mass[mm]);

      // Global particle flow muons
      // isPFIsolationValid() not defined in NanoAnalyzer but implicitly implemented
      if (Muon_isPFcand[mm] && Muon_isGlobal[mm]) {

        h_p_reco->Fill(muon.P());
        h_pt_reco_b4->Fill(Muon_pt[mm]);
        h_eta_reco_b4->Fill(Muon_eta[mm]);
        h_phi_reco->Fill(Muon_phi[mm]);
        h_relPFIso_mu->Fill(Muon_pfRelIso04_all[mm]);
        h_dxy_mu->Fill(Muon_dxy[mm]);
        h_SIP3d_mu_b4->Fill(Muon_sip3d[mm]);

        // Muon selection | according to report, reliso=1.0
        if (abs(Muon_sip3d[mm]) < 4. && abs(Muon_dxy[mm]) < 0.5 && abs(Muon_dz[mm]) < 1. && Muon_pfRelIso04_all[mm] < 0.4) { // 0.4 or 1.0

          if (Muon_pt[mm] > 5. && abs(Muon_eta[mm]) < 2.4) {
            good_flag=1;
          /*  if ((run==163071&&event==11989798) || (run==172992&&event==1153485608)|| (run==173657&&event==34442568))
               cout << run << "," << event << "," << mm << "," << nMuon << "," << abs(Muon_sip3d[mm]) << "," << abs(Muon_dxy[mm]) << "," <<  Muon_pfRelIso04_all[mm] << "," << Muon_pt[mm] << "," << abs(Muon_eta[mm]) << "," << Muon_charge[mm] << "," << good_flag << endl;
              cout << Muon_pfRelIso04_all[mm]<<endl; } */
            goodMuons.push_back( make_pair(mm, Muon_pt[mm]) );
            nGoodMuons++;

          }

        } // Muon selection

      } // Global particle flow muons

    } // Loop over all muons/event


    h_ngmu->Fill(nGoodMuons);


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

                  goodElectrons.push_back( make_pair(ee, Electron_pt[ee]) );
                  nGoodElectrons++;

                }
              }

              else if (Electron_isEE[ee]) {

                if (Electron_pfRelIso03_all[ee] < 0.4) {

                  goodElectrons.push_back( make_pair(ee, Electron_pt[ee]) );
                  nGoodElectrons++;

                }
              }  // EB/EE
            }  // dxy, dz
          }  // lost hits, SIP 3D
        }  // pT, SC eta
      }  // Particle flow electrons
    } // Loop over all electrons/event


    h_nge->Fill(nGoodElectrons);


    /////////////////////////// Z to 2mu ////////////////////////////

    // if (nGoodMuons >= 2) {

    //   int m1 = goodMuons.at(0).first;
    //   int m2 = goodMuons.at(1).first;
     
    //   if (Muon_charge[m1] + Muon_charge[m2] == 0) { 

    //     muon1.SetPtEtaPhiM(Muon_pt[m1], Muon_eta[m1], Muon_phi[m1], Muon_mass[m1]);
    //     muon2.SetPtEtaPhiM(Muon_pt[m2], Muon_eta[m2], Muon_phi[m2], Muon_mass[m2]);

    //     float SQM = (0.105658) * (0.105658);
    //     float S1 = sqrt( (muon1.P()*muon1.P() + SQM) * (muon2.P()*muon2.P() + SQM) );
    //     float S2 = muon1.Px() * muon2.Px() + muon1.Py() * muon2.Py() + muon1.Pz() * muon2.Pz();
    //     float S = sqrt(2. * (SQM + S1 - S2));

    //     h_mZ_2mu->Fill(S);

    //     // pT and eta after the cuts
    //     for (int mm = 0; mm < nGoodMuons; mm++) {

    //       h_pt_after_Zto2mu->Fill(Muon_pt[goodMuons.at(mm).first]);
    //       h_eta_after_Zto2mu->Fill(Muon_eta[goodMuons.at(mm).first]);

    //     }

    //   } // charge 1 + charge 2 = 0

    // } // nGoodMuons >= 2

    /////////////////////////// Z to 2e /////////////////////////////

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

    ////////////////////////// ZZ to 4mu ////////////////////////////
    bool combo1234=false,combo1324=false,combo1423=false;
    
    if (nGoodMuons >= 4) {

      int m1 = goodMuons.at(0).first;       
      int m2 = goodMuons.at(1).first;      
      int m3 = goodMuons.at(2).first;
      int m4 = goodMuons.at(3).first;
      if (Muon_charge[m1] + Muon_charge[m2] + Muon_charge[m3] + Muon_charge[m4]== 0) { 

        muon1.SetPtEtaPhiM(Muon_pt[m1], Muon_eta[m1], Muon_phi[m1], Muon_mass[m1]);
        muon2.SetPtEtaPhiM(Muon_pt[m2], Muon_eta[m2], Muon_phi[m2], Muon_mass[m2]);
        muon3.SetPtEtaPhiM(Muon_pt[m3], Muon_eta[m3], Muon_phi[m3], Muon_mass[m3]);
        muon4.SetPtEtaPhiM(Muon_pt[m4], Muon_eta[m4], Muon_phi[m4], Muon_mass[m4]);

	// First combination: Combine muon 1234
	if (Muon_charge[m1] + Muon_charge[m2] == 0) {
      combo1234=true;
	  eZ12 = (sqrt(muon1.P() * muon1.P() + sqm1)) +
		 (sqrt(muon2.P() * muon2.P() + sqm1));
	
	  pxZ12 = muon1.Px() + muon2.Px();
	  pyZ12 = muon1.Py() + muon2.Py();
	  pzZ12 = muon1.Pz() + muon2.Pz();

	  if (Muon_charge[m3] + Muon_charge[m4] == 0) {
            eZ34 = (sqrt(muon3.P() * muon3.P() + sqm1)) +
		   (sqrt(muon4.P() * muon4.P() + sqm1));

            pxZ34 = muon3.Px() + muon4.Px();
	    pyZ34 = muon3.Py() + muon4.Py();
	    pzZ34 = muon3.Pz() + muon4.Pz();
        
	    // Calculate p4
	    pZ12 = sqrt((pxZ12 * pxZ12) + (pyZ12 * pyZ12) + (pzZ12 * pzZ12));
	    pZ34 = sqrt((pxZ34 * pxZ34) + (pyZ34 * pyZ34) + (pzZ34 * pzZ34));

	    pTZ12 = sqrt((pxZ12 * pxZ12) + (pyZ12 * pyZ12));   //P_T2 = P_x2 + P_y2
	    pTZ34 = sqrt((pxZ34 * pxZ34) + (pyZ34 * pyZ34));

	    mZ12 = sqrt((eZ12 * eZ12) - (pZ12 * pZ12));   // m2 = E2-p2
	    mZ34 = sqrt((eZ34 * eZ34) - (pZ34 * pZ34));

	    if (mZ12 > 0.) h_mZ12_4mu->Fill(mZ12);
	    if (mZ34 > 0.) h_mZ34_4mu->Fill(mZ34);
              
          }
        }

	dZ12 = abs( mZ12 - mZ );
	dZ34 = abs( mZ34 - mZ );

	// take the smallest difference between mass
	// to use for 4muon combination
	dZc1 = (dZ12 < dZ34) ? dZ12 : dZ34; 
 
	// Second combination: Combine muon 1324
	if (Muon_charge[m1] + Muon_charge[m3] == 0) {
      combo1324=true;
	  eZ13 = (sqrt(muon1.P() * muon1.P() + sqm1)) +
		 (sqrt(muon3.P() * muon3.P() + sqm1));

	  pxZ13 = muon1.Px() + muon3.Px();
	  pyZ13 = muon1.Py() + muon3.Py();
	  pzZ13 = muon1.Pz() + muon3.Pz();

	  if (Muon_charge[m2] + Muon_charge[m4] == 0) {

            eZ24 = (sqrt(muon2.P() * muon2.P() + sqm1)) +
		   (sqrt(muon4.P() * muon4.P() + sqm1));

            pxZ24 = muon2.Px() + muon4.Px();
	    pyZ24 = muon2.Py() + muon4.Py();
	    pzZ24 = muon2.Pz() + muon4.Pz();

	    // Calculate p4
	    pZ13 = sqrt((pxZ13 * pxZ13) + (pyZ13 * pyZ13) + (pzZ13 * pzZ13));
	    pZ24 = sqrt((pxZ24 * pxZ24) + (pyZ24 * pyZ24) + (pzZ24 * pzZ24));

	    pTZ13 = sqrt((pxZ13 * pxZ13) + (pyZ13 * pyZ13));
	    pTZ24 = sqrt((pxZ24 * pxZ24) + (pyZ24 * pyZ24));

	    mZ13 = sqrt((eZ13 * eZ13) - (pZ13 * pZ13));
	    mZ24 = sqrt((eZ24 * eZ24) - (pZ24 * pZ24));

	    if (mZ13 > 0.) h_mZ13_4mu->Fill(mZ13);
            if (mZ24 > 0.) h_mZ24_4mu->Fill(mZ24);
	  }
	}

	dZ13 = abs( mZ13 - mZ );
	dZ24 = abs( mZ24 - mZ );

	dZc2 = (dZ13 < dZ24) ? dZ13 : dZ24; 

	// Third combination: Combine muon 1423
	if (Muon_charge[m1] + Muon_charge[m4] == 0) {
      combo1423=true;
	  eZ14 = (sqrt(muon1.P() * muon1.P() + sqm1)) +
		 (sqrt(muon4.P() * muon4.P() + sqm1));

	  pxZ14 = muon1.Px() + muon4.Px();
	  pyZ14 = muon1.Py() + muon4.Py();
	  pzZ14 = muon1.Pz() + muon4.Pz();

	  if (Muon_charge[m2] + Muon_charge[m3] == 0) {

            eZ23 = (sqrt(muon2.P() * muon2.P() + sqm1)) +
		   (sqrt(muon3.P() * muon3.P() + sqm1));
       
	    pxZ23 = muon2.Px() + muon3.Px();
	    pyZ23 = muon2.Py() + muon3.Py();
	    pzZ23 = muon2.Pz() + muon3.Pz();

	    // Calculate p4
	    pZ14 = sqrt((pxZ14 * pxZ14) + (pyZ14 * pyZ14) + (pzZ14 * pzZ14));
	    pZ23 = sqrt((pxZ23 * pxZ23) + (pyZ23 * pyZ23) + (pzZ23 * pzZ23));

	    pTZ14 = sqrt((pxZ14 * pxZ14) + (pyZ14 * pyZ14));
	    pTZ23 = sqrt((pxZ23 * pxZ23) + (pyZ23 * pyZ23));

	    mZ14 = sqrt((eZ14 * eZ14) - (pZ14 * pZ14));
	    mZ23 = sqrt((eZ23 * eZ23) - (pZ23 * pZ23));
  
            if (mZ14 > 0.) h_mZ14_4mu->Fill(mZ14);
	    if (mZ23 > 0.) h_mZ23_4mu->Fill(mZ23);

          }
        }

	dZ14 = abs( mZ14 - mZ );
	dZ23 = abs( mZ23 - mZ );
  
	dZc3 = (dZ14 < dZ23) ? dZ14 : dZ23;
  
  
  /*if ((run==163071&&event==11989798) || (run==172992&&event==1153485608) || (run==173657&&event==34442568))
   { cout<<dZc1<<","<<dZc2<<","<<dZc3<<","<<combo1234<<","<<combo1324<<","<<combo1423<<endl;
   cout << " " << endl; }*/


	bool ptZadaug = false;


/*####################### Logic Overhaul starts from here ################################*/


    if (!combo1324) {   // This means combo 13-24 has like charges and isn't possible 

	if (dZc1 < dZc3) {       // This means combo 12-34 is what we need even though combo 14-23 also has opposite charges, but 12-34 is closer in mass

	  if (dZ12 < dZ34) {
	    eZa  = eZ12;     
	    pxZa = pxZ12;
	    pyZa = pyZ12;
	    pzZa = pzZ12;
	    pTZa = pTZ12;
	    mZa  = mZ12;

	    if (muon1.Pt() > 20. and muon2.Pt() > 10.)
	    ptZadaug = true;

	    eZb  = eZ34;
	    pxZb = pxZ34;
	    pyZb = pyZ34;
	    pzZb = pzZ34;
	    pTZb = pTZ34;
	    mZb  = mZ34;
 
	  }
	      
          else {

	    eZa  = eZ34;
	    pxZa = pxZ34;
	    pyZa = pyZ34;
	    pzZa = pzZ34;
	    pTZa = pTZ34;
	    mZa  = mZ34;

	    if (muon3.Pt() > 20. and muon4.Pt() > 10.)
	    ptZadaug = true;

	    eZb  = eZ12;
	    pxZb = pxZ12;
	    pyZb = pyZ12;
	    pzZb = pzZ12;
	    pTZb = pTZ12;
	    mZb  = mZ12;

          }
        }
    else {  // This is for the case when 14-23 is closer in mass
    
       if (dZ14 < dZ23) {

	    eZa  = eZ14;
	    pxZa = pxZ14;
	    pyZa = pyZ14;
	    pzZa = pzZ14;
	    pTZa = pTZ14;
	    mZa  = mZ14;

	    if (muon1.Pt() > 20. and muon4.Pt() > 10.)
	    ptZadaug = true;

	    eZb  = eZ23;
	    pxZb = pxZ23;
	    pyZb = pyZ23;
	    pzZb = pzZ23;
	    pTZb = pTZ23;
	    mZb  = mZ23;

	  }

	  else {

	    eZa  = eZ23;
	    pxZa = pxZ23;
	    pyZa = pyZ23;
	    pzZa = pzZ23;
	    pTZa = pTZ23;
	    mZa  = mZ23;

	    if (muon2.Pt() > 20. and muon3.Pt() > 10.)
	    ptZadaug = true;

	    eZb  = eZ14;
	    pxZb = pxZ14;
	    pyZb = pyZ14;
	    pzZb = pzZ14;
	    pTZb = pTZ14;
	    mZb  = mZ14;

	  }
    } 
    }

    else if (!combo1234) {   // This means combo 12-34 does not work because of like charges

	 if (dZc2 < dZc3) {  // This means out of 13-24 and 14-23, 13-24 is closer in mass

	  if (dZ13 < dZ24) {    //Out of muon pairs 13 and 24, 13 is closer in mass

	    eZa  = eZ13;
	    pxZa = pxZ13;
	    pyZa = pyZ13;
	    pzZa = pzZ13;
	    pTZa = pTZ13;
	    mZa  = mZ13;

	    if (muon1.Pt() > 20. and muon3.Pt() > 10.)
	    ptZadaug = true;

	    eZb  = eZ24;
	    pxZb = pxZ24;
	    pyZb = pyZ24;
	    pzZb = pzZ24;
	    pTZb = pTZ24;
	    mZb  = mZ24;

	  }
	      
          else {   // out of muon pairs 13 and 24, 24 is closer in mass

	    eZa  = eZ24;
	    pxZa = pxZ24;
	    pyZa = pyZ24;
	    pzZa = pzZ24;
	    pTZa = pTZ24;
	    mZa  = mZ24;

	    if (muon2.Pt() > 20. and muon4.Pt() > 10.)
	    ptZadaug = true;

	    eZb  = eZ13;
	    pxZb = pxZ13;
	    pyZb = pyZ13;
	    pzZb = pzZ13;
	    pTZb = pTZ13;
	    mZb  = mZ13;

	  }
	}
    else {
        if (dZ14 < dZ23) {

	    eZa  = eZ14;
	    pxZa = pxZ14;
	    pyZa = pyZ14;
	    pzZa = pzZ14;
	    pTZa = pTZ14;
	    mZa  = mZ14;

	    if (muon1.Pt() > 20. and muon4.Pt() > 10.)
	    ptZadaug = true;

	    eZb  = eZ23;
	    pxZb = pxZ23;
	    pyZb = pyZ23;
	    pzZb = pzZ23;
	    pTZb = pTZ23;
	    mZb  = mZ23;

	  }

	  else {

	    eZa  = eZ23;
	    pxZa = pxZ23;
	    pyZa = pyZ23;
	    pzZa = pzZ23;
	    pTZa = pTZ23;
	    mZa  = mZ23;

	    if (muon2.Pt() > 20. and muon3.Pt() > 10.)
	    ptZadaug = true;

	    eZb  = eZ14;
	    pxZb = pxZ14;
	    pyZb = pyZ14;
	    pzZb = pzZ14;
	    pTZb = pTZ14;
	    mZb  = mZ14;

	  }
    }
    }

    else {    // This only leaves the case where combo 1423 is not possible, now we check against the combos 1324 and 1234

	 if (dZc1 < dZc2) {    // for the case where combo 12-34 works
  
	  if (dZ12 < dZ34) {    // In combo 12-34, muon pair 12 is closer in mass
	    eZa  = eZ12;     
	    pxZa = pxZ12;
	    pyZa = pyZ12;
	    pzZa = pzZ12;
	    pTZa = pTZ12;
	    mZa  = mZ12;

	    if (muon1.Pt() > 20. and muon2.Pt() > 10.)
	    ptZadaug = true;

	    eZb  = eZ34;
	    pxZb = pxZ34;
	    pyZb = pyZ34;
	    pzZb = pzZ34;
	    pTZb = pTZ34;
	    mZb  = mZ34;
 
	  }
	      
          else {   // Same but now muon pair 34 is closer in mass

	    eZa  = eZ34;
	    pxZa = pxZ34;
	    pyZa = pyZ34;
	    pzZa = pzZ34;
	    pTZa = pTZ34;
	    mZa  = mZ34;

	    if (muon3.Pt() > 20. and muon4.Pt() > 10.)
	    ptZadaug = true;

	    eZb  = eZ12;
	    pxZb = pxZ12;
	    pyZb = pyZ12;
	    pzZb = pzZ12;
	    pTZb = pTZ12;
	    mZb  = mZ12;

          }
     }
     else   // for the case where combo 13-24 works
     {
       if (dZ13 < dZ24) {    //In combo 13-24, muon pair 13 is closer in mass

	    eZa  = eZ13;
	    pxZa = pxZ13;
	    pyZa = pyZ13;
	    pzZa = pzZ13;
	    pTZa = pTZ13;
	    mZa  = mZ13;

	    if (muon1.Pt() > 20. and muon3.Pt() > 10.)
	    ptZadaug = true;

	    eZb  = eZ24;
	    pxZb = pxZ24;
	    pyZb = pyZ24;
	    pzZb = pzZ24;
	    pTZb = pTZ24;
	    mZb  = mZ24;

	  }
	      
          else {   // Same but now muon pair 24 is closer in mass

	    eZa  = eZ24;
	    pxZa = pxZ24;
	    pyZa = pyZ24;
	    pzZa = pzZ24;
	    pTZa = pTZ24;
	    mZa  = mZ24;

	    if (muon2.Pt() > 20. and muon4.Pt() > 10.)
	    ptZadaug = true;

	    eZb  = eZ13;
	    pxZb = pxZ13;
	    pyZb = pyZ13;
	    pzZb = pzZ13;
	    pTZb = pTZ13;
	    mZb  = mZ13;

	  }

     }
	}
/* ######################## Logic Overhaul Ends here ############################33 */
	if (ptZadaug) {

	  if (mZa > 40. && mZa < 120.) {

	    if (mZb > 12. && mZb < 120.) {

	      h_mZa_4mu->Fill(mZa);
	      h_mZb_4mu->Fill(mZb);

	      // 4 vector
	      p4Za.SetPxPyPzE(pxZa, pyZa, pzZa, eZa);
	      p4Zb.SetPxPyPzE(pxZb, pyZb, pzZb, eZb);

	      p4H = p4Za + p4Zb;

	      mass4mu = p4H.M();
	      pt_4mu = p4H.Pt();
	      eta_4mu = p4H.Eta();
	      phi_4mu = p4H.Phi();

	      px4mu = p4H.Px();
	      py4mu = p4H.Py();
	      pz4mu = p4H.Pz();
	      E4mu = p4H.E();

	      pt_mu1 = muon1.Pt();
	      pt_mu2 = muon2.Pt();
	      pt_mu3 = muon3.Pt();
	      pt_mu4 = muon4.Pt();

	      eta_mu1 = muon1.Eta();
	      eta_mu2 = muon2.Eta();
	      eta_mu3 = muon3.Eta();
	      eta_mu4 = muon4.Eta();

	      phi_mu1 = muon1.Phi();
	      phi_mu2 = muon2.Phi();
	      phi_mu3 = muon3.Phi();
	      phi_mu4 = muon4.Phi();

	      cas_mu1 = Muon_charge[m1];
	      cas_mu2 = Muon_charge[m2];
	      cas_mu3 = Muon_charge[m3];
	      cas_mu4 = Muon_charge[m4];

	      px_mu1 = muon1.Px();
	      px_mu2 = muon2.Px();
	      px_mu3 = muon3.Px();
	      px_mu4 = muon4.Px();

	      py_mu1 = muon1.Py();
	      py_mu2 = muon2.Py();
	      py_mu3 = muon3.Py();
	      py_mu4 = muon4.Py();
		
	      pz_mu1 = muon1.Pz();
	      pz_mu2 = muon2.Pz();
	      pz_mu3 = muon3.Pz();
	      pz_mu4 = muon4.Pz();

	      E_mu1 = sqrt((muon1.P() * muon1.P()) + sqm1);
	      E_mu2 = sqrt((muon2.P() * muon2.P()) + sqm1);
	      E_mu3 = sqrt((muon3.P() * muon3.P()) + sqm1);
	      E_mu4 = sqrt((muon4.P() * muon4.P()) + sqm1);
      //cout << run << ", " << event << "," << aa << ", " << mass4mu << ", " << pxZa  << ", "<< pyZa << ", "<<  pzZa << ", " << pxZb  << ", "<< pyZb  << ", "<<  pzZb <<  "," << muon1.Pt() << "," << muon2.Pt() << "," << muon3.Pt() << "," << muon4.Pt() << endl;
	   // cout << run << ", " << event << "," << aa << ", " << mass4mu << endl;
          if (mass4mu > 70.) {

	        h_m1_m4mu->Fill(mass4mu);
	        h_m2_m4mu->Fill(mass4mu);
	        h_m3_m4mu->Fill(mass4mu);
	        h_m4_m4mu->Fill(mass4mu);
         // cout << run << ", "<< event << ", "<< mass4mu << "," << mZa << "," <<mZb << "," << pt_4mu <<  "," << muon1.Pt() << "," << muon2.Pt() << "," << muon3.Pt() << "," << muon4.Pt() << endl;
	        

	        for (int i = 0; i < nGoodMuons; i++) {

	          h_relPFIso_mu_after->Fill(Muon_pfRelIso04_all[goodMuons.at(i).first]);
	          h_pt_after->Fill(goodMuons.at(i).second);
	          h_eta_after->Fill(Muon_eta[goodMuons.at(i).first]);

                }
	      }
	    } // end of mZb
          } // end of mZa
        } // end of ptZadaug
      } // total charge
    } // nGoodMuons >= 4


    ////////////////////////// ZZ to 4e /////////////////////////////

    if (nGoodElectrons >= 4) {

      int e1 = goodElectrons.at(0).first;       
      int e2 = goodElectrons.at(1).first;      
      int e3 = goodElectrons.at(2).first;
      int e4 = goodElectrons.at(3).first;
     
      if (Electron_charge[e1] + Electron_charge[e2] + Electron_charge[e3] + Electron_charge[e4]== 0) { 

        ele1.SetPtEtaPhiM(Electron_pt[e1], Electron_eta[e1], Electron_phi[e1], Electron_mass[e1]);
        ele2.SetPtEtaPhiM(Electron_pt[e2], Electron_eta[e2], Electron_phi[e2], Electron_mass[e2]);
        ele3.SetPtEtaPhiM(Electron_pt[e3], Electron_eta[e3], Electron_phi[e3], Electron_mass[e3]);
        ele4.SetPtEtaPhiM(Electron_pt[e4], Electron_eta[e4], Electron_phi[e4], Electron_mass[e4]);

	// First combination: Combine ele 1234
  bool combo1234=false,combo1324=false,combo1423=false;   // again we must verify which combo is actually being formed
    
	if (Electron_charge[e1] + Electron_charge[e2] == 0) {
     
    combo1234=true; 
	  eZ12 = (sqrt(ele1.P() * ele1.P() + sqm1)) +
		 (sqrt(ele2.P() * ele2.P() + sqm1));
	
	  pxZ12 = ele1.Px() + ele2.Px();
	  pyZ12 = ele1.Py() + ele2.Py();
	  pzZ12 = ele1.Pz() + ele2.Pz();

	  if (Electron_charge[e3] + Electron_charge[e4] == 0) {
            eZ34 = (sqrt(ele3.P() * ele3.P() + sqm1)) +
		   (sqrt(ele4.P() * ele4.P() + sqm1));

            pxZ34 = ele3.Px() + ele4.Px();
	    pyZ34 = ele3.Py() + ele4.Py();
	    pzZ34 = ele3.Pz() + ele4.Pz();

	    // Calculate p4
	    pZ12 = sqrt((pxZ12 * pxZ12) + (pyZ12 * pyZ12) + (pzZ12 * pzZ12));
	    pZ34 = sqrt((pxZ34 * pxZ34) + (pyZ34 * pyZ34) + (pzZ34 * pzZ34));

	    pTZ12 = sqrt((pxZ12 * pxZ12) + (pyZ12 * pyZ12));
	    pTZ34 = sqrt((pxZ34 * pxZ34) + (pyZ34 * pyZ34));

	    mZ12 = sqrt((eZ12 * eZ12) - (pZ12 * pZ12));
	    mZ34 = sqrt((eZ34 * eZ34) - (pZ34 * pZ34));

	    if (mZ12 > 0.) h_mZ12_4e->Fill(mZ12);
	    if (mZ34 > 0.) h_mZ34_4e->Fill(mZ34);
              
          }
        }

	dZ12 = abs( mZ12 - mZ );
	dZ34 = abs( mZ34 - mZ );

	// take the smallest difference between mass
	// to use for 4electron combination
	dZc1 = (dZ12 < dZ34) ? dZ12 : dZ34; 
 
	// Second combination: Combine electron 1324
	if (Electron_charge[e1] + Electron_charge[e3] == 0) {
    
    combo1324=true;
	  eZ13 = (sqrt(ele1.P() * ele1.P() + sqm1)) +
		 (sqrt(ele3.P() * ele3.P() + sqm1));

	  pxZ13 = ele1.Px() + ele3.Px();
	  pyZ13 = ele1.Py() + ele3.Py();
	  pzZ13 = ele1.Pz() + ele3.Pz();

	  if (Electron_charge[e2] + Electron_charge[e4] == 0) {

            eZ24 = (sqrt(ele2.P() * ele2.P() + sqm1)) +
		   (sqrt(ele4.P() * ele4.P() + sqm1));

            pxZ24 = ele2.Px() + ele4.Px();
	    pyZ24 = ele2.Py() + ele4.Py();
	    pzZ24 = ele2.Pz() + ele4.Pz();

	    // Calculate p4
	    pZ13 = sqrt((pxZ13 * pxZ13) + (pyZ13 * pyZ13) + (pzZ13 * pzZ13));
	    pZ24 = sqrt((pxZ24 * pxZ24) + (pyZ24 * pyZ24) + (pzZ24 * pzZ24));

	    pTZ13 = sqrt((pxZ13 * pxZ13) + (pyZ13 * pyZ13));
	    pTZ24 = sqrt((pxZ24 * pxZ24) + (pyZ24 * pyZ24));

	    mZ13 = sqrt((eZ13 * eZ13) - (pZ13 * pZ13));
	    mZ24 = sqrt((eZ24 * eZ24) - (pZ24 * pZ24));

	    if (mZ13 > 0.) h_mZ13_4e->Fill(mZ13);
            if (mZ24 > 0.) h_mZ24_4e->Fill(mZ24);
	  }
	}

	dZ13 = abs( mZ13 - mZ );
	dZ24 = abs( mZ24 - mZ );

	dZc2 = (dZ13 < dZ24) ? dZ13 : dZ24; 

	// Third combination: Combine electron 1423
  	if (Electron_charge[e1] + Electron_charge[e4] == 0) {
    
    combo1423=true;
	  eZ14 = (sqrt(ele1.P() * ele1.P() + sqm1)) +
		 (sqrt(ele4.P() * ele4.P() + sqm1));

	  pxZ14 = ele1.Px() + ele4.Px();
	  pyZ14 = ele1.Py() + ele4.Py();
	  pzZ14 = ele1.Pz() + ele4.Pz();

	  if (Electron_charge[e2] + Electron_charge[e3] == 0) {

            eZ23 = (sqrt(ele2.P() * ele2.P() + sqm1)) +
		   (sqrt(ele3.P() * ele3.P() + sqm1));
       
	    pxZ23 = ele2.Px() + ele3.Px();
	    pyZ23 = ele2.Py() + ele3.Py();
	    pzZ23 = ele2.Pz() + ele3.Pz();

	    // Calculate p4
	    pZ14 = sqrt((pxZ14 * pxZ14) + (pyZ14 * pyZ14) + (pzZ14 * pzZ14));
	    pZ23 = sqrt((pxZ23 * pxZ23) + (pyZ23 * pyZ23) + (pzZ23 * pzZ23));

	    pTZ14 = sqrt((pxZ14 * pxZ14) + (pyZ14 * pyZ14));
	    pTZ23 = sqrt((pxZ23 * pxZ23) + (pyZ23 * pyZ23));

	    mZ14 = sqrt((eZ14 * eZ14) - (pZ14 * pZ14));
	    mZ23 = sqrt((eZ23 * eZ23) - (pZ23 * pZ23));
  
            if (mZ14 > 0.) h_mZ14_4e->Fill(mZ14);
	    if (mZ23 > 0.) h_mZ23_4e->Fill(mZ23);

          }
        }

	dZ14 = abs( mZ14 - mZ );
	dZ23 = abs( mZ23 - mZ );

	dZc3 = (dZ14 < dZ23) ? dZ14 : dZ23;

	bool ptZadaug = false;
  
  if (!combo1324) {    // This means combo 13-24 is not possible, so we must only check between 12-34 and 14-23
	
  if (dZc1 < dZc3) {   // the first possibility is that combo 12-34 is closer in mass

	  if (dZ12 < dZ34) {  // Now we examine the sub-possibility, where Z12 is closer in mass
	    eZa  = eZ12;     
	    pxZa = pxZ12;
	    pyZa = pyZ12;
	    pzZa = pzZ12;
	    pTZa = pTZ12;
	    mZa  = mZ12;

	    if (ele1.Pt() > 20. and ele2.Pt() > 10.)
	    ptZadaug = true;

	    eZb  = eZ34;
	    pxZb = pxZ34;
	    pyZb = pyZ34;
	    pzZb = pzZ34;
	    pTZb = pTZ34;
	    mZb  = mZ34;
 
	  }
	      
          else { // The sub-possibility where Z34 is closer in mass

	    eZa  = eZ34;
	    pxZa = pxZ34;
	    pyZa = pyZ34;
	    pzZa = pzZ34;
	    pTZa = pTZ34;
	    mZa  = mZ34;

	    if (ele3.Pt() > 20. and ele4.Pt() > 10.)
	    ptZadaug = true;

	    eZb  = eZ12;
	    pxZb = pxZ12;
	    pyZb = pyZ12;
	    pzZb = pzZ12;
	    pTZb = pTZ12;
	    mZb  = mZ12;

          }
        }
     else {  // we move to the case where combo 14-23 is closer in mass
      
      if (dZ14 < dZ23) {

	    eZa  = eZ14;
	    pxZa = pxZ14;
	    pyZa = pyZ14;
	    pzZa = pzZ14;
	    pTZa = pTZ14;
	    mZa  = mZ14;

	    if (ele1.Pt() > 20. and ele4.Pt() > 10.)
	    ptZadaug = true;

	    eZb  = eZ23;
	    pxZb = pxZ23;
	    pyZb = pyZ23;
	    pzZb = pzZ23;
	    pTZb = pTZ23;
	    mZb  = mZ23;

	  }

	  else {

	    eZa  = eZ23;
	    pxZa = pxZ23;
	    pyZa = pyZ23;
	    pzZa = pzZ23;
	    pTZa = pTZ23;
	    mZa  = mZ23;

	    if (ele2.Pt() > 20. and ele3.Pt() > 10.)
	    ptZadaug = true;

	    eZb  = eZ14;
	    pxZb = pxZ14;
	    pyZb = pyZ14;
	    pzZb = pzZ14;
	    pTZb = pTZ14;
	    mZb  = mZ14;

	  }

     }   
  } // end of the check for the case where only combo 12-34 and combo 14-23 are possible.      
  
  else if (!combo1234) {  // that means combo 12-34 is not possible, we examine the remaining two

	  if (dZc2 < dZc3) {  // This is the case where combo 13-24 is closer in mass

	  if (dZ13 < dZ24) {

	    eZa  = eZ13;
	    pxZa = pxZ13;
	    pyZa = pyZ13;
	    pzZa = pzZ13;
	    pTZa = pTZ13;
	    mZa  = mZ13;

	    if (ele1.Pt() > 20. and ele3.Pt() > 10.)
	    ptZadaug = true;

	    eZb  = eZ24;
	    pxZb = pxZ24;
	    pyZb = pyZ24;
	    pzZb = pzZ24;
	    pTZb = pTZ24;
	    mZb  = mZ24;

	  }
	      
          else {

	    eZa  = eZ24;
	    pxZa = pxZ24;
	    pyZa = pyZ24;
	    pzZa = pzZ24;
	    pTZa = pTZ24;
	    mZa  = mZ24;

	    if (ele2.Pt() > 20. and ele4.Pt() > 10.)
	    ptZadaug = true;

	    eZb  = eZ13;
	    pxZb = pxZ13;
	    pyZb = pyZ13;
	    pzZb = pzZ13;
	    pTZb = pTZ13;
	    mZb  = mZ13;

	  }
	}
  else {
    if (dZ14 < dZ23) {

	    eZa  = eZ14;
	    pxZa = pxZ14;
	    pyZa = pyZ14;
	    pzZa = pzZ14;
	    pTZa = pTZ14;
	    mZa  = mZ14;

	    if (ele1.Pt() > 20. and ele4.Pt() > 10.)
	    ptZadaug = true;

	    eZb  = eZ23;
	    pxZb = pxZ23;
	    pyZb = pyZ23;
	    pzZb = pzZ23;
	    pTZb = pTZ23;
	    mZb  = mZ23;

	  }

	  else {

	    eZa  = eZ23;
	    pxZa = pxZ23;
	    pyZa = pyZ23;
	    pzZa = pzZ23;
	    pTZa = pTZ23;
	    mZa  = mZ23;

	    if (ele2.Pt() > 20. and ele3.Pt() > 10.)
	    ptZadaug = true;

	    eZb  = eZ14;
	    pxZb = pxZ14;
	    pyZb = pyZ14;
	    pzZb = pzZ14;
	    pTZb = pTZ14;
	    mZb  = mZ14;

	  }
  }
  }
  else {    // This only leaves the case where combo 14-23 is not possible, we examine 12-34 and 13-24

	if (dZc1 < dZc2) {  // first possibility, combo 12-34 is closer in mass

	  if (dZ12 < dZ34) {  // Now we examine the sub-possibility, where Z12 is closer in mass
	    eZa  = eZ12;     
	    pxZa = pxZ12;
	    pyZa = pyZ12;
	    pzZa = pzZ12;
	    pTZa = pTZ12;
	    mZa  = mZ12;

	    if (ele1.Pt() > 20. and ele2.Pt() > 10.)
	    ptZadaug = true;

	    eZb  = eZ34;
	    pxZb = pxZ34;
	    pyZb = pyZ34;
	    pzZb = pzZ34;
	    pTZb = pTZ34;
	    mZb  = mZ34;
 
	  }
	      
          else { // The sub-possibility where Z34 is closer in mass

	    eZa  = eZ34;
	    pxZa = pxZ34;
	    pyZa = pyZ34;
	    pzZa = pzZ34;
	    pTZa = pTZ34;
	    mZa  = mZ34;

	    if (ele3.Pt() > 20. and ele4.Pt() > 10.)
	    ptZadaug = true;

	    eZb  = eZ12;
	    pxZb = pxZ12;
	    pyZb = pyZ12;
	    pzZb = pzZ12;
	    pTZb = pTZ12;
	    mZb  = mZ12;

          }
	}
  else {  // Now examine the case where combo 13-24 is closer in mass
   if (dZ13 < dZ24) {

	    eZa  = eZ13;
	    pxZa = pxZ13;
	    pyZa = pyZ13;
	    pzZa = pzZ13;
	    pTZa = pTZ13;
	    mZa  = mZ13;

	    if (ele1.Pt() > 20. and ele3.Pt() > 10.)
	    ptZadaug = true;

	    eZb  = eZ24;
	    pxZb = pxZ24;
	    pyZb = pyZ24;
	    pzZb = pzZ24;
	    pTZb = pTZ24;
	    mZb  = mZ24;

	  }
	      
          else {

	    eZa  = eZ24;
	    pxZa = pxZ24;
	    pyZa = pyZ24;
	    pzZa = pzZ24;
	    pTZa = pTZ24;
	    mZa  = mZ24;

	    if (ele2.Pt() > 20. and ele4.Pt() > 10.)
	    ptZadaug = true;

	    eZb  = eZ13;
	    pxZb = pxZ13;
	    pyZb = pyZ13;
	    pzZb = pzZ13;
	    pTZb = pTZ13;
	    mZb  = mZ13;

	  }
  }
  }

	if (ptZadaug) {

	  if (mZa > 40. && mZa < 120.) {

	    if (mZb > 12. && mZb < 120.) {

	      h_mZa_4e->Fill(mZa);
	      h_mZb_4e->Fill(mZb);

	      // 4 vector
	      p4Za.SetPxPyPzE(pxZa, pyZa, pzZa, eZa);
	      p4Zb.SetPxPyPzE(pxZb, pyZb, pzZb, eZb);

	      p4H = p4Za + p4Zb;

	      mass4e = p4H.M();
	      pt_4e = p4H.Pt();
	      eta_4e = p4H.Eta();
	      phi_4e = p4H.Phi();

	      px4e = p4H.Px();
	      py4e = p4H.Py();
	      pz4e = p4H.Pz();
	      E4e = p4H.E();

	      pt_e1 = ele1.Pt();
	      pt_e2 = ele2.Pt();
	      pt_e3 = ele3.Pt();
	      pt_e4 = ele4.Pt();

	      eta_e1 = ele1.Eta();
	      eta_e2 = ele2.Eta();
	      eta_e3 = ele3.Eta();
	      eta_e4 = ele4.Eta();

	      phi_e1 = ele1.Phi();
	      phi_e2 = ele2.Phi();
	      phi_e3 = ele3.Phi();
	      phi_e4 = ele4.Phi();

	      cas_e1 = Electron_charge[e1];
	      cas_e2 = Electron_charge[e2];
	      cas_e3 = Electron_charge[e3];
	      cas_e4 = Electron_charge[e4];

	      px_e1 = ele1.Px();
	      px_e2 = ele2.Px();
	      px_e3 = ele3.Px();
	      px_e4 = ele4.Px();

	      py_e1 = ele1.Py();
	      py_e2 = ele2.Py();
	      py_e3 = ele3.Py();
	      py_e4 = ele4.Py();
		
	      pz_e1 = ele1.Pz();
	      pz_e2 = ele2.Pz();
	      pz_e3 = ele3.Pz();
	      pz_e4 = ele4.Pz();

	      E_e1 = sqrt((ele1.P() * ele1.P()) + sqm1);
	      E_e2 = sqrt((ele2.P() * ele2.P()) + sqm1);
	      E_e3 = sqrt((ele3.P() * ele3.P()) + sqm1);
	      E_e4 = sqrt((ele4.P() * ele4.P()) + sqm1);

	      if (mass4e > 70.) {

	        h_m1_m4e->Fill(mass4e);
	        h_m2_m4e->Fill(mass4e);
	        h_m3_m4e->Fill(mass4e);
	        h_m4_m4e->Fill(mass4e);

	        for (int i = 0; i < nGoodElectrons; i++) {

	          h_relPFIso_e_after->Fill(Electron_pfRelIso03_all[goodElectrons.at(i).first]);
	          h_pt_after->Fill(goodElectrons.at(i).second);
	          h_eta_after->Fill(Electron_eta[goodElectrons.at(i).first]);

                }
	      }
	    } // end of mZb
          } // end of mZa
        } // end of ptZadaug
      } // total charge 
    } // nGoodElectrons >= 4


    ///////////////////////// ZZ to 2mu2e ///////////////////////////

    if (nGoodMuons >= 2 && nGoodElectrons >= 2) {

      int m1 = goodMuons.at(0).first;
      int m2 = goodMuons.at(1).first;
      int e1 = goodElectrons.at(0).first;
      int e2 = goodElectrons.at(1).first;
     
      if (Muon_charge[m1] + Muon_charge[m2] + Electron_charge[e1] + Electron_charge[e2]== 0) { 

        muon1.SetPtEtaPhiM(Muon_pt[m1], Muon_eta[m1], Muon_phi[m1], Muon_mass[m1]);
        muon2.SetPtEtaPhiM(Muon_pt[m2], Muon_eta[m2], Muon_phi[m2], Muon_mass[m2]);
        ele1.SetPtEtaPhiM(Electron_pt[e1], Electron_eta[e1], Electron_phi[e1], Electron_mass[e1]);
        ele2.SetPtEtaPhiM(Electron_pt[e2], Electron_eta[e2], Electron_phi[e2], Electron_mass[e2]);
	  
        // For case 2mu2e, there is only 1 combination
	if (Muon_charge[m1] + Muon_charge[m2] == 0) {

	  eZ12 = (sqrt(muon1.P() * muon1.P() + sqm1)) + (sqrt(muon2.P() * muon2.P() + sqm1));
	
	  pxZ12 = muon1.Px() + muon2.Px();
	  pyZ12 = muon1.Py() + muon2.Py();
	  pzZ12 = muon1.Pz() + muon2.Pz();

	  if (Electron_charge[e1] + Electron_charge[e2] == 0) {
        
            eZ34 = (sqrt(ele1.P() * ele1.P() + sqme)) + (sqrt(ele2.P() * ele2.P() + sqme));

	    pxZ34 = ele1.Px() + ele2.Px();
	    pyZ34 = ele1.Py() + ele2.Py();
	    pzZ34 = ele1.Pz() + ele2.Pz();

            // Calculate the momentum and invariant mass
            // Muons = 1 2, electrons = 3 4
            pZ12 = sqrt((pxZ12 * pxZ12) + (pyZ12 * pyZ12) + (pzZ12 * pzZ12));
	    pZ34 = sqrt((pxZ34 * pxZ34) + (pyZ34 * pyZ34) + (pzZ34 * pzZ34));

	    pTZ12 = sqrt((pxZ12 * pxZ12) + (pyZ12 * pyZ12));
	    pTZ34 = sqrt((pxZ34 * pxZ34) + (pyZ34 * pyZ34));

	    mZ12 = sqrt((eZ12 * eZ12) - (pZ12 * pZ12));
	    mZ34 = sqrt((eZ34 * eZ34) - (pZ34 * pZ34));

	    if (mZ12 > 0.) h_mZmu_2mu2e->Fill(mZ12);
	    if (mZ34 > 0.) h_mZe_2mu2e->Fill(mZ34);

	  }
	}

	dZ12 = abs(mZ12 - mZ); // mu
	dZ34 = abs(mZ34 - mZ); // e

	bool ptZadaug = false;

	if (dZ12 < dZ34) {

	  eZa  = eZ12;
	  pxZa = pxZ12;
	  pyZa = pyZ12;
	  pzZa = pzZ12;
	  pTZa = pTZ12;
	  mZa  = mZ12;

	  if (muon1.Pt() > 20. && muon2.Pt() > 10.) {ptZadaug = true;}

	  eZb  = eZ34;
	  pxZb = pxZ34;
	  pyZb = pyZ34;
	  pzZb = pzZ34;
	  pTZb = pTZ34;
	  mZb  = mZ34;

	}
	  
        else {

	  eZa  = eZ34;
	  pxZa = pxZ34;
	  pyZa = pyZ34;
	  pzZa = pzZ34;
	  pTZa = pTZ34;
	  mZa  = mZ34;

	  if (ele1.Pt() > 20. && ele2.Pt() > 10.) {ptZadaug = true;}

	  eZb  = eZ12;
	  pxZb = pxZ12;
	  pyZb = pyZ12;
	  pzZb = pzZ12;
	  pTZb = pTZ12;
	  mZb  = mZ12;

	}

	if (ptZadaug) {
	  if (mZa > 40. && mZa < 120.) {
	    if (mZb > 12. && mZb < 120.) {
              h_mZa_2mu2e->Fill(mZa);
	      h_mZb_2mu2e->Fill(mZb);

	      // Now combine these 2 muons and 2 electrons
 
	      // Calculate 4 lepton: 2muons 2electrons

	      p4Za.SetPxPyPzE(pxZa, pyZa, pzZa, eZa);
	      p4Zb.SetPxPyPzE(pxZb, pyZb, pzZb, eZb);

	      p4H = p4Za + p4Zb;

	      mass2mu2e = p4H.M();
	      pt_2mu2e = p4H.Pt();
	      eta_2mu2e = p4H.Eta();
	      phi_2mu2e = p4H.Phi();
		
	      px2mu2e = p4H.Px();
	      py2mu2e = p4H.Py();
	      pz2mu2e = p4H.Pz();
	      E2mu2e = p4H.E();

	      pt_2mu1 = muon1.Pt();
	      pt_2mu2 = muon2.Pt();
	      pt_2e1 = ele1.Pt();
	      pt_2e2 = ele2.Pt();

	      eta_2mu1 = muon1.Eta();
	      eta_2mu2 = muon2.Eta();
	      eta_2e1 = ele1.Eta();
	      eta_2e2 = ele2.Eta();

	      phi_2mu1 = muon1.Phi();
	      phi_2mu2 = muon2.Phi();
	      phi_2e1 = ele1.Phi();
	      phi_2e2 = ele2.Phi();

	      cas_2mu1 = Muon_charge[m1];
	      cas_2mu2 = Muon_charge[m2];
	      cas_2e1 = Electron_charge[e1];
	      cas_2e2 = Electron_charge[e2];

	      px_2mu1 = muon1.Px();
	      px_2mu2 = muon2.Px();
	      px_2e1 = ele1.Px();
	      px_2e2 = ele2.Px();

	      py_2mu1 = muon1.Py();
	      py_2mu2 = muon2.Py();
	      py_2e1 = ele1.Py();
	      py_2e2 = ele2.Py();
		
	      pz_2mu1 = muon1.Pz();
	      pz_2mu2 = muon2.Pz();
	      pz_2e1 = ele1.Pz();
	      pz_2e2 = ele2.Pz();

	      E_2mu1 = sqrt((muon1.P() * muon1.P()) + sqm1);
	      E_2mu2 = sqrt((muon2.P() * muon2.P()) + sqm1);
	      E_2e1 = ele1.E();
	      E_2e2 = ele2.E();

	      if (mass2mu2e > 70.) {

		h_m1_m2mu2e->Fill(mass2mu2e);
		h_m2_m2mu2e->Fill(mass2mu2e);
		h_m3_m2mu2e->Fill(mass2mu2e);
		h_m4_m2mu2e->Fill(mass2mu2e);

		for (int mm = 0; mm < nGoodMuons; mm++) {

		  h_relPFIso_2mu_after->Fill(Muon_pfRelIso04_all[goodMuons.at(mm).first]);
		  h_pt_after_2mu2e->Fill(Muon_pt[goodMuons.at(mm).first]);
		  h_eta_after_2mu2e->Fill(Muon_eta[goodMuons.at(mm).first]);

		}

		for (int ee = 0; ee < nGoodElectrons; ee++) {
	      
		  h_relPFIso_2e_after->Fill(Electron_pfRelIso03_all[goodElectrons.at(ee).first]);
		  h_pt_e_after_2mu2e->Fill(Electron_pt[goodElectrons.at(ee).first]);
		  h_eta_e_after_2mu2e->Fill(Electron_eta[goodElectrons.at(ee).first]);

		}	    
	  } 
	} // mZb
  } // mZa
} // ptZadaug
} // end of total charge
} // nGoodMuons >= 2 && nGoodElectrons >= 2



  } // nevent
  

  
  // out file
  TFile fout((dir_out+"/"+name_out_file+".root").c_str(),"RECREATE");
  
  // ---------------------------------------------------------------------

  TDirectory *dir = fout.mkdir("Events");
  dir->cd(); 

  // Muon selection
  h_p_reco->Write();
  h_pt_reco_b4->Write();
  h_eta_reco_b4->Write();
  h_phi_reco->Write();
  h_relPFIso_mu->Write();
  h_dxy_mu->Write();
  h_SIP3d_mu_b4->Write();
  h_ngmu->Write(); 

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

  // Z to 2mu
  h_mZ_2mu->Write();
  h_pt_after_Zto2mu->Write();
  h_eta_after_Zto2mu->Write();

  // Z to 2e
  h_mZ_2e->Write();
  h_pt_e_after_Zto2e->Write();
  h_eta_e_after_Zto2e->Write();

  // ZZ to 4mu

  h_mZ12_4mu->Write();
  h_mZ34_4mu->Write();
  h_mZ13_4mu->Write();
  h_mZ24_4mu->Write();
  h_mZ14_4mu->Write();
  h_mZ23_4mu->Write();
  h_mZa_4mu->Write();
  h_mZb_4mu->Write();
  h_m1_m4mu->Write();
  h_m2_m4mu->Write();
  h_m3_m4mu->Write();
  h_m4_m4mu->Write(); 

  // ZZ to 4e

  h_mZ12_4e->Write();
  h_mZ34_4e->Write();
  h_mZ13_4e->Write();
  h_mZ24_4e->Write();
  h_mZ14_4e->Write();
  h_mZ23_4e->Write();
  h_mZa_4e->Write();
  h_mZb_4e->Write();
  h_m1_m4e->Write();
  h_m2_m4e->Write();
  h_m3_m4e->Write();
  h_m4_m4e->Write();

  // ZZ to 2mu 2e

  h_mZmu_2mu2e->Write();
  h_mZe_2mu2e->Write();
  h_mZa_2mu2e->Write();
  h_mZb_2mu2e->Write();
  h_m1_m2mu2e->Write();
  h_m2_m2mu2e->Write();
  h_m3_m2mu2e->Write();
  h_m4_m2mu2e->Write();
  h_relPFIso_2mu_after->Write();
  h_pt_after_2mu2e->Write();
  h_eta_after_2mu2e->Write();
  h_relPFIso_2e_after->Write();
  h_pt_e_after_2mu2e->Write();
  h_eta_e_after_2mu2e->Write();

  fout.Close();
  gROOT->ProcessLine(".q");

} // end 