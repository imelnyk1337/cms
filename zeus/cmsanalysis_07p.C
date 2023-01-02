#include "Riostream.h"
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TChain.h"
#include "TH1F.h"
#include "TH2D.h"
#include "TLorentzVector.h"
#include "TLatex.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TStopwatch.h"
#include "TStyle.h"
#include "TLine.h"
#include "TMath.h"
#include <string>

using namespace std;

const int maxlepton_n = 50; // 128 is the maximum number but it's only one event -> putting to 50 makes it running faster
const int maxjet_n = 50; 
const int pdgmuon = 13;
const int pdgelec = 11;
const int pdgtau = 15;
const double muonmass = 0.105658; // muon mass in GeV
const float PI = 3.14159265359; // just PI ;)

//definition
void cmsanalysisv8_path (std::string fin_path, int chain, std::string legend_title, int mcdata, std::string png_dir); //mcdata is 0 for data and 1 for mc, chain is 1 if one needs to use multiple files, 0 for one file only

//main
void cmsanalysis_07p () {

    TStopwatch tall;
    tall.Start();

    cout<<"-------------------------------"<<endl;
    cout<<" CMS nanoAODv8 - ANALYSIS "<<endl;
    cout<<"-------------------------------"<<endl;

    cout<<"-------------------------------"<<endl;
    cout<<" ATLAS DATA - 2012 "<<endl;
    cout<<"-------------------------------"<<endl;

    cmsanalysisv8_path("/nfs/dust/zeus/group/schwenzr/ZEUStoCMSfiles_test/zeustocms_07p.root", 0," ZEUS DATA", 0,"zeusdata");
    

/*
    cout<<"-------------------------------"<<endl;
    cout<<" ATLAS MC - 2012 "<<endl;
    cout<<"-------------------------------"<<endl;

    cmsanalysisv8_path("/nfs/dust/cms/user/lolivi/ATLAS_CMS_MAP/cms_atlas2012_mc_new.root", 0," ATLAS MC - 2012 ", 1,"atlas2012mc");

    cout<<"-------------------------------"<<endl;
    cout<<" ATLAS DATA - 2016 "<<endl;
    cout<<"-------------------------------"<<endl;

    //cmsanalysisv8_path("/nfs/dust/cms/user/arbal/ATLAS2016Data/", 1," ATLAS DATA - 2016 ", 0,"atlas2016data");

    cout<<"-------------------------------"<<endl;
    cout<<" CMS DATA - 2011 - OLD "<<endl;
    cout<<"-------------------------------"<<endl;
    
    //cmsanalysisv8_path("/nfs/dust/cms/user/geiser/nanoAOD/2011/Data/Data11_DoubleMu.root", 0," CMS DOUBLE MUON DATA - 2011 - OLD ",0,"cms2011data_old");

    cout<<"-------------------------------"<<endl;
    cout<<" CMS MC - 2011 - OLD "<<endl;
    cout<<"-------------------------------"<<endl;
    
    //cmsanalysisv8_path("/nfs/dust/cms/user/geiser/nanoAOD/2011/MC/MC11_DY50.root", 0," CMS DOUBLE MUON MC - 2011 - OLD ",1,"cms2011mc_old");

    cout<<"-------------------------------"<<endl;
    cout<<" CMS DOUBLE MUON DATA RUN A - 2011 - NEW "<<endl;
    cout<<"-------------------------------"<<endl;

    //cmsanalysisv8_path("/nfs/dust/cms/user/yangq2/NanoAODplus_v1/output/data/2011/DoubleMu/CMS_Run2011A_DoubleMu_AOD_12Oct2013-v1.root", 0," CMS DOUBLE MUON DATA RUN A - 2011 - NEW ",0,"cms2011runA_new");

    cout<<"-------------------------------"<<endl;
    cout<<" CMS DOUBLE MUON DATA RUN B - 2011 - NEW "<<endl;
    cout<<"-------------------------------"<<endl;

    //cmsanalysisv8_path("/nfs/dust/cms/user/yangq2/NanoAODplus_v1/output/data/2011/DoubleMu/CMS_Run2011B_DoubleMu_AOD_12Oct2013-v1.root", 0," CMS DOUBLE MUON DATA RUN B - 2011 - NEW ",0,"cms2011runB_new");

    cout<<"-------------------------------"<<endl;
    cout<<" CMS MUONIA DATA RUN A - 2010 - NEW "<<endl;
    cout<<"-------------------------------"<<endl;

    cmsanalysisv8_path("/nfs/dust/cms/user/yangq2/NanoAODplus_v1/output/data/2010/MuOnia/CMS_Run2010A_MuOnia_AOD_Apr21ReReco-v1.root", 0," CMS MUONIA DATA RUN A - 2010 - NEW ",0,"cms2010muonia_runA");

    cout<<"-------------------------------"<<endl;
    cout<<" CMS MUONIA DATA RUN B - 2010 - NEW "<<endl;
    cout<<"-------------------------------"<<endl;

    cmsanalysisv8_path("/nfs/dust/cms/user/yangq2/NanoAODplus_v1/output/data/2010/MuOnia/CMS_Run2010B_MuOnia_AOD_Apr21ReReco-v1.root", 0," CMS MUONIA DATA RUN B - 2010 - NEW ",0,"cms2010muonia_runB");


    cout<<"-------------------------------"<<endl;
    cout<<" CMS DOUBLE MUON DATA - 2016B Ver1 "<<endl;
    cout<<"-------------------------------"<<endl;

    //cmsanalysisv8_path("/nfs/dust/cms/user/geiser/nanoAODv8/DoubleMuon/2016Bver1/", 1," CMS DOUBLE MUON DATA - 2016B Ver1 ",0,"cms2016B1data");

    cout<<"-------------------------------"<<endl;
    cout<<" CMS DOUBLE MUON DATA - 2016B Ver2 "<<endl;
    cout<<"-------------------------------"<<endl;

    //cmsanalysisv8_path("/nfs/dust/cms/user/geiser/nanoAODv8/DoubleMuon/2016Bver2/", 1," CMS DOUBLE MUON DATA - 2016B Ver2 ",0,"cms2016B2data");

    cout<<"-------------------------------"<<endl;
    cout<<" CMS DOUBLE MUON DATA - 2016C "<<endl;
    cout<<"-------------------------------"<<endl;

    //cmsanalysisv8_path("/nfs/dust/cms/user/geiser/nanoAODv8/DoubleMuon/2016C/", 1," CMS DOUBLE MUON DATA - 2016C ",0,"cms2016Cdata");

    cout<<"-------------------------------"<<endl;
    cout<<" CMS DOUBLE MUON DATA - 2016D "<<endl;
    cout<<"-------------------------------"<<endl;

    //cmsanalysisv8_path("/nfs/dust/cms/user/geiser/nanoAODv8/DoubleMuon/2016D/", 1," CMS DOUBLE MUON DATA - 2016D ",0,"cms2016Ddata");

    cout<<"-------------------------------"<<endl;
    cout<<" CMS DOUBLE MUON DATA - 2016E "<<endl;
    cout<<"-------------------------------"<<endl;

    //cmsanalysisv8_path("/nfs/dust/cms/user/geiser/nanoAODv8/DoubleMuon/2016E/", 1," CMS DOUBLE MUON DATA - 2016E ",0,"cms2016Edata");

    cout<<"-------------------------------"<<endl;
    cout<<" CMS DOUBLE MUON DATA - 2016F "<<endl;
    cout<<"-------------------------------"<<endl;

    //cmsanalysisv8_path("/nfs/dust/cms/user/geiser/nanoAODv8/DoubleMuon/2016F/", 1," CMS DOUBLE MUON DATA - 2016F ",0,"cms2016Fdata");

    cout<<"-------------------------------"<<endl;
    cout<<" CMS DOUBLE MUON DATA - 2016FHIPM "<<endl;
    cout<<"-------------------------------"<<endl;

    //cmsanalysisv8_path("/nfs/dust/cms/user/geiser/nanoAODv8/DoubleMuon/2016FHIPM/", 1," CMS DOUBLE MUON DATA - 2016FHIPM ",0,"cms2016FHIPMdata");

    cout<<"-------------------------------"<<endl;
    cout<<" CMS DOUBLE MUON DATA - 2016G "<<endl;
    cout<<"-------------------------------"<<endl;

    //cmsanalysisv8_path("/nfs/dust/cms/user/geiser/nanoAODv8/DoubleMuon/2016G/", 1," CMS DOUBLE MUON DATA - 2016G ",0,"cms2016Gdata");

    cout<<"-------------------------------"<<endl;
    cout<<" CMS DOUBLE MUON DATA - 2016H "<<endl;
    cout<<"-------------------------------"<<endl;

    //cmsanalysisv8_path("/nfs/dust/cms/user/geiser/nanoAODv8/DoubleMuon/2016H/", 1," CMS DOUBLE MUON DATA - 2016H ",0,"cms2016Hdata");
*/

    tall.Stop(); // global timer
    cout << "- Global Timer: ";
    tall.Print();

    cout << "-------------------------------" << endl;
    cout << " END OF FILE " << endl;
    cout << "-------------------------------" << endl;

    return;

}

//implementation
void cmsanalysisv8_path (std::string fin_path, int chain, std::string legend_title, int mcdata, std::string png_dir) {

    //timer
    TStopwatch t;
    t.Start();

    ifstream check_in(fin_path.c_str()); //check if input exists
    if(!check_in) {
        cout<<fin_path<<" does not exist \n";
        return; 
    }

    //CMS ntuples - "nanoAODv8" [https://cms-nanoaod-integration.web.cern.ch/integration/master-106X/mc102X_doc.html] (Only the converted ones, needs to be expanded)
    Int_t PV_npvs, TrigObj_id;
    Int_t Muon_charge[maxlepton_n], Electron_charge[maxlepton_n], Tau_charge[maxlepton_n];
    static UInt_t nMuon, nElectron = 0, nTau = 0, run, luminosityBlock;
    ULong64_t event;
    Float_t PV_z, PV_chi2;
    Float_t Muon_pt[maxlepton_n], Electron_pt[maxlepton_n], Tau_pt[maxlepton_n], Muon_eta[maxlepton_n], Electron_eta[maxlepton_n], Tau_eta[maxlepton_n];
    Float_t Muon_phi[maxlepton_n], Electron_phi[maxlepton_n], Tau_phi[maxlepton_n], Muon_dz[maxlepton_n], Electron_dz[maxlepton_n], Tau_dz[maxlepton_n];
    Float_t Muon_dxy[maxlepton_n], Electron_dxy[maxlepton_n], Tau_dxy[maxlepton_n];
    Float_t Muon_dxyErr[maxlepton_n], Electron_dxyErr[maxlepton_n], Tau_dxyErr[maxlepton_n];
    Float_t Muon_ip3d[maxlepton_n], Electron_ip3d[maxlepton_n], Muon_sip3d[maxlepton_n], Electron_sip3d[maxlepton_n];
    Float_t Muon_pfRelIso03_all[maxlepton_n], Electron_pfRelIso03_all[maxlepton_n], Tau_pfRelIso03_all[maxlepton_n];
    Float_t Muon_pfRelIso03_chg[maxlepton_n], Electron_pfRelIso03_chg[maxlepton_n], Tau_pfRelIso03_chg[maxlepton_n];
    Float_t Muon_pfRelIso04_all[maxlepton_n];
    Bool_t Muon_isGlobal[maxlepton_n], Muon_isTracker[maxlepton_n], Muon_softId[maxlepton_n], Muon_isPFcand[maxlepton_n];

    //input addressing
    TChain *tin = new TChain("Events");
    cout << "- Reading files -> "<<fin_path<<endl;
    if (chain==0) tin->Add(fin_path.c_str());
    if (chain==1) tin->Add((fin_path + std::string("*.root")).c_str());
    int nentries = (int)tin->GetEntries();

    if (nentries==0) {
        cout<<"- \"Events\" Tree not found \n";
        tin = new TChain("Event");    
    }

    if (chain==0) tin->Add(fin_path.c_str());
    if (chain==1) tin->Add((fin_path + std::string("*.root")).c_str());
    nentries = (int)tin->GetEntries();

    if (nentries==0) {
        cout<<"- \"Events\" Tree not found \n";
        return;   
    }

    // could create problems
    tin->SetBranchStatus("Muon_cleanmask",0);
    tin->SetBranchStatus("Electron_cleanmask",0);
    tin->SetBranchStatus("Tau_cleanmask",0);

    // general
    tin->SetBranchAddress("event", &event);
    tin->SetBranchAddress("run", &run);
    tin->SetBranchAddress("luminosityBlock", &luminosityBlock);
    tin->SetBranchAddress("PV_npvs",&PV_npvs);
    tin->SetBranchAddress("PV_z",&PV_z);
    tin->SetBranchAddress("PV_chi2",&PV_chi2);
    tin->SetBranchAddress("TrigObj_id",&TrigObj_id);

    // muons
    tin->SetBranchAddress("nMuon", &nMuon);
    tin->SetBranchAddress("Muon_pt",Muon_pt);
    tin->SetBranchAddress("Muon_eta",Muon_eta);
    tin->SetBranchAddress("Muon_phi",Muon_phi);
    tin->SetBranchAddress("Muon_dz",Muon_dz);
    tin->SetBranchAddress("Muon_charge",Muon_charge);
    tin->SetBranchAddress("Muon_dxy",Muon_dxy);
    tin->SetBranchAddress("Muon_dxyErr",Muon_dxyErr);
    tin->SetBranchAddress("Muon_ip3d",Muon_ip3d);
    tin->SetBranchAddress("Muon_sip3d",Muon_sip3d);
    tin->SetBranchAddress("Muon_pfRelIso03_chg",Muon_pfRelIso03_chg);
    tin->SetBranchAddress("Muon_pfRelIso03_all",Muon_pfRelIso03_all);
    tin->SetBranchAddress("Muon_pfRelIso04_all",Muon_pfRelIso04_all);

    //muons branches only on CMS (setting ATLAS to true)
    TBranch *bglobalmu = tin->GetBranch("Muon_isGlobal");
    if (bglobalmu==NULL) for(int i=0;i<maxlepton_n;i++) Muon_isGlobal[i]=true;
    else tin->SetBranchAddress("Muon_isGlobal",Muon_isGlobal);

    TBranch *bsoftmu = tin->GetBranch("Muon_softId");
    if (bsoftmu==NULL) for(int i=0;i<maxlepton_n;i++) Muon_softId[i]=true;
    else tin->SetBranchAddress("Muon_softId",Muon_softId);

    TBranch *btrackmu = tin->GetBranch("Muon_isTracker");
    if (btrackmu==NULL) for(int i=0;i<maxlepton_n;i++) Muon_isTracker[i]=true;
    else tin->SetBranchAddress("Muon_isTracker",Muon_isTracker);

    TBranch *bpfcand = tin->GetBranch("Muon_isPFcand");
    if (bpfcand==NULL) for(int i=0;i<maxlepton_n;i++) Muon_isPFcand[i]=true;
    else tin->SetBranchAddress("Muon_isPFcand",Muon_isPFcand);

    //electrons
    tin->SetBranchAddress("nElectron", &nElectron); 
    tin->SetBranchAddress("Electron_pt",Electron_pt);
    tin->SetBranchAddress("Electron_eta",Electron_eta);
    tin->SetBranchAddress("Electron_phi",Electron_phi);
    tin->SetBranchAddress("Electron_dz",Electron_dz);
    tin->SetBranchAddress("Electron_charge",Electron_charge);
    tin->SetBranchAddress("Electron_dxy",Electron_dxy);
    tin->SetBranchAddress("Electron_dxyErr",Electron_dxyErr);
    tin->SetBranchAddress("Electron_ip3d",Electron_ip3d);
    tin->SetBranchAddress("Electron_sip3d",Electron_sip3d);
    tin->SetBranchAddress("Electron_pfRelIso03_chg",Electron_pfRelIso03_chg);
    tin->SetBranchAddress("Electron_pfRelIso03_all",Electron_pfRelIso03_all);

    //taus
    tin->SetBranchAddress("nTau", &nTau); 
    tin->SetBranchAddress("Tau_pt",Tau_pt);
    tin->SetBranchAddress("Tau_eta",Tau_eta);
    tin->SetBranchAddress("Tau_phi",Tau_phi);
    tin->SetBranchAddress("Tau_dz",Tau_dz);
    tin->SetBranchAddress("Tau_charge",Tau_charge);
    tin->SetBranchAddress("Tau_dxy",Tau_dxy);

    TLorentzVector p4mu1, p4mu2,p4dimuon; //lorentz vectors for muons

    ifstream rootdirec((std::string("rootfiles/") + png_dir + std::string("/hfile.root")).c_str()); //check if directories already exists
    if (!rootdirec) gSystem->Exec((std::string("mkdir -p rootfiles/" + png_dir)).c_str()); //creates the directory

    TFile *fout = new TFile((std::string("rootfiles/") + png_dir + std::string("/hfile.root")).c_str(),"RECREATE");

    //TH1F -> Maximum Precision of 7 digits per histogram
    TH1F **hdimu = new TH1F*[4]; //double muon distributions
    TH1F **hmu = new TH1F*[5]; //single muon distributions
    TH1F **hz1 = new TH1F*[4]; //distributions of Z boson -> custom cuts
    TH1F **hz2 = new TH1F*[4]; //distributions of Z boson -> atlas open data cuts
    TH1F **hz3 = new TH1F*[4]; //distributions of Z boson -> cms open data cuts (fabian)
    TH1F **hz4 = new TH1F*[4]; //distributions of Z boson -> higgs analysis (paula)
    TH1F **hjpsi1 = new TH1F*[4]; //distributions of J/Psi -> custom cuts
    TH1F **hjpsi2 = new TH1F*[4]; //distributions of J/Psi -> global muons cms open data
    TH1F **hjpsi3 = new TH1F*[4]; //distributions of J/Psi -> tracker muons cms open data
    TH1F **hjpsi4 = new TH1F*[4]; //distributions of J/Psi -> zeus cuts (raphael)

    //distributions
    char hdimutitle[100]; //title of dimuon histograms 
    char hdimuname[100]; //name of dimuon histograms
    char hdimutitles[4][100] = {" pT Distribution - Double Muon -"," Rapidity Distribution - Double Muon -"," #phi Distribution - Double Muon -","  M_{#mu #mu} Distribution - Double Muon -"};
    char hdimuxaxis[4][100] = {" pT [GeV/c] ","Rapidity","#phi"," M_{#mu #mu} [GeV/c^{2}]"};
    char hdimuyaxis[4][100] = {"Events / 1 GeV","Events / 0.1 ", "Events / 0.1 ","Events / 0.1 GeV"};
    double hdimurngmin[4] = {-0.5, -3.05, -3.55, -0.05};
    double hdimurngmax[4] = {120.5, 3.05, 3.55, 200.05};
    int hdimubins[4] = {121, 61, 71, 2001};

    char hmutitle[100]; //title of single muon histograms 
    char hmuname[100]; //name of single muon histograms
    char hmutitles[5][100] = {" pT Distribution - Single Muon -"," Pseudorapidity Distribution - Single Muon -"," #phi Distribution - Single Muon -", " pfRelIso03_chg Distribution - Z -", " pfRelIso03_all Distribution - Z -"};
    char hmuxaxis[5][100] = {" pT [GeV/c] ","#eta","#phi","pfRelIso03_chg","pfRelIso03_all"};
    char hmuyaxis[5][100] = {"Events / 1 GeV","Events / 0.1 ", "Events / 0.1 ","Events / 0.01","Events / 0.01"};
    double hmurngmin[5] = {-0.5, -3.05, -3.55, -0.105, -0.105};
    double hmurngmax[5] = {120.5, 3.05, 3.55, 3.005, 3.005};
    int hmubins[5] = {121, 61, 71, 311, 311};

    char hztitle[100]; //title of z histograms 
    char hztitles[4][100] = {" pT Distribution - Z -"," Rapidity Distribution - Z -"," #phi Distribution - Z -","  M_{#mu #mu} - Z -"};
    char hzxaxis[4][100] = {"pT [GeV/c]","Rapidity","#phi","M_{#mu #mu} [GeV/c^{2}]"};
    char hzyaxis[4][100] = {"Events / 1 GeV","Events / 0.1", "Events / 0.1","Events / 0.1 GeV"};
    double hzrngmin[4] = {-0.5, -3.05, -3.55, 59.95};
    double hzrngmax[4] = {120.5, 3.05, 3.55, 120.05};
    int hzbins[4] = {121, 61, 71, 601};

    char hz1name[100]; //name of z histograms -> cut 1
    char hz2name[100]; //name of z histograms -> cut 2
    char hz3name[100]; //name of z histograms -> cut 3
    char hz4name[100]; //name of z histograms -> cut 4

    char hjpsititle[100]; //title of j/psi histograms 
    char hjpsititles[4][100] = {" pT Distribution - J/#Psi -"," Rapidity Distribution - J/#Psi -"," #phi Distribution - J/#Psi -","  M_{#mu #mu} - J/#Psi -"};
    char hjpsixaxis[4][100] = {" pT [GeV/c] ","Rapidity","#phi"," M_{#mu #mu} [GeV/c^{2}]"};
    char hjpsiyaxis[4][100] = {"Events / 1 GeV","Events / 0.1", "Events / 0.1","Events / 0.01 GeV"};
    double hjpsirngmin[4] = {-0.5, -3.05, -3.55, 1.995};
    double hjpsirngmax[4] = {120.5, 3.05, 3.55, 5.005};
    int hjpsibins[4] = {121, 61, 71, 301};

    char hjpsi1name[100]; //name of j/psi histograms -> cut 1
    char hjpsi2name[100]; //name of j/psi histograms -> cut 2
    char hjpsi3name[100]; //name of j/psi histograms -> cut 3
    char hjpsi4name[100]; //name of j/psi histograms -> cut 4

TH1F* hjpsipt = new TH1F("hjpsipt", "p_{T} - J/Psi - ZEUS DATA", 50, 0, 10); 

    for(int i=0;i<5;i++) { //double muon distributions

        if (i<4) { 
            
            //double muon distributions
            sprintf(hdimuname,"hdimu%d",i);
            strcpy(hdimutitle,hdimutitles[i]); //assigning 
            strncat(hdimutitle,legend_title.c_str(),50); //adding legend to title
            hdimu[i] = new TH1F(hdimuname,hdimutitle,hdimubins[i],hdimurngmin[i],hdimurngmax[i]); //TH1F declaration
            hdimu[i]->GetXaxis()->SetTitle(hdimuxaxis[i]);
            hdimu[i]->GetYaxis()->SetTitle(hdimuyaxis[i]);
            if (mcdata==0) hdimu[i]->SetLineColor(kBlue);
            if (mcdata==1) hdimu[i]->SetLineColor(kRed);

            //Z distributions
            sprintf(hz1name,"hz%d_cut1",i);
            sprintf(hz2name,"hz%d_cut2",i);
            sprintf(hz3name,"hz%d_cut3",i);
            sprintf(hz4name,"hz%d_cut4",i);
            strcpy(hztitle,hztitles[i]); //assigning 
            strncat(hztitle,legend_title.c_str(),50); //adding legend to title
            hz1[i] = new TH1F(hz1name,hztitle,hzbins[i],hzrngmin[i],hzrngmax[i]); //TH1F declaration
            hz1[i]->GetXaxis()->SetTitle(hzxaxis[i]);
            hz1[i]->GetYaxis()->SetTitle(hzyaxis[i]);
            if (mcdata==0) hz1[i]->SetLineColor(kBlue);
            if (mcdata==1) hz1[i]->SetLineColor(kRed);
            hz2[i] = new TH1F(hz2name,hztitle,hzbins[i],hzrngmin[i],hzrngmax[i]); //TH1F declaration
            hz2[i]->GetXaxis()->SetTitle(hzxaxis[i]);
            hz2[i]->GetYaxis()->SetTitle(hzyaxis[i]);
            if (mcdata==0) hz2[i]->SetLineColor(kBlue);
            if (mcdata==1) hz2[i]->SetLineColor(kRed);
            hz3[i] = new TH1F(hz3name,hztitle,hzbins[i],hzrngmin[i],hzrngmax[i]); //TH1F declaration
            hz3[i]->GetXaxis()->SetTitle(hzxaxis[i]);
            hz3[i]->GetYaxis()->SetTitle(hzyaxis[i]);
            if (mcdata==0) hz3[i]->SetLineColor(kBlue);
            if (mcdata==1) hz3[i]->SetLineColor(kRed);
            hz4[i] = new TH1F(hz4name,hztitle,hzbins[i],hzrngmin[i],hzrngmax[i]); //TH1F declaration
            hz4[i]->GetXaxis()->SetTitle(hzxaxis[i]);
            hz4[i]->GetYaxis()->SetTitle(hzyaxis[i]);
            if (mcdata==0) hz4[i]->SetLineColor(kBlue);
            if (mcdata==1) hz4[i]->SetLineColor(kRed);

            //jpsi distributions
            sprintf(hjpsi1name,"hjpsi%d_cut1",i);
            sprintf(hjpsi2name,"hjpsi%d_cut2",i);
            sprintf(hjpsi3name,"hjpsi%d_cut3",i);
            sprintf(hjpsi4name,"hjpsi%d_cut4",i);
            strcpy(hjpsititle,hjpsititles[i]); //assigning 
            strncat(hjpsititle,legend_title.c_str(),50); //adding legend to title
            hjpsi1[i] = new TH1F(hjpsi1name,hjpsititle,hjpsibins[i],hjpsirngmin[i],hjpsirngmax[i]); //TH1F declaration
            hjpsi1[i]->GetXaxis()->SetTitle(hjpsixaxis[i]);
            hjpsi1[i]->GetYaxis()->SetTitle(hjpsiyaxis[i]);
            if (mcdata==0) hjpsi1[i]->SetLineColor(kBlue);
            if (mcdata==1) hjpsi1[i]->SetLineColor(kRed);
            hjpsi2[i] = new TH1F(hjpsi2name,hjpsititle,hjpsibins[i],hjpsirngmin[i],hjpsirngmax[i]); //TH1F declaration
            hjpsi2[i]->GetXaxis()->SetTitle(hjpsixaxis[i]);
            hjpsi2[i]->GetYaxis()->SetTitle(hjpsiyaxis[i]);
            if (mcdata==0) hjpsi2[i]->SetLineColor(kBlue);
            if (mcdata==1) hjpsi2[i]->SetLineColor(kRed);
            hjpsi3[i] = new TH1F(hjpsi3name,hjpsititle,hjpsibins[i],hjpsirngmin[i],hjpsirngmax[i]); //TH1F declaration
            hjpsi3[i]->GetXaxis()->SetTitle(hjpsixaxis[i]);
            hjpsi3[i]->GetYaxis()->SetTitle(hjpsiyaxis[i]);
            if (mcdata==0) hjpsi3[i]->SetLineColor(kBlue);
            if (mcdata==1) hjpsi3[i]->SetLineColor(kRed);
            hjpsi4[i] = new TH1F(hjpsi4name,hjpsititle,hjpsibins[i],hjpsirngmin[i],hjpsirngmax[i]); //TH1F declaration
            hjpsi4[i]->GetXaxis()->SetTitle(hjpsixaxis[i]);
            hjpsi4[i]->GetYaxis()->SetTitle(hjpsiyaxis[i]);
            if (mcdata==0) hjpsi4[i]->SetLineColor(kBlue);
            if (mcdata==1) hjpsi4[i]->SetLineColor(kRed);

        }

        //single muon distributions
        sprintf(hmuname,"hmu%d",i);
        strcpy(hmutitle,hmutitles[i]); //assigning
        strncat(hmutitle,legend_title.c_str(),50); //adding legend to title
        hmu[i] = new TH1F(hmuname,hmutitle,hmubins[i],hmurngmin[i],hmurngmax[i]); //TH1F declaration
        hmu[i]->GetXaxis()->SetTitle(hmuxaxis[i]);
        hmu[i]->GetYaxis()->SetTitle(hmuyaxis[i]);
        if (mcdata==0) hmu[i]->SetLineColor(kBlue);
        if (mcdata==1) hmu[i]->SetLineColor(kRed);

    }

    //nentries = 1e6; //test

    //reading tin
    cout<<"- Entries # "<<nentries<<endl;

    int numlepmax = 0; //reads data and finds maximum of lep_n so that we can put an upper limit
    int numjetmax = 0;
    int singleentries = 0;
    int doubleentries = 0;
    int zentries = 0;
    int jpsientries = 0;

    for (int i=0;i<nentries;i++) {

        tin->GetEntry(i); //reading
        if (i%((int)(nentries/10))==0) cout<<"- Entry # "<<i<<endl;
        if ((int)nMuon>numlepmax) numlepmax = nMuon;
        if ((int)nElectron>numlepmax) numlepmax = nElectron;
        if ((int)nTau>numlepmax) numlepmax = nTau;

        for (int m1=0; m1<(int)nMuon && m1<maxlepton_n; m1++) { //loop on first muon

            if (m1 == maxlepton_n - 1) cout<<"- WARNING: event # "<<i<<" exceeded maximum size of arrays "<<maxlepton_n<<" < "<<nMuon<<endl;
            singleentries++;

            /*
            if (event==54469157) {
                cout<<" | event | run | luminosityBlock | nMuon | Muon_pt | Muon_eta | Muon_phi | Muon_pfRelIso03_chg | Muon_pfRelIso03_all | PV_npvs | TrigObj_id | Muon_charge | PV_z |"<<endl;
                cout<<" | "<<event<<" | "<<run<<" | "<<luminosityBlock<<" | "<<nMuon<<" | "<<Muon_pt[m1]<<" | "<<Muon_eta[m1]<<" | "<<Muon_phi[m1]<<" | "<<Muon_pfRelIso03_chg[m1]<<" | "<<Muon_pfRelIso03_all[m1]<<" | "<<PV_npvs<<" | "<<TrigObj_id<<" | "<<Muon_charge[m1]<<" | "<<PV_z<<" | "<<endl;
                //if (m1==nMuon-1) return;                
            }
            */

            //filling lorentz vector of first muon
            p4mu1.SetPtEtaPhiM(Muon_pt[m1], Muon_eta[m1],Muon_phi[m1],muonmass); //GeV units

            hmu[0]->Fill(Muon_pt[m1]);
            hmu[1]->Fill(Muon_eta[m1]);
            hmu[2]->Fill(Muon_phi[m1]);
            //hmu[3]->Fill(Muon_pfRelIso03_chg[m1]);
            //hmu[4]->Fill(Muon_pfRelIso03_all[m1]);

            for (int m2=m1+1; m2<(int)nMuon && m2<maxlepton_n; m2++) { //loop on second muon

                doubleentries++;

		float dphimus = Muon_phi[m1] - Muon_phi[m2];
   		        if (dphimus > PI) dphimus = dphimus - 2.*PI;  
      	                if (dphimus < -PI) dphimus = dphimus + 2.*PI;

			if ( (abs( abs(dphimus) - PI) < 5.*PI/100.) && (abs(Muon_eta[m1] + Muon_eta[m2]) < 5.*PI/100.)) continue; 

                //filling lorentz vector of second muon
                p4mu2.SetPtEtaPhiM(Muon_pt[m2], Muon_eta[m2], Muon_phi[m2],muonmass);
                p4dimuon = p4mu1 + p4mu2; //dimuon lorentz vector

                //filling double muon histograms
                hdimu[0]->Fill(p4dimuon.Pt());
                hdimu[1]->Fill(p4dimuon.Rapidity());
                hdimu[2]->Fill(p4dimuon.Phi());
                hdimu[3]->Fill(p4dimuon.M());

                //filling Z histograms -> custom cuts 
                if ( ((Muon_pt[m1]>5 && Muon_pt[m2]>25) || (Muon_pt[m1]>25 && Muon_pt[m2]>5)) 
                    && (Bool_t)Muon_isGlobal[m1] && (Bool_t)Muon_isGlobal[m2] 
                    && (int)(Muon_charge[m1] + Muon_charge[m2])==0) { //Z cuts #1 -> custom

                    if (p4dimuon.M()>80 && p4dimuon.M()<100) {
                        zentries++;
                        hz1[0]->Fill(p4dimuon.Pt());
                        hz1[1]->Fill(p4dimuon.Rapidity());
                        hz1[2]->Fill(p4dimuon.Phi());
                        hmu[3]->Fill(Muon_pfRelIso03_chg[m1]); //Relative Isolation for Z
                        hmu[4]->Fill(Muon_pfRelIso03_all[m1]);
                        hmu[3]->Fill(Muon_pfRelIso03_chg[m2]);
                        hmu[4]->Fill(Muon_pfRelIso03_all[m2]);
                    }

                    hz1[3]->Fill(p4dimuon.M());
                }

                //filling Z histograms -> ATLAS 2012 cuts 
                if (Muon_pt[m1]>25 && Muon_pt[m2]>25 
                    && Muon_pfRelIso03_all[m1]<0.15 && Muon_pfRelIso03_all[m2]<0.15 
                    && Muon_pfRelIso03_chg[m1]<0.15 && Muon_pfRelIso03_chg[m2]<0.15 
                    && (int)(Muon_charge[m1] + Muon_charge[m2])==0) { //Z cuts #2 -> ATLAS
                    
                    if (p4dimuon.M()>80 && p4dimuon.M()<100) {
                        hz2[0]->Fill(p4dimuon.Pt());
                        hz2[1]->Fill(p4dimuon.Rapidity());
                        hz2[2]->Fill(p4dimuon.Phi());    
                    }

                    hz2[3]->Fill(p4dimuon.M());
                }

                //filling Z histograms -> Fabian Cuts (CMS open data)
                if (Muon_pt[m1]>20 && Muon_pt[m2]>20 
                    && Muon_pfRelIso03_all[m1]<0.15 && Muon_pfRelIso03_all[m2]<0.15 
                    && TMath::Abs(Muon_eta[m1])<=2.1 && TMath::Abs(Muon_eta[m2])<=2.1
                    && TMath::Abs(Muon_dxy[m1])<=0.2 && TMath::Abs(Muon_dxy[m2])<=0.2
                    && (Bool_t)Muon_isGlobal[m1] && (Bool_t)Muon_isGlobal[m2]
                    && (int)(Muon_charge[m1] + Muon_charge[m2])==0) { //Z cuts #3 -> CMS Fabian (Open data)
                    
                    if (p4dimuon.M()>80 && p4dimuon.M()<100) {
                        hz3[0]->Fill(p4dimuon.Pt());
                        hz3[1]->Fill(p4dimuon.Rapidity());
                        hz3[2]->Fill(p4dimuon.Phi());
                    }

                    hz3[3]->Fill(p4dimuon.M());
                }

                if (Muon_pt[m1]>5 && Muon_pt[m2]>5 
                    && Muon_pfRelIso04_all[m1]<0.4 && Muon_pfRelIso04_all[m2]<0.4 
                    && TMath::Abs(Muon_eta[m1])<=2.4 && TMath::Abs(Muon_eta[m2])<=2.4
                    && TMath::Abs(Muon_dxy[m1])<=0.5 && TMath::Abs(Muon_dxy[m2])<=0.5
                    && TMath::Abs(Muon_dz[m1])<=1 && TMath::Abs(Muon_dz[m2])<=1
                    && TMath::Abs(Muon_sip3d[m1])<=4 && TMath::Abs(Muon_sip3d[m2])<=4
                    && (Bool_t)Muon_isGlobal[m1] && (Bool_t)Muon_isGlobal[m2]
                    && (Bool_t)Muon_isPFcand[m1] && (Bool_t)Muon_isPFcand[m2]
                    && (int)(Muon_charge[m1] + Muon_charge[m2])==0) { //Z cuts #4 -> Higgs Analysis (Paula)
                    
                    if (p4dimuon.M()>80 && p4dimuon.M()<100) {
                        hz4[0]->Fill(p4dimuon.Pt());
                        hz4[1]->Fill(p4dimuon.Rapidity());
                        hz4[2]->Fill(p4dimuon.Phi());
                    }

                    hz4[3]->Fill(p4dimuon.M());
                }

                //filling J/Psi histograms -> Custom cuts
                if ( ((Muon_pt[m1]>5 && Muon_pt[m2]>25) || (Muon_pt[m1]>25 && Muon_pt[m2]>5)) 
                    && (Bool_t)Muon_softId[m1] && (Bool_t)Muon_softId[m2] 
                    && (Bool_t)Muon_isGlobal[m1] && (Bool_t)Muon_isGlobal[m2] 
                    && (int)(Muon_charge[m1] + Muon_charge[m2])==0) { //J/Psi cuts -> neutral
                    
                    if (p4dimuon.M()>2.8 && p4dimuon.M()<3.4) {
                        jpsientries++;
                        hjpsi1[0]->Fill(p4dimuon.Pt());
                        hjpsi1[1]->Fill(p4dimuon.Rapidity());
                        hjpsi1[2]->Fill(p4dimuon.Phi());
                    }

                    hjpsi1[3]->Fill(p4dimuon.M());
                }

                //filling J/Psi histograms -> Global Muons (from Fabian)
                if ((Bool_t)Muon_isGlobal[m1] && (Bool_t)Muon_isGlobal[m2] 
                    && TMath::Abs(p4dimuon.Rapidity()) < 1.2 //1.6/2.4
                    && (int)(Muon_charge[m1] + Muon_charge[m2])==0) { //J/Psi cuts -> neutral

                    if (p4dimuon.M()>2.8 && p4dimuon.M()<3.4) {
                        hjpsi2[0]->Fill(p4dimuon.Pt());
                        hjpsi2[1]->Fill(p4dimuon.Rapidity());
                        hjpsi2[2]->Fill(p4dimuon.Phi());
                    }
                    hjpsi2[3]->Fill(p4dimuon.M());
                }

                //filling J/Psi histograms -> Tracker Muons (from Fabian)
                if ((Bool_t)Muon_isTracker[m1] && (Bool_t)Muon_isTracker[m2] 
                    && ((TMath::Abs(Muon_eta[m1])<=1.3 && Muon_pt[m1]>3.3) || (TMath::Abs(Muon_eta[m1])>1.3 && TMath::Abs(Muon_eta[m1])<=2.2 && p4mu1.P()>2.9) || (TMath::Abs(Muon_eta[m1])>2.2 && TMath::Abs(Muon_eta[m1])<=2.3 && Muon_pt[m1]>0.8)) 
                    && ((TMath::Abs(Muon_eta[m2])<=1.3 && Muon_pt[m2]>3.3) || (TMath::Abs(Muon_eta[m2])>1.3 && TMath::Abs(Muon_eta[m2])<=2.2 && p4mu2.P()>2.9) || (TMath::Abs(Muon_eta[m2])>2.2 && TMath::Abs(Muon_eta[m2])<=2.3 && Muon_pt[m2]>0.8)) 
                    && TMath::Abs(p4dimuon.Rapidity()) < 1.2 //1.6/2.4
                    && (int)(Muon_charge[m1] + Muon_charge[m2])==0) { //J/Psi cuts -> neutral
                    
                    if (p4dimuon.M()>2.8 && p4dimuon.M()<3.4) {
                        hjpsi3[0]->Fill(p4dimuon.Pt());
                        hjpsi3[1]->Fill(p4dimuon.Rapidity());
                        hjpsi3[2]->Fill(p4dimuon.Phi());
                    }
                    hjpsi3[3]->Fill(p4dimuon.M());
                }

                if (Muon_pt[m1]>1 && Muon_pt[m2]>1
                    && (Bool_t)Muon_softId[m1] && (Bool_t)Muon_softId[m2] 
                    && (Bool_t)Muon_isGlobal[m1] && (Bool_t)Muon_isGlobal[m2] 
                    && (int)(Muon_charge[m1] + Muon_charge[m2])==0) { //J/Psi cuts -> Zeus (raphael)


                    
                    if (p4dimuon.M()>2.8 && p4dimuon.M()<3.4) {
			//hjpsipt->Fill(p4dimuon.Pt());
                        hjpsi4[0]->Fill(p4dimuon.Pt());
                        hjpsi4[1]->Fill(p4dimuon.Rapidity());
                        hjpsi4[2]->Fill(p4dimuon.Phi());
/*
if (abs(p4dimuon.Rapidity())<0.05) {
                cout<<" | event | run | luminosityBlock | nMuon | Muon_pt | Muon_eta | Muon_phi | Muon_pfRelIso03_chg | Muon_pfRelIso03_all | PV_npvs | PV_chi2 | Muon_charge | PV_z |"<<endl;
                cout<<" | "<<event<<" | "<<run<<" | "<<luminosityBlock<<" | "<<nMuon<<" | "<<Muon_pt[m1]<<" | "<<Muon_eta[m1]<<" | "<<Muon_phi[m1]<<" | "<<Muon_pfRelIso03_chg[m1]<<" | "<<Muon_pfRelIso03_all[m1]<<" | "<<PV_npvs<<" | "<<PV_chi2<<" | "<<Muon_charge[m1]<<" | "<<PV_z<<" | "<<endl;
cout<<" | event | run | luminosityBlock | nMuon | Muon_pt | Muon_eta | Muon_phi | Muon_pfRelIso03_chg | Muon_pfRelIso03_all | PV_npvs | PV_chi2 | Muon_charge | PV_z |"<<endl;
                cout<<" | "<<event<<" | "<<run<<" | "<<luminosityBlock<<" | "<<nMuon<<" | "<<Muon_pt[m2]<<" | "<<Muon_eta[m2]<<" | "<<Muon_phi[m2]<<" | "<<Muon_pfRelIso03_chg[m2]<<" | "<<Muon_pfRelIso03_all[m2]<<" | "<<PV_npvs<<" | "<<PV_chi2<<" | "<<Muon_charge[m2]<<" | "<<PV_z<<" | "<<endl;
cout << p4dimuon.Pt() << " " << p4dimuon.Rapidity() << " " << p4dimuon.Phi() << " " << p4dimuon.M() << endl;
                //if (m1==nMuon-1) return;                
            }
*/

                    }

                    hjpsi4[3]->Fill(p4dimuon.M());
                }

            }//end loop on second muon
        }//end loop on first muon
    }//end loop on entries

    cout<<"- Maximum Array Size = "<<numlepmax<<endl;
    cout<<"- Single Muon Entries = "<<singleentries<<endl;
    cout<<"- Double Muon Entries = "<<doubleentries<<endl;
    cout<<"- Z Entries = "<<zentries<<endl;
    cout<<"- J/Psi Entries = "<<jpsientries<<endl;
    cout<<"- Peak Eta (BinLowEdge,BinUpEdge,Entries) = "<<hmu[1]->GetBinLowEdge(hmu[1]->GetMaximumBin())<<" , "<<hmu[1]->GetBinLowEdge(hmu[1]->GetMaximumBin()) + hmu[1]->GetBinWidth(hmu[1]->GetMaximumBin())<<" , "<<hmu[1]->GetBinContent(hmu[1]->GetMaximumBin())<<endl;
    cout<<"- Peak Phi (BinLowEdge,BinUpEdge,Entries) = "<<hmu[2]->GetBinLowEdge(hmu[2]->GetMaximumBin())<<" , "<<hmu[2]->GetBinLowEdge(hmu[2]->GetMaximumBin()) + hmu[2]->GetBinWidth(hmu[2]->GetMaximumBin())<<" , "<<hmu[2]->GetBinContent(hmu[2]->GetMaximumBin())<<endl;

    //plotting and saving canvas -> maybe put in a different code
    gROOT->SetBatch(kTRUE); //no pop up
    gStyle->SetOptStat("nemruo"); //name, entries, mean, rms, underflow, overflow
    //gStyle->SetOptStat(0);

    TCanvas** cmu = new TCanvas*[5]; //array of canvas for single muons
    TCanvas** cdimu = new TCanvas*[4]; //array of canvas for double muons
    TCanvas** cz1 = new TCanvas*[4]; //array of canvas for Z
    TCanvas** cz2 = new TCanvas*[4]; //array of canvas for Z
    TCanvas** cz3 = new TCanvas*[4]; //array of canvas for Z
    TCanvas** cz4 = new TCanvas*[4]; //array of canvas for Z
    TCanvas** cjpsi1 = new TCanvas*[4]; //array of canvas for J/Psi
    TCanvas** cjpsi2 = new TCanvas*[4]; //array of canvas for J/Psi
    TCanvas** cjpsi3 = new TCanvas*[4]; //array of canvas for J/Psi
    TCanvas** cjpsi4 = new TCanvas*[4]; //array of canvas for J/Psi
    char ctitle[50]; //title of canvas

    TLine *lz_l = new TLine(80,0,80,hz1[3]->GetMaximum()/2.);
    TLine *lz_r = new TLine(100,0,100,hz1[3]->GetMaximum()/2.);
    TLine *ljpsi_l = new TLine(2.8,0,2.8,hjpsi1[3]->GetMaximum()/2.);
    TLine *ljpsi_r = new TLine(3.4,0,3.4,hjpsi1[3]->GetMaximum()/2.);

    ifstream direc((std::string("plots/") + png_dir + std::string("/hmu0.png")).c_str()); //check if directories already exists
    if (!direc) gSystem->Exec((std::string("mkdir -p plots/" + png_dir)).c_str()); //creates the directory (plots has to already exist)

    for (int i=0;i<5;i++){ //delete histograms avoid memory leaks

        if (i<4) {

            //double muon distributions
            cdimu[i]=new TCanvas();
            sprintf(ctitle,"/hdimu%d.png",i); //name of file
            hdimu[i]->Draw("hist");
            if (i==3) { //full dimuon mass spectrum
                cdimu[i]->SetLogy();
                cdimu[i]->SetLogx();
            }
            cdimu[i]->SaveAs((std::string("plots/") + png_dir+ctitle).c_str());
            cdimu[i]->Close();
            //delete hdimu[i];

            //Z distributions
            cz1[i]=new TCanvas();
            sprintf(ctitle,"/hz%d_cuts1.png",i); //name of file
            hz1[i]->Draw("hist");
            if (i==3) {
                lz_l->Draw("same");
                lz_r->Draw("same");
            }
            cz1[i]->SaveAs((std::string("plots/") + png_dir + ctitle).c_str());
            cz1[i]->Close();
            //delete hz1[i];

            cz2[i]=new TCanvas();
            sprintf(ctitle,"/hz%d_cuts2.png",i); //name of file
            hz2[i]->Draw("hist");
            if (i==3) {
                lz_l->Draw("same");
                lz_r->Draw("same");
            }
            cz2[i]->SaveAs((std::string("plots/") + png_dir + ctitle).c_str());
            cz2[i]->Close();
            //delete hz2[i];

            cz3[i]=new TCanvas();
            sprintf(ctitle,"/hz%d_cuts3.png",i); //name of file
            hz3[i]->Draw("hist");
            if (i==3) {
                lz_l->Draw("same");
                lz_r->Draw("same");
            }
            cz3[i]->SaveAs((std::string("plots/") + png_dir + ctitle).c_str());
            cz3[i]->Close();
            //delete hz3[i];

            cz4[i]=new TCanvas();
            sprintf(ctitle,"/hz%d_cuts4.png",i); //name of file
            hz4[i]->Draw("hist");
            if (i==3) {
                lz_l->Draw("same");
                lz_r->Draw("same");
            }
            cz4[i]->SaveAs((std::string("plots/") + png_dir + ctitle).c_str());
            cz4[i]->Close();
            //delete hz4[i];

            //j/psi distributions
            cjpsi1[i]=new TCanvas();
            sprintf(ctitle,"/hjpsi%d_cuts1.png",i); //name of file
            hjpsi1[i]->Draw("hist");
            if (i==3) {
                ljpsi_l->Draw("same");
                ljpsi_r->Draw("same");
            }
            cjpsi1[i]->SaveAs((std::string("plots/") + png_dir + ctitle).c_str());
            cjpsi1[i]->Close();
            //delete hjpsi1[i];

            cjpsi2[i]=new TCanvas();
            sprintf(ctitle,"/hjpsi%d_cuts2.png",i); //name of file
            hjpsi2[i]->Draw("hist");
            if (i==3) {
                ljpsi_l->Draw("same");
                ljpsi_r->Draw("same");
            }
            cjpsi2[i]->SaveAs((std::string("plots/") + png_dir + ctitle).c_str());
            cjpsi2[i]->Close();
            //delete hjpsi2[i];

            cjpsi3[i]=new TCanvas();
            sprintf(ctitle,"/hjpsi%d_cuts3.png",i); //name of file
            hjpsi3[i]->Draw("hist");
            if (i==3) {
                ljpsi_l->Draw("same");
                ljpsi_r->Draw("same");
            }
            cjpsi3[i]->SaveAs((std::string("plots/") + png_dir + ctitle).c_str());
            cjpsi3[i]->Close();
            //delete hjpsi3[i];

            cjpsi4[i]=new TCanvas();
            sprintf(ctitle,"/hjpsi%d_cuts4.png",i); //name of file
            hjpsi4[i]->Draw("hist");
	    //hjpsipt->Draw("hist");
            if (i==3) {
                ljpsi_l->Draw("same");
                ljpsi_r->Draw("same");
            }
            cjpsi4[i]->SaveAs((std::string("plots/") + png_dir + ctitle).c_str());
            cjpsi4[i]->Close();
            //delete hjpsi4[i];

        }

        //single muon distributions
        cmu[i]=new TCanvas();
        sprintf(ctitle,"/hmu%d.png",i); //name of file
        hmu[i]->Draw("hist");
        if (i==3 || i==4) cmu[i]->SetLogy(); //pfreliso
        cmu[i]->SaveAs((std::string("plots/") + png_dir + ctitle).c_str());
        cmu[i]->Close();
        //delete hmu[i];
    }

    fout->Write();
    fout->Close();

    //performance
    cout<<"- Timer: ";
    t.Stop();
    t.Print();

    return;

}
