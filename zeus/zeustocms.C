// system include files
#include <memory>
#include <iostream>
#include "stdlib.h"
// for histogramming
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TStyle.h"
#include "TLorentzVector.h"
#include "TTree.h" //
#include "TMath.h"
#include "TFile.h"
//#include "Math/VectorUtil_cint.h"

#include <vector>	
#include <algorithm>
#include <utility>
#include <stdexcept>

#include "TMatrixT.h"
#include "TMath.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TGaxis.h"
#include "TChain.h"
//#include "Math/VectorUtil.h"
//#include "Math/GenVector/VectorUtil.h"

#define PiConstant TMath::Pi()
#define DoublePiConstant TMath::TwoPi()

using namespace std;

void zeustocms() {

	printf("=========================== Start =====================================\n");
	//Maximum of Nmu. We need this to declare our arrays before we know each Nmu.
	const Int_t Nmax = 20;
	const Int_t trackNMax = 300;

	//Declare constants
	
	//ZEUS variables
	Int_t Runnr, Eventnr, Ntrkvtx, Trk_ntracks, Trk_id[Nmax]; 
	Float_t Xvtx, Yvtx, Zvtx, Chivtx;

	Int_t Nmu;
	Int_t Mucharge[Nmax], Muqual[Nmax], Mutrid[Nmax], Muzufid[Nmax], Muchid[Nmax], Muvtxid[Nmax], Muvtxfl[Nmax], Mujetid_a[Nmax]; 
	Float_t Mupt[Nmax], Muth[Nmax], Muph[Nmax], Muperr[Nmax], Mup[Nmax][3], Muz[Nmax], Mudxy[Nmax], Mudz[Nmax], Muisol[Nmax][10];

	Float_t Bspt_x, Bspt_y, Bspt_z, Bspt_xer, Bspt_yer, Bspt_zer, Bspt_dxdz, Bspt_dydz;


	//Electron variables
	Int_t Emncand, Emtrknr[Nmax]; // electric charge
	Float_t Empos[Nmax][3], Emcalpos[Nmax][3], Emtrkq[Nmax], Emetneartrk[Nmax][2], Emdcabeam[Nmax], Emdca[Nmax], Emtrkth[Nmax], Emenin[Nmax], Emtrkph[Nmax];
	Float_t Empt[Nmax], Trk_pca[trackNMax][3], Emcalene[Nmax], Emtrkp[Nmax], Emprob[Nmax], Emph[Nmax];
	Float_t Emxel[Nmax], Emyel[Nmax], Emq2el[Nmax], Emxda[Nmax], Emyda[Nmax], Emq2da[Nmax], Emxda_cell[Nmax], Emyda_cell[Nmax], Emq2da_cell[Nmax], Emfemc[Nmax]; 
	

	//CMS variables
	Int_t run, event, luminosityBlock, nPVtx, nMuon;

	nPVtx = 100;
	
	//Get the demanded file
	TChain in_chain("orange");
	in_chain.Add("/pnfs/desy.de/dphep/online/zeus/z/ntup/07p/v08b/data/root/data_07p*.root");	//good ntuple file to check code fast: data_07p_61792_61792_01.root
	// in_chain.Add("/pnfs/desy.de/dphep/online/zeus/z/ntup/06p/v08b/data/root/data_*.root");	//good ntuple file to check code fast: data_07p_61792_61792_01.root
	// in_chain.Add("/pnfs/desy.de/dphep/online/zeus/z/ntup/06e/v08b/data/root/data_*.root");	//good ntuple file to check code fast: data_07p_61792_61792_01.root
	// in_chain.Add("/pnfs/desy.de/dphep/online/zeus/z/ntup/05/v08b/data/root/data_*.root");	//good ntuple file to check code fast: data_07p_61792_61792_01.root

  	in_chain.SetBranchAddress("Runnr", &Runnr);
  	in_chain.SetBranchAddress("Eventnr", &Eventnr);
	in_chain.SetBranchAddress("Trk_ntracks", &Trk_ntracks);
	in_chain.SetBranchAddress("Trk_id", &Trk_id);

	//
	in_chain.SetBranchAddress("Xvtx", &Xvtx);
  	in_chain.SetBranchAddress("Yvtx", &Yvtx);
	in_chain.SetBranchAddress("Zvtx", &Zvtx);
  	in_chain.SetBranchAddress("Chivtx", &Chivtx);
	in_chain.SetBranchAddress("Ntrkvtx", &Ntrkvtx);
	
	//Muon data
	in_chain.SetBranchAddress("Nmu", &Nmu);
  	in_chain.SetBranchAddress("Mucharge", &Mucharge);
	in_chain.SetBranchAddress("Mupt", &Mupt);
  	in_chain.SetBranchAddress("Muth", &Muth);
	in_chain.SetBranchAddress("Muph", &Muph);
	in_chain.SetBranchAddress("Muqual", &Muqual);
  	in_chain.SetBranchAddress("Muperr", &Muperr);
	in_chain.SetBranchAddress("Mup", &Mup);
  	in_chain.SetBranchAddress("Muzufid", &Muzufid);
	in_chain.SetBranchAddress("Muchid", &Muchid);
	in_chain.SetBranchAddress("Muz", &Muz);
  	in_chain.SetBranchAddress("Mudxy", &Mudxy);
	in_chain.SetBranchAddress("Mudz", &Mudz);
  	in_chain.SetBranchAddress("Muisol", &Muisol);
	in_chain.SetBranchAddress("Muvtxid", &Muvtxid);
	in_chain.SetBranchAddress("Muvtxfl", &Muvtxfl);
  	in_chain.SetBranchAddress("Mujetid_a", &Mujetid_a);

	// Electron data
	in_chain.SetBranchAddress("Emprob", &Emprob);
	in_chain.SetBranchAddress("Emncand", &Emncand);
	in_chain.SetBranchAddress("Emcalpos", &Emcalpos);
	in_chain.SetBranchAddress("Empos", &Empos);
	in_chain.SetBranchAddress("Emetneartrk", &Emetneartrk);
	in_chain.SetBranchAddress("Emdcabeam", &Emdcabeam);
	in_chain.SetBranchAddress("Emdca", &Emdca);
	in_chain.SetBranchAddress("Emtrkth", &Emtrkth);
	in_chain.SetBranchAddress("Emenin", &Emenin);
	in_chain.SetBranchAddress("Emtrknr", &Emtrknr);
	in_chain.SetBranchAddress("Emtrkph", &Emtrkph);
	in_chain.SetBranchAddress("Emph", &Emph);
	in_chain.SetBranchAddress("Empt", &Empt);
	in_chain.SetBranchAddress("Trk_pca", &Trk_pca);
	in_chain.SetBranchAddress("Emcalene", &Emcalene);
	in_chain.SetBranchAddress("Emtrkp", &Emtrkp);
	in_chain.SetBranchAddress("Emxel", &Emxel);
	in_chain.SetBranchAddress("Emyel", &Emyel);
	in_chain.SetBranchAddress("Emq2el", &Emq2el);
	in_chain.SetBranchAddress("Emxda", &Emxda);
	in_chain.SetBranchAddress("Emyda", &Emyda);
	in_chain.SetBranchAddress("Emq2da", &Emq2da);
	in_chain.SetBranchAddress("Emxda_cell", &Emxda_cell);
	in_chain.SetBranchAddress("Emyda_cell", &Emyda_cell);
	in_chain.SetBranchAddress("Emq2da_cell", &Emq2da_cell);
	in_chain.SetBranchAddress("Emfemc", &Emfemc);
	

	//
	in_chain.SetBranchAddress("Bspt_x",&Bspt_x);
  	in_chain.SetBranchAddress("Bspt_y",&Bspt_y);
	in_chain.SetBranchAddress("Bspt_z",&Bspt_z);
 	in_chain.SetBranchAddress("Bspt_xer",&Bspt_xer);
  	in_chain.SetBranchAddress("Bspt_yer",&Bspt_yer);
	in_chain.SetBranchAddress("Bspt_zer",&Bspt_zer);
	in_chain.SetBranchAddress("Bspt_dxdz",&Bspt_dxdz);
  	in_chain.SetBranchAddress("Bspt_dydz",&Bspt_dydz);



	//Create a new file and in_chain for mapped values (CMS)
	TFile rootfilemapped("/nfs/dust/zeus/group/imelnyk/zeustocms.root","recreate");
	TTree *Events = new TTree("Events", "tree mapped");
	Events->Branch("run", &run,"run/I");
  	Events->Branch("event", &event,"event/I");
	Events->Branch("luminosityBlock", &luminosityBlock,"luminosityBlock/I");

	//
	Events->Branch("PV_x", &PV_x,"PV_x/F");
  	Events->Branch("PV_y", &PV_y,"PV_y/F");
	Events->Branch("PV_z", &PV_z,"PV_z/F");
  	Events->Branch("PV_chi2", &PV_chi2,"PV_chi2/F");

	//
	Events->Branch("nPVtx", &nPVtx,"nPVtx/I");
  	Events->Branch("PVtx_x", PVtx_x,"PVtx_x[nPVtx]/F");
	Events->Branch("PVtx_y", PVtx_y,"PVtx_y[nPVtx]/F");
	Events->Branch("PVtx_z", PVtx_x,"PVtx_z[nPVtx]/F");
	Events->Branch("PVtx_isValid", PVtx_isValid,"PVtx_isValid[nPVtx]/B");
	Events->Branch("PVtx_isGood", PVtx_isGood,"PVtx_isGood[nPVtx]/B");
	Events->Branch("PVtx_isMain", PVtx_isMain,"PVtx_isMain[nPVtx]/B");
	Events->Branch("PVtx_isFake", PVtx_isFake,"PVtx_isFake[nPVtx]/B");
	Events->Branch("PVtx_ntrkfit", PVtx_ntrkfit,"PVtx_ntrkfit[nPVtx]/I");
	Events->Branch("PVtx_chi2", PVtx_chi2,"PVtx_chi2[nPVtx]/I");
	Events->Branch("PVtx_Id", PVtx_Id,"PVtx_Id[nPVtx]/I");
	Events->Branch("PVtx_Rho", PVtx_Rho,"PVtx_Rho[nPVtx]/I");
	
	//Muon data
	Events->Branch("nMuon", &nMuon,"nMuon/I");
  	Events->Branch("Muon_charge", Muon_charge,"Muon_charge[nMuon]/I");
	Events->Branch("Muon_pt", Muon_pt,"Muon_pt[nMuon]/F");
  	Events->Branch("Muon_eta", Muon_eta,"Muon_eta[nMuon]/F");
	Events->Branch("Muon_phi", Muon_phi,"Muon_phi[nMuon]/F");
	Events->Branch("Muon_mass", Muon_mass,"Muon_mass[nMuon]/F");
	Events->Branch("Muon_looseId", Muon_looseId,"Muon_looseId[nMuon]/B");
	Events->Branch("Muon_softId", Muon_softId,"Muon_softId[nMuon]/B");
  	Events->Branch("Muon_mediumId", Muon_mediumId,"Muon_mediumId[nMuon]/B");
	Events->Branch("Muon_tightId", Muon_tightId,"Muon_tightId[nMuon]/B");
	Events->Branch("Muon_trkIdx", Muon_trkIdx,"Muon_trkIdx[nMuon]/I");
	Events->Branch("Muon_ptErr", Muon_ptErr,"Muon_ptErr[nMuon]/F");
  	Events->Branch("Muon_isPFcand", Muon_isPFcand,"Muon_isPFcand[nMuon]/B");
	Events->Branch("Muon_isTracker", Muon_isTracker,"Muon_isTracker[nMuon]/B");
	Events->Branch("Muon_isGlobal", Muon_isGlobal,"Muon_isGlobal[nMuon]/B");
	Events->Branch("Muon_isStandAlone", Muon_isStandAlone,"Muon_isStandAlone[nMuon]/B");
	Events->Branch("Muon_pdgId", Muon_pdgId,"pdgId[nMuon]/I");
  	Events->Branch("Muon_chi2", Muon_chi2,"Muon_chi2[nMuon]/I");
	Events->Branch("Muon_z", Muon_z,"Muon_z[nMuon]/F");
	Events->Branch("Muon_dxy", Muon_dxy,"Muon_dxy[nMuon]/F");
	Events->Branch("Muon_dz", Muon_dz,"Muon_dz[nMuon]/F");
	Events->Branch("Muon_pfreliso03_all", Muon_pfreliso03_all,"Muon_pfreliso03_all[nMuon]/F");
	Events->Branch("Muon_vtxIdx", Muon_vtxIdx,"Muon_vtxIdx[nMuon]/F");
	Events->Branch("Muon_vtxFlag", Muon_vtxFlag,"Muon_vtxFlag[nMuon]/F");
	Events->Branch("Muon_jetIdx", Muon_jetIdx,"Muon_jetIdx[nMuon]/F");

	// Electron data
	Events->Branch("nElectron", &nElectron, "nElectron/I");
	Events->Branch("Electron_SCeta", Electron_SCeta, "Electron_SCeta[nElectron]/F");
	Events->Branch("Electron_charge", Electron_charge, "Electron_charge[nElectron]/I");
	Events->Branch("Electron_convDcot", Electron_convDcot, "Electron_convDcot[nElectron]/F");
	Events->Branch("Electron_convDist", Electron_convDist, "Electron_convDist[nElectron]/F");
	Events->Branch("Electron_convVeto", Electron_convVeto, "Electron_convVeto[nElectron]/B");
	Events->Branch("Electron_convVetoOld", Electron_convVetoOld, "Electron_convVetoOld[nElectron]/B");
	Events->Branch("Electron_cutBased", Electron_cutBased, "Electron_cutBased[nElectron]/I");
	Events->Branch("Electron_deltaEtaSC", Electron_deltaEtaSC, "Electron_deltaEtaSC[nElectron]/F");
	Events->Branch("Electron_deltaEtaSCtr", Electron_deltaEtaSCtr, "Electron_deltaEtaSCtr[nElectron]/F");
	Events->Branch("Electron_deltaPhiSC", Electron_deltaPhiSC, "Electron_deltaPhiSC[nElectron]/F");
	Events->Branch("Electron_deltaPhiSCtr", Electron_deltaPhiSCtr, "Electron_deltaPhiSCtr[nElectron]/F");
	Events->Branch("Electron_dr03EcalRecHitSumEt", Electron_dr03EcalRecHitSumEt, "Electron_dr03EcalRecHitSumEt[nElectron]/F");
	Events->Branch("Electron_dr03EcalRecHitSumEtOld", Electron_dr03EcalRecHitSumEtOld, "Electron_dr03EcalRecHitSumEtOld[nElectron]/F");
	Events->Branch("Electron_dr03HcalDepth1TowerSumEt", Electron_dr03HcalDepth1TowerSumEt, "Electron_dr03HcalDepth1TowerSumEt[nElectron]/F");
	Events->Branch("Electron_dr03HcalDepth1TowerSumEtOld", Electron_dr03HcalDepth1TowerSumEtOld, "Electron_dr03HcalDepth1TowerSumEtOld[nElectron]/F");
	Events->Branch("Electron_dr03HcalTowerSumEt", Electron_dr03HcalTowerSumEt, "Electron_dr03HcalTowerSumEt[nElectron]/F");
	Events->Branch("Electron_dr03TkSumPt", Electron_dr03TkSumPt, "Electron_dr03TkSumPt[nElectron]/F");
	Events->Branch("Electron_dr03TkSumPtOld", Electron_dr03TkSumPtOld, "Electron_dr03TkSumPtOld[nElectron]/F");
	Events->Branch("Electron_dxy", Electron_dxy, "Electron_dxy[nElectron]/F");
	Events->Branch("Electron_dxyErr", Electron_dxyErr, "Electron_dxyErr[nElectron]/F");
	Events->Branch("Electron_dz", Electron_dz, "Electron_dz[nElectron]/F");	
	Events->Branch("Electron_dzErr", Electron_dzErr, "Electron_dzErr[nElectron]/F");
	Events->Branch("Electron_eInvMinusPInv", Electron_eInvMinusPInv, "Electron_eInvMinusPInv[nElectron]/F");
	Events->Branch("Electron_eInvMinusPInvOld", Electron_eInvMinusPInvOld, "Electron_eInvMinusPInvOld[nElectron]/F");
	Events->Branch("Electron_eta", Electron_eta, "Electron_eta[nElectron]/F");
	Events->Branch("Electron_genPartIdx", Electron_genPartIdx, "Electron_genPartIdx[nElectron]/I");
	Events->Branch("Electron_hoe", Electron_hoe, "Electron_hoe[nElectron]/F");
	Events->Branch("Electron_ip3d", Electron_ip3d, "Electron_ip3d[nElectron]/F");
	Events->Branch("Electron_isEB", Electron_isEB, "Electron_isEB[nElectron]/B");
	Events->Branch("Electron_isEE", Electron_isEE, "Electron_isEE[nElectron]/B");
	Events->Branch("Electron_isNano", Electron_isNano, "Electron_isNano[nElectron]/B");
	Events->Branch("Electron_isPFcand", Electron_isPFcand, "Electron_isPFcand[nElectron]/B");
	Events->Branch("Electron_lostHits", Electron_lostHits, "Electron_lostHits[nElectron]/C");
	Events->Branch("Electron_mass", Electron_mass, "Electron_mass[nElectron]/F");
	Events->Branch("Electron_nNano", Electron_nNano, "Electron_nNano[nElectron]/I");
	Events->Branch("Electron_pfRelIso03_all", Electron_pfRelIso03_all, "Electron_pfRelIso03_all[nElectron]/F");
	Events->Branch("Electron_pfRelIso03_chg", Electron_pfRelIso03_chg, "Electron_pfRelIso03_chg[nElectron]/F");
	Events->Branch("Electron_phi", Electron_phi, "Electron_phi[nElectron]/F");
	Events->Branch("Electron_pt", Electron_pt, "Electron_pt[nElectron]/F");
	Events->Branch("Electron_sieie", Electron_sieie, "Electron_sieie[nElectron]/F");
	Events->Branch("Electron_sieieR1", Electron_sieieR1, "Electron_sieieR1[nElectron]/F");
	Events->Branch("Electron_simId", Electron_simId, "Electron_simId[nElectron]/I");
	Events->Branch("Electron_sip3d", Electron_sip3d, "Electron_sip3d[nElectron]/F");
	Events->Branch("Electron_tightCharge", Electron_tightCharge, "Electron_tightCharge[nElectron]/I");
	Events->Branch("Electron_vtxIdx", Electron_vtxIdx, "Electron_vtxIdx[nElectron]/I");
	Events->Branch("Electron_x", Electron_x, "Electron_x[nElectron]/F");
	Events->Branch("Electron_y", Electron_y, "Electron_y[nElectron]/F");
	Events->Branch("Electron_z", Electron_z, "Electron_z[nElectron]/F");
	Events->Branch("Emprob_ZEUS", Emprob_ZEUS, "Emprob_ZEUS[nElectron]/F");

	Events->Branch("Emxel_ZEUS", Emxel_ZEUS, "Emxel_ZEUS[nElectron]/F");
	Events->Branch("Emyel_ZEUS", Emyel_ZEUS, "Emyel_ZEUS[nElectron]/F");
	Events->Branch("Emq2el_ZEUS", Emq2el_ZEUS, "Emq2el_ZEUS[nElectron]/F");
	Events->Branch("Emxda_ZEUS", Emxda_ZEUS, "Emxda_ZEUS[nElectron]/F");
	Events->Branch("Emyda_ZEUS", Emyda_ZEUS, "Emyda_ZEUS[nElectron]/F");
	Events->Branch("Emq2da_ZEUS", Emq2da_ZEUS, "Emq2da_ZEUS[nElectron]/F");
	Events->Branch("Emxda_cell_ZEUS", Emxda_cell_ZEUS, "Emxda_cell_ZEUS[nElectron]/F");
	Events->Branch("Emyda_cell_ZEUS", Emyda_cell_ZEUS, "Emyda_cell_ZEUS[nElectron]/F");
	Events->Branch("Emq2da_cell_ZEUS", Emq2da_cell_ZEUS, "Emq2da_cell_ZEUS[nElectron]/F");


	//
	Events->Branch("Bsp_x", &Bsp_x);
  	Events->Branch("Bsp_y", &Bsp_y);
	Events->Branch("Bsp_z", &Bsp_z);
 	Events->Branch("Bsp_widthx", &Bsp_widthx);
  	Events->Branch("Bsp_widthy", &Bsp_widthy);
	Events->Branch("Bsp_sigmaz", &Bsp_sigmaz);
	Events->Branch("Bsp_dxdz", &Bsp_dxdz);
  	Events->Branch("Bsp_dydz", &Bsp_dydz);

	//Get length of in_chain
	Int_t nin_chain = in_chain.GetEntries();

	//Our guesses for good Muqual parameters
	Int_t Muqual_looseId = 0;
	Int_t Muqual_softId = 3;
	Int_t Muqual_mediumId = 4;
	Int_t Muqual_tightId = 6;
	
	//Put the corresponding data into the new in_chain
	// Event loop
	for (int j = 0;  j < nin_chain; j++) {

		//Flag
		if (j % 10000 == 0) cout << "Event number: " << j << " / " << nin_chain << "  ---  " << float(j)/float(nin_chain)*100 << "%" << endl;
		

		in_chain.GetEntry(j);
		
		run = Runnr;
		event = Eventnr;
		luminosityBlock = 0;
		PV_x = Xvtx;
		PV_y = Yvtx;
		PV_z = Zvtx;
		PV_chi2 = Chivtx;

		nPVtx = 1;
		PVtx_x[0] = Xvtx;
		PVtx_y[0] = Yvtx;
		PVtx_z[0] = Zvtx;
		if(Chivtx>=0) {
			PVtx_isValid[0] = true;
			PVtx_isFake[0] = false;
		}
		else {
			PVtx_isValid[0] = false;
			PVtx_isFake[0] = true;
		}
		PVtx_isGood[0] = true;//some cuts!!!
		PVtx_isMain[0] = true;
		PVtx_ntrkfit[0] = Ntrkvtx;
		PVtx_chi2[0] = Chivtx;
		PVtx_Id[0] = 1;
		PVtx_Rho[0] = sqrt((Xvtx-Bspt_x)*(Xvtx-Bspt_x) + (Yvtx-Bspt_y)*(Yvtx-Bspt_y));

		//Give a warning when Nmu is greater than our chosen Nmax from teh beginning
		if((Nmu > Nmax) || (Emncand > Nmax) || (Trk_ntracks > trackNMax)) { // 
			cout << "WARNING! Nmu is greater than your selected maximal length!" << endl;
		}

		nMuon = Nmu;
		//Loop over the Muon arrays
		for (int i=0; i<nMuon; i++) {
			Muon_charge[i] = Mucharge[i];
			Muon_pt[i] = Mupt[i];
			Muon_eta[i] = -log(tan(Muth[i]/2));
			Muon_phi[i] = Muph[i];
			Muon_mass[i] = 0.105658;
			
			//Check Muqual
			//Set them first all to false (default)
			Muon_tightId[i] = false;
			Muon_mediumId[i] = false;
			Muon_softId[i] = false;
			Muon_looseId[i] = false;
			
			if(Muqual[i] == Muqual_tightId) {
				Muon_tightId[i] = true;
			}
			if(Muqual[i] > Muqual_mediumId) {
				Muon_mediumId[i] = true;
			}
			if(Muqual[i] > Muqual_softId) {
				Muon_softId[i] = true;
			}
			if(Muqual[i] > Muqual_looseId) {
				Muon_looseId[i] = true;
			}

			Muon_trkIdx[i] = Mutrid[i];
			Muon_ptErr[i] = Muperr[i]*sin(Muth[i]);

			if(Muzufid[i] > 0) {
				Muon_isPFcand[i] = true;
			}
			else {
				Muon_isPFcand[i] = false;
			}

			Muon_isTracker[i] = true;
			
			if(Muqual[i] > 0) {
				Muon_isGlobal[i] = true;
			}
			else {
				Muon_isGlobal[i] = false;
			}

			if(Muqual[i] == 6) {
				Muon_isStandAlone[i] = true;
			}
			else {
				Muon_isStandAlone[i] = false;
			}

			Muon_pdgId[i] = -13*Mucharge[i];
			Muon_chi2[i] = Muchid[i];
			Muon_z[i] = Muz[i];
			Muon_dxy[i] = Mudxy[i];
			Muon_dz[i] = Mudz[i];
			Muon_pfreliso03_all[i] = 1; //needs to be fixed!!!
			Muon_vtxIdx[i] = Muvtxid[i];
			Muon_vtxFlag[i] = Muvtxfl[i];
			Muon_jetIdx[i] = Mujetid_a[i];

			//calculate and fill Dimu structure
		}
		nElectron = Emncand;

		for (int iElectron = 0; iElectron < nElectron; iElectron++) {

			Float_t Electron_SCeta_x = Emcalpos[iElectron][0] - PV_x;
			Float_t Electron_SCeta_y = Emcalpos[iElectron][1] - PV_y;
			Float_t Electron_SCeta_z = Emcalpos[iElectron][2] - PV_z;
			Float_t R = sqrt(pow(Electron_SCeta_x, 2) + pow(Electron_SCeta_y, 2));
			Float_t Theta = atan2(R, Electron_SCeta_z);
			Electron_SCeta[iElectron] = -log(tan(Theta / 2));			
			Electron_convDcot[iElectron] = -999.; // needed further investigation
			Electron_convDist[iElectron] = -999.; // needed further invessigation
			Electron_convVeto[iElectron] = false; // needed further investigation
			Electron_convVetoOld[iElectron] = false; // needed further investigation

			if (Emtrknr[iElectron] == 0) { // electron doesn't have assosiation with trk

				Electron_eta[iElectron] = Electron_SCeta[iElectron]; // take from cal
				Electron_dr03TkSumPt[iElectron] = -999.;
				Electron_dr03TkSumPtOld[iElectron] = -999.;
				Electron_eInvMinusPInv[iElectron] = -999.;
				Electron_eInvMinusPInvOld[iElectron] = -999.;
				if (Emph[iElectron] > PiConstant) {
					Electron_phi[iElectron] = Emph[iElectron] - DoublePiConstant;
				}
				else if (Emph[iElectron] < -PiConstant) {
					Electron_phi[iElectron] = DoublePiConstant - abs(Emph[iElectron]);
				}
				else {
					Electron_phi[iElectron] = Emph[iElectron]; // take from cal
				}
				Electron_pfRelIso03_chg[iElectron] = -999.;
				Electron_charge[iElectron] = 0;

			}
			else {

				Int_t emTrackID = Emtrknr[iElectron];
				// ********************** Loop over all tracks *************************************************
				// ********************** Check track id and match it with electron track **********************
				for (int nTrk = 0; nTrk < Trk_ntracks; nTrk++) {

					if (Trk_id[nTrk] == emTrackID) {
						Electron_x[iElectron] = Trk_pca[nTrk][0];
						Electron_y[iElectron] = Trk_pca[nTrk][1];
						Electron_z[iElectron] = Trk_pca[nTrk][2];
					}

				}

				Electron_charge[iElectron] = Emtrkq[iElectron]; // needed modification: put 0 if electron is not from track (e.g. from calorimetr)
				Electron_eta[iElectron] = -log(tan(Emtrkth[iElectron] / 2));
				Electron_dr03TkSumPt[iElectron] = Emetneartrk[iElectron][1];
				Electron_dr03TkSumPtOld[iElectron] = Emetneartrk[iElectron][1];
				Electron_eInvMinusPInv[iElectron] = (1/Emcalene[iElectron]) - (1/Emtrkp[iElectron]);
				Electron_eInvMinusPInvOld[iElectron] = (1/Emcalene[iElectron]) - (1/Emtrkp[iElectron]); // ????
				if (Emtrkph[iElectron] > PiConstant) {
					Electron_phi[iElectron] = Emtrkph[iElectron] - DoublePiConstant;
				}
				else if (Emtrkph[iElectron] < -PiConstant) {
					Electron_phi[iElectron] = DoublePiConstant - abs(Emtrkph[iElectron]);
				}
				else {
					Electron_phi[iElectron] = Emtrkph[iElectron]; // take from cal
				}
				Electron_pfRelIso03_chg[iElectron] = Emetneartrk[iElectron][1] / Empt[iElectron];
			}


			if (Emprob[iElectron] < 0.001) {
				Electron_cutBased[iElectron] = 0;
			} // very bad
			else {
				if (Emtrknr[iElectron] == 0) {
					Electron_cutBased[iElectron] = 2;
				}
				else {
					Electron_cutBased[iElectron] = 4; // cut: Electron_cutBased == 4 -- from track
				}
			}
		
			 // very good // cut: Electron_cutBased == 4
			 // needed further investigation
			Electron_deltaEtaSC[iElectron] = Electron_SCeta[iElectron] - Electron_eta[iElectron]; // needed further investigation
			Electron_deltaEtaSCtr[iElectron] = -999.; // needed further investigation
			Electron_deltaPhiSC[iElectron] = -999.; // needed further investigation
			Electron_deltaPhiSCtr[iElectron] = -999.; // needed further investigation
			Electron_dr03EcalRecHitSumEt[iElectron] = -999.; // needed further investigation
			Electron_dr03EcalRecHitSumEtOld[iElectron] = -999.; // needed further investigation
			Electron_dr03HcalDepth1TowerSumEt[iElectron] = -999.; // needed further investigation
			Electron_dr03HcalDepth1TowerSumEtOld[iElectron] = -999.; // needed further investigation
			Electron_dr03HcalTowerSumEt[iElectron] = -999.; // needed further investigation
			Electron_dxy[iElectron] = Emdcabeam[iElectron];
			Electron_dxyErr[iElectron] = 0.2;
			Electron_dz[iElectron] = sqrt(pow(Emdca[iElectron], 2) - pow(Emdcabeam[iElectron], 2));
			Electron_dzErr[iElectron] = 0.4;
			Electron_genPartIdx[iElectron] = -1; // needed further investigation
			Electron_hoe[iElectron] = (1/Emfemc[iElectron]) - 1; 
			Electron_ip3d[iElectron] = Emdca[iElectron];
			Electron_isEB[iElectron] = bool((Empos[iElectron][2] < Float_t(200.)) && (Empos[iElectron][2] > Float_t(-120.)));
			Electron_isEE[iElectron] = bool(!((Empos[iElectron][2] < Float_t(200.)) && (Empos[iElectron][2] > Float_t(-120.))));
			Electron_isNano[iElectron] = true;
			Electron_isPFcand[iElectron] = true; // needed further investigation
			Electron_lostHits[iElectron] = '0'; // needed further investigation
			Electron_mass[iElectron] = 0.000511; // in GeV taken from Particle Data Group
			Electron_nNano[iElectron] = Emncand;
			Electron_pfRelIso03_all[iElectron] = Emenin[iElectron] / Empt[iElectron];
			Electron_pt[iElectron] = Empt[iElectron];
			Electron_sieie[iElectron] = -999.; // needed furhter investigation
			Electron_sieieR1[iElectron] = -999.; // needed further investigation
			Electron_simId[iElectron] = -1;
			Electron_sip3d[iElectron] = -999.; // needed further investigation
			Electron_tightCharge[iElectron] = 1; // needed further investigation
			Electron_vtxIdx[iElectron] = -1;
			Emprob_ZEUS[iElectron] = Emprob[iElectron];
			Emxel_ZEUS[iElectron] = Emxel[iElectron];
			Emyel_ZEUS[iElectron] = Emyel[iElectron];
			Emq2el_ZEUS[iElectron] = Emq2el[iElectron];
			Emxda_ZEUS[iElectron] = Emxda[iElectron];
			Emyda_ZEUS[iElectron] = Emyda[iElectron];
			Emq2da_ZEUS[iElectron] = Emq2da[iElectron];
			Emxda_cell_ZEUS[iElectron] = Emxda_cell[iElectron];
			Emyda_cell_ZEUS[iElectron] = Emyda_cell[iElectron];
			Emq2da_cell_ZEUS[iElectron] = Emq2da_cell[iElectron];


			// std::cout << "   ***** Inside electron loop: " << iElectron << std::endl;
		}

		Bsp_x = Bspt_x;
		Bsp_y = Bspt_y;
		Bsp_z = Bspt_z;
		Bsp_widthx = Bspt_xer;
		Bsp_widthy = Bspt_yer;
		Bsp_sigmaz = Bspt_zer;
		Bsp_dxdz = Bspt_dxdz;
		Bsp_dydz = Bspt_dydz;
	
		Events->Fill();
	}

Events->Write();
printf("=========================== End =====================================\n");




///////////////////////////////////Plot some histograms//////////////////////////////////////////////

/*
	//Create histograms
  	TH1D *hMuon_charge = new TH1D("hMuon_charge", "Muon charge;Charge;Counts", 20, -2, 2);	hMuon_charge->SetDirectory(0);
  	TH1D *hMuon_pt = new TH1D("hMuon_pt", "Muon p_{T};Transverse momentum; Counts", 50, 0, 20);	hMuon_pt->SetDirectory(0);
  	TH1D *hMuon_eta = new TH1D("hMuon_eta", "Muon #eta;#eta;Counts", 100, -7, 7);	hMuon_eta->SetDirectory(0);
	TH1D *hMuon_phi = new TH1D("hMuon_phi", "Muon #phi;#phi;Counts", 100, 0, 7);	hMuon_phi->SetDirectory(0);

	//Loop over the entries of the in_chain
	for(int i=0; i<nin_chain; i++) {
		Events->GetEntry(i);
		
		//Loop over each Muon Array
		for(int j=0; j<nMuon; j++) {
 			hMuon_charge->Fill(Muon_charge[j]);
			hMuon_pt->Fill(Muon_pt[j]);
    			hMuon_eta->Fill(Muon_eta[j]);
    			hMuon_phi->Fill(Muon_phi[j]);
  		}
	}

	//Plot them together on a Canvas
  	TCanvas *cMu1 = new TCanvas("cMu1","Mapping from ZEUS to CMS", 600, 400);
  	cMu1->Divide(1,1);
	TCanvas *cMu2 = new TCanvas("cMu2","Mapping from ZEUS to CMS", 600, 400);
  	cMu2->Divide(1,1);
	TCanvas *cMu3 = new TCanvas("cMu3","Mapping from ZEUS to CMS", 600, 400);
  	cMu3->Divide(1,1);

  	//cMu1->cd(1);
  	//hMuon_charge->Draw();

  	cMu1->cd(1);
 	hMuon_pt->Draw();

  	cMu2->cd(1);
  	hMuon_eta->Draw();

  	cMu3->cd(1);
  	hMuon_phi->Draw();

  	cMu1->Draw();
	cMu1->SaveAs("Muon_pt_zeustocms.pdf");
	cMu2->Draw();
	cMu2->SaveAs("Muon_eta_zeustocms.pdf");
	cMu3->Draw();
	cMu3->SaveAs("Muon_phi_zeustocms.pdf");
*/

//in_chain.Close();
//rootfilemapped.Close();
}

int main() {

	zeustocms();

	return 0;
}
