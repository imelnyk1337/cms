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
#include "TTree.h" 
#include "TMath.h"
#include "TFile.h"
// #include "Math/VectorUtil_cint.h"

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
// #include "Math/VectorUtil.h"
// #include "Math/GenVector/VectorUtil.h"

#define makeMuons true
#define makeElectrons true
#define makeJets true


using namespace std;


void zeustocms_t() {

	std::cout << "=========================== Start conversion =====================================" << std::endl;
	// Maximum of Nmu. We need this to declare our arrays before we know each Nmu.
	const Int_t Nmax = 20;

	// ********************** Common variables **********************

	Int_t Runnr, Eventnr, Ntrkvtx, Trk_ntracks, Trk_id[Nmax];
	Float_t Xvtx, Yvtx, Zvtx, Chivtx;

	// ==================================================================================================================================== 
	// ========================================================== ZEUS variables ==========================================================
	// ====================================================================================================================================

	// ********************** Muons variables *************************
	Int_t Nmu;
	Int_t Mucharge[Nmax], Muqual[Nmax], Mutrid[Nmax], Muzufid[Nmax], Muchid[Nmax], Muvtxid[Nmax], Muvtxfl[Nmax], Mujetid_a[Nmax];
	Float_t Mupt[Nmax], Muth[Nmax], Muph[Nmax], Muperr[Nmax], Mup[Nmax][3], Muz[Nmax], Mudxy[Nmax], Mudz[Nmax], Muisol[Nmax][10];

	Float_t Bspt_x, Bspt_y, Bspt_z, Bspt_xer, Bspt_yer, Bspt_zer, Bspt_dxdz, Bspt_dydz;	


	// ********************** Electron variables **********************
	Int_t Emncand, Emtrkq[Nmax], Emtrknr[Nmax]; // electric charge
	Float_t Emcalpos[Nmax][3], Emetneartrk[Nmax][2], Emdcabeam[Nmax], Emdca[Nmax], Emtrkth[Nmax], Emenin[Nmax],
							Emtrkph[Nmax], Empt[Nmax], Trk_pca[Nmax][3], Emcalene[Nmax], Emtrkp[Nmax];
	
	// ====================================================================================================================================
	// ========================================================== CMS variables ===========================================================
	// ====================================================================================================================================


	// ********************** Common variables **********************

	Int_t run, event, luminosityBlock, nPVtx, nMuon, nElec;

	nPVtx = 100;
	nMuon = 100;

	Float_t PV_x, PV_y, PV_z, PV_chi2;

	Float_t PVtx_x[nPVtx], PVtx_y[nPVtx], PVtx_z[nPVtx], PVtx_chi2[nPVtx], PVtx_Rho[nPVtx];
	
	Bool_t PVtx_isValid[nPVtx], PVtx_isGood[nPVtx], PVtx_isMain[nPVtx], PVtx_isFake[nPVtx];
	
	Int_t PVtx_ntrkfit[nPVtx], PVtx_Id[nPVtx];

	// ********************** Muons variables *************************

	Int_t Muon_charge[nMuon], Muon_trkIdx[nMuon], Muon_pdgId[nMuon];

	Bool_t Muon_looseId[nMuon], Muon_softId[nMuon], Muon_mediumId[nMuon], Muon_tightId[nMuon], Muon_isPFcand[nMuon], Muon_isTracker[nMuon], Muon_isGlobal[nMuon], Muon_isStandAlone[nMuon];					

	Float_t Muon_pt[nMuon], Muon_eta[nMuon], Muon_phi[nMuon], Muon_mass[nMuon], Muon_ptErr[nMuon], Muon_z[nMuon], Muon_dxy[nMuon], Muon_dz[nMuon], Muon_pfreliso03_all[nMuon], Muon_vtxIdx[nMuon], Muon_vtxFlag[nMuon], Muon_jetIdx[nMuon], Muon_chi2[nMuon];

	Float_t Bsp_x, Bsp_y, Bsp_z, Bsp_widthx, Bsp_widthy, Bsp_sigmaz, Bsp_dxdz, Bsp_dydz;


	// ********************** Electron variables **********************
	UInt_t nElectron, Electron_nNano[nElectron];
	nElectron = 200;
	Int_t Electron_charge[nElectron], Electron_cutBased[nElectron], Electron_simId[nElectron], Electron_tightCharge[nElectron], Electron_vtxIdx[nElectron];
	UChar_t Electron_lostHits[nElectron];
	Bool_t Electron_convVeto[nElectron], Electron_convVetoOld[nElectron], Electron_isEB[nElectron], Electron_isEE[nElectron], Electron_isNano[nElectron], Electron_isPFcand[nElectron];  

	Float_t Electron_SCeta[nElectron], Electron_convDcot[nElectron], Electron_convDist[nElectron],
			Electron_deltaEtaSC[nElectron], Electron_deltaEtaSCtr[nElectron],
			Electron_deltaPhiSC[nElectron], Electron_deltaPhiSCtr[nElectron], Electron_dr03EcalRecHitSumEt[nElectron],
			Electron_dr03EcalRecHitSumEtOld[nElectron], Electron_dr03HcalDepth1TowerSumEt[nElectron], Electron_dr03HcalDepth1TowerSumEtOld[nElectron],
			Electron_dr03HcalTowerSumEt[nElectron], Electron_dr03TkSumPt[nElectron], Electron_dr03TkSumPtOld[nElectron], Electron_dxy[nElectron],
			Electron_dxyErr[nElectron], Electron_dz[nElectron], Electron_dzErr[nElectron], Electron_eInvMinusPInv[nElectron],
			Electron_eInvMinusPInvOld[nElectron], Electron_eta[nElectron], Electron_genPartIdx[nElectron], Electron_hoe[nElectron],
			Electron_ip3d[nElectron], Electron_mass[nElectron], Electron_pfRelIso03_all[nElectron], Electron_pfRelIso03_chg[nElectron],
			Electron_ph[nElectron], Electron_pt[nElectron], Electron_sieie[nElectron],
			Electron_sieieR1[nElectron], Electron_sip3d[nElectron], Electron_x[nElectron], Electron_y[nElectron], Electron_z[nElectron];




	// Get the demanded file
	TChain inputChain("orange");
	inputChain.Add("/pnfs/desy.de/dphep/online/zeus/z/ntup/07p/v08b/data/root/data_07p_61792_61792_01.root");	//good ntuple file to check code fast: data_07p_61792_61792_01.root

	// ********************** Linking common data **********************************

  	inputChain.SetBranchAddress("Runnr", &Runnr);
  	inputChain.SetBranchAddress("Eventnr", &Eventnr);
	inputChain.SetBranchAddress("Trk_ntracks", &Trk_ntracks);
	inputChain.SetBranchAddress("Trk_id", &Trk_id);

	inputChain.SetBranchAddress("Xvtx", &Xvtx);
  	inputChain.SetBranchAddress("Yvtx", &Yvtx);
	inputChain.SetBranchAddress("Zvtx", &Zvtx);
  	inputChain.SetBranchAddress("Chivtx", &Chivtx);
	inputChain.SetBranchAddress("Ntrkvtx", &Ntrkvtx);
	
	// ********************** Linking muons data from ZEUS *************************
	inputChain.SetBranchAddress("Nmu", &Nmu);
  	inputChain.SetBranchAddress("Mucharge", &Mucharge);
	inputChain.SetBranchAddress("Mupt", &Mupt);
  	inputChain.SetBranchAddress("Muth", &Muth);
	inputChain.SetBranchAddress("Muph", &Muph);
	inputChain.SetBranchAddress("Muqual", &Muqual);
  	inputChain.SetBranchAddress("Muperr", &Muperr);
	inputChain.SetBranchAddress("Mup", &Mup);
  	inputChain.SetBranchAddress("Muzufid", &Muzufid);
	inputChain.SetBranchAddress("Muchid", &Muchid);
	inputChain.SetBranchAddress("Muz", &Muz);
  	inputChain.SetBranchAddress("Mudxy", &Mudxy);
	inputChain.SetBranchAddress("Mudz", &Mudz);
  	inputChain.SetBranchAddress("Muisol", &Muisol);
	inputChain.SetBranchAddress("Muvtxid", &Muvtxid);
	inputChain.SetBranchAddress("Muvtxfl", &Muvtxfl);
  	inputChain.SetBranchAddress("Mujetid_a", &Mujetid_a);

	//
	inputChain.SetBranchAddress("Bspt_x", &Bspt_x);
  	inputChain.SetBranchAddress("Bspt_y", &Bspt_y);
	inputChain.SetBranchAddress("Bspt_z", &Bspt_z);
 	inputChain.SetBranchAddress("Bspt_xer", &Bspt_xer);
  	inputChain.SetBranchAddress("Bspt_yer", &Bspt_yer);
	inputChain.SetBranchAddress("Bspt_zer", &Bspt_zer);
	inputChain.SetBranchAddress("Bspt_dxdz", &Bspt_dxdz);
  	inputChain.SetBranchAddress("Bspt_dydz", &Bspt_dydz);

	// ********************** Linking electron data from ZEUS **********************
	inputChain.SetBranchAddress("Emncand", &Emncand);
	inputChain.SetBranchAddress("Emcalpos", &Emcalpos);
	inputChain.SetBranchAddress("Emetneartrk", &Emetneartrk);
	inputChain.SetBranchAddress("Emdcabeam", &Emdcabeam);
	inputChain.SetBranchAddress("Emdca", &Emdca);
	inputChain.SetBranchAddress("Emtrkth", &Emtrkth);
	inputChain.SetBranchAddress("Emenin", &Emenin);
	inputChain.SetBranchAddress("Emtrknr", &Emtrknr);
	inputChain.SetBranchAddress("Emtrkph", &Emtrkph);
	inputChain.SetBranchAddress("Empt", &Empt);
	inputChain.SetBranchAddress("Trk_pca", &Trk_pca);
	inputChain.SetBranchAddress("Emcalene", &Emcalene);
	inputChain.SetBranchAddress("Emtrkp", &Emtrkp);

	// ==================================================================================================================================== 
	// ==================================================================================================================================== 
	// ==================================================================================================================================== 


	// Create a new file and inputChain for mapped values (CMS)
	TFile outputFile("/nfs/dust/zeus/group/imelnyk/zeustocms_t.root", "recreate");
	TTree* Events = new TTree("Events", "tree mapped");
	Events->Branch("run", &run, "run/I");
  	Events->Branch("event", &event, "event/I");
	Events->Branch("luminosityBlock", &luminosityBlock, "luminosityBlock/I");

	//
	Events->Branch("PV_x", &PV_x, "PV_x/F");
  	Events->Branch("PV_y", &PV_y, "PV_y/F");
	Events->Branch("PV_z", &PV_z, "PV_z/F");
  	Events->Branch("PV_chi2", &PV_chi2, "PV_chi2/F");

	//
	Events->Branch("nPVtx", &nPVtx, "nPVtx/I");
  	Events->Branch("PVtx_x", PVtx_x, "PVtx_x[nPVtx]/F");
	Events->Branch("PVtx_y", PVtx_y, "PVtx_y[nPVtx]/F");
	Events->Branch("PVtx_z", PVtx_x, "PVtx_z[nPVtx]/F");
	Events->Branch("PVtx_isValid", PVtx_isValid, "PVtx_isValid[nPVtx]/B");
	Events->Branch("PVtx_isGood", PVtx_isGood, "PVtx_isGood[nPVtx]/B");
	Events->Branch("PVtx_isMain", PVtx_isMain, "PVtx_isMain[nPVtx]/B");
	Events->Branch("PVtx_isFake", PVtx_isFake, "PVtx_isFake[nPVtx]/B");
	Events->Branch("PVtx_ntrkfit", PVtx_ntrkfit, "PVtx_ntrkfit[nPVtx]/I");
	Events->Branch("PVtx_chi2", PVtx_chi2, "PVtx_chi2[nPVtx]/I");
	Events->Branch("PVtx_Id", PVtx_Id, "PVtx_Id[nPVtx]/I");
	Events->Branch("PVtx_Rho", PVtx_Rho, "PVtx_Rho[nPVtx]/I");
	
	// ********************** Muon data CMS **********************
	Events->Branch("nMuon", &nMuon, "nMuon/I");
  	Events->Branch("Muon_charge", Muon_charge, "Muon_charge[nMuon]/I");
	Events->Branch("Muon_pt", Muon_pt, "Muon_pt[nMuon]/F");
  	Events->Branch("Muon_eta", Muon_eta, "Muon_eta[nMuon]/F");
	Events->Branch("Muon_phi", Muon_phi, "Muon_phi[nMuon]/F");
	Events->Branch("Muon_mass", Muon_mass, "Muon_mass[nMuon]/F");
	Events->Branch("Muon_looseId", Muon_looseId, "Muon_looseId[nMuon]/B");
	Events->Branch("Muon_softId", Muon_softId, "Muon_softId[nMuon]/B");
  	Events->Branch("Muon_mediumId", Muon_mediumId, "Muon_mediumId[nMuon]/B");
	Events->Branch("Muon_tightId", Muon_tightId, "Muon_tightId[nMuon]/B");
	Events->Branch("Muon_trkIdx", Muon_trkIdx, "Muon_trkIdx[nMuon]/I");
	Events->Branch("Muon_ptErr", Muon_ptErr, "Muon_ptErr[nMuon]/F");
  	Events->Branch("Muon_isPFcand", Muon_isPFcand, "Muon_isPFcand[nMuon]/B");
	Events->Branch("Muon_isTracker", Muon_isTracker, "Muon_isTracker[nMuon]/B");
	Events->Branch("Muon_isGlobal", Muon_isGlobal, "Muon_isGlobal[nMuon]/B");
	Events->Branch("Muon_isStandAlone", Muon_isStandAlone, "Muon_isStandAlone[nMuon]/B");
	Events->Branch("Muon_pdgId", Muon_pdgId, "pdgId[nMuon]/I");
  	Events->Branch("Muon_chi2", Muon_chi2, "Muon_chi2[nMuon]/I");
	Events->Branch("Muon_z", Muon_z, "Muon_z[nMuon]/F");
	Events->Branch("Muon_dxy", Muon_dxy, "Muon_dxy[nMuon]/F");
	Events->Branch("Muon_dz", Muon_dz, "Muon_dz[nMuon]/F");
	Events->Branch("Muon_pfreliso03_all", Muon_pfreliso03_all, "Muon_pfreliso03_all[nMuon]/F");
	Events->Branch("Muon_vtxIdx", Muon_vtxIdx, "Muon_vtxIdx[nMuon]/F");
	Events->Branch("Muon_vtxFlag", Muon_vtxFlag, "Muon_vtxFlag[nMuon]/F");
	Events->Branch("Muon_jetIdx", Muon_jetIdx, "Muon_jetIdx[nMuon]/F");


	// ********************** Electron data CMS **********************
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
	Events->Branch("Electron_ph", Electron_ph, "Electron_ph[nElectron]/F");
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

	//
	Events->Branch("Bsp_x", &Bsp_x);
  	Events->Branch("Bsp_y", &Bsp_y);
	Events->Branch("Bsp_z", &Bsp_z);
 	Events->Branch("Bsp_widthx", &Bsp_widthx);
  	Events->Branch("Bsp_widthy", &Bsp_widthy);
	Events->Branch("Bsp_sigmaz", &Bsp_sigmaz);
	Events->Branch("Bsp_dxdz", &Bsp_dxdz);
  	Events->Branch("Bsp_dydz", &Bsp_dydz);


	// ==================================================================================================================================== 
	// ==================================================================================================================================== 
	// ==================================================================================================================================== 

	// ====================================================================================================================================
	// ========================================================== Conversion ==============================================================
	// ====================================================================================================================================

	// Get length of inputChain
	Int_t nEntriesInput = inputChain.GetEntries();
	Int_t testSize = 100000;

	// Our guesses for good Muqual parameters
	Int_t Muqual_looseId = 0;
	Int_t Muqual_softId = 3;
	Int_t Muqual_mediumId = 4;
	Int_t Muqual_tightId = 6;
	
	// Put the corresponding data into the new inputChain

	// ********************** Event loop **********************
	for (int iEvt = 0; iEvt < 100; iEvt++) {

		// Flag
		if (iEvt % 10000 == 0) {
			std::cout << "Event number: " << iEvt << " / " << nEntriesInput << "  ---  " << float(iEvt) / float(nEntriesInput) * 100 << "%" << std::endl;
		}
		
		inputChain.GetEntry(iEvt);
		std::cout << "--- Inside event loop: " << iEvt + 1 << " out of 50" << std::endl;
		
		// ********************** Common variables **********************
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

		if (Chivtx >= 0) {
			PVtx_isValid[0] = true;
			PVtx_isFake[0] = false;
		}
		else {
			PVtx_isValid[0] = false;
			PVtx_isFake[0] = true;
		}
		PVtx_isGood[0] = true; // some cuts!!!
		PVtx_isMain[0] = true;
		PVtx_ntrkfit[0] = Ntrkvtx;
		PVtx_chi2[0] = Chivtx;
		PVtx_Id[0] = 1;
		PVtx_Rho[0] = sqrt((Xvtx-Bspt_x)*(Xvtx-Bspt_x) + (Yvtx-Bspt_y)*(Yvtx-Bspt_y));

		// Give a warning when Nmu is greater than our chosen Nmax from teh beginning
		if ((Nmu > Nmax) || (Emncand > Nmax)) { //  
			std::cout << "////////////////////////////////////////////////////////////////////////////////////////////////" << std::endl;
			std::cout << "//////////////////// WARNING! Nmu VARIABLE IS OVER THER LENHTH OF AN ARRAY! ////////////////////" << std::endl;
			std::cout << "////////////////////////////////////////////////////////////////////////////////////////////////" << std::endl;
		}

		nMuon = Nmu;
		// ********************** Moun loop **********************
		for (int iMuon = 0; iMuon < nMuon; ++iMuon) {
			// for (djklagj;ldskj) loop over tracks


			Muon_charge[iMuon] = Mucharge[iMuon];
			Muon_pt[iMuon] = Mupt[iMuon];
			Muon_eta[iMuon] = -log(tan(Muth[iMuon] / 2));
			Muon_phi[iMuon] = Muph[iMuon];
			Muon_mass[iMuon] = 0.105658;
			
			// Check Muqual
			// Set them first all to false (default)
			Muon_tightId[iMuon] = false;
			Muon_mediumId[iMuon] = false;
			Muon_softId[iMuon] = false;
			Muon_looseId[iMuon] = false;
			
			if (Muqual[iMuon] == Muqual_tightId) {
				Muon_tightId[iMuon] = true;
			}
			if (Muqual[iMuon] > Muqual_mediumId) {
				Muon_mediumId[iMuon] = true;
			}
			if (Muqual[iMuon] > Muqual_softId) {
				Muon_softId[iMuon] = true;
			}
			if (Muqual[iMuon] > Muqual_looseId) {
				Muon_looseId[iMuon] = true;
			}

			Muon_trkIdx[iMuon] = Mutrid[iMuon];
			Muon_ptErr[iMuon] = Muperr[iMuon] * sin(Muth[iMuon]);

			if (Muzufid[iMuon] > 0) {
				Muon_isPFcand[iMuon] = true;
			}
			else {
				Muon_isPFcand[iMuon] = false;
			}

			Muon_isTracker[iMuon] = true;
			
			if (Muqual[iMuon] > 0) {
				Muon_isGlobal[iMuon] = true;
			}
			else {
				Muon_isGlobal[iMuon] = false;
			}

			if (Muqual[iMuon] == 6) {
				Muon_isStandAlone[iMuon] = true;
			}
			else {
				Muon_isStandAlone[iMuon] = false;
			}

			Muon_pdgId[iMuon] = -13*Mucharge[iMuon];
			Muon_chi2[iMuon] = Muchid[iMuon];
			Muon_z[iMuon] = Muz[iMuon];
			Muon_dxy[iMuon] = Mudxy[iMuon];
			Muon_dz[iMuon] = Mudz[iMuon];
			Muon_pfreliso03_all[iMuon] = 1; // needs to be fixed!!!
			Muon_vtxIdx[iMuon] = Muvtxid[iMuon];
			Muon_vtxFlag[iMuon] = Muvtxfl[iMuon];
			Muon_jetIdx[iMuon] = Mujetid_a[iMuon];

			// calculate and fill Dimu structure
			std::cout << "   ***** Inside muon loop: " << iMuon << std::endl;
		}

		nElectron = Emncand;
		// // ********************** Electron loop **********************
		for (int iElectron = 0; iElectron < nElectron; iElectron++) {

			
			Float_t Electron_SCeta_x = abs(Emcalpos[iElectron][0] - PV_x);
			Float_t Electron_SCeta_y = abs(Emcalpos[iElectron][1] - PV_y);
			Float_t Electron_SCeta_z = abs(Emcalpos[iElectron][2] - PV_z);
			Electron_SCeta[iElectron] = sqrt(pow(Electron_SCeta_x, 2) + pow(Electron_SCeta_y, 2) + pow(Electron_SCeta_z, 2));
			Electron_charge[iElectron] = Emtrkq[iElectron]; // needed modification: put 0 if electron is not from track (e.g. from calorimetr)
			Electron_convDcot[iElectron] = -999.; // needed further investigation
			Electron_convDist[iElectron] = -999.; // needed further invessigation
			Electron_convVeto[iElectron] = false; // needed further investigation
			Electron_convVetoOld[iElectron] = false; // needed further investigation
			Electron_cutBased[iElectron] = 0; // needed further investigation
			Electron_deltaEtaSC[iElectron] = -999.; // needed further investigation
			Electron_deltaEtaSCtr[iElectron] = -999.; // needed further investigation
			Electron_deltaPhiSC[iElectron] = -999.; // needed further investigation
			Electron_deltaPhiSCtr[iElectron] = -999.; // needed further investigation
			Electron_dr03EcalRecHitSumEt[iElectron] = -999.; // needed further investigation
			Electron_dr03EcalRecHitSumEtOld[iElectron] = -999.; // needed further investigation
			Electron_dr03HcalDepth1TowerSumEt[iElectron] = -999.; // needed further investigation
			Electron_dr03HcalDepth1TowerSumEtOld[iElectron] = -999.; // needed further investigation
			Electron_dr03HcalTowerSumEt[iElectron] = -999.; // needed further investigation
			Electron_dr03TkSumPt[iElectron] = Emetneartrk[iElectron][1];
			Electron_dr03TkSumPtOld[iElectron] = Emetneartrk[iElectron][1];
			Electron_dxy[iElectron] = Emdcabeam[iElectron];
			Electron_dxyErr[iElectron] = 0.2;
			Electron_dzErr[iElectron] = 0.4;
			Electron_eInvMinusPInv[iElectron] = (1/Emcalene[iElectron]) - (1/Emtrkp[iElectron]);
			Electron_eInvMinusPInvOld[iElectron] = -999.; // needed further investigation
			Electron_eta[iElectron] = -log(tan(Emtrkth[iElectron] / 2));
			Electron_genPartIdx[iElectron] = -1; // needed further investigation
			Electron_hoe[iElectron] = -999.; // needed further investigation
			Electron_ip3d[iElectron] = Emdca[iElectron];
			Electron_isEB[iElectron] = false; // needed further investigation
			Electron_isEE[iElectron] = false; // needed further investigation
			Electron_isNano[iElectron] = true;
			Electron_isPFcand[iElectron] = false; // needed further investigation
			Electron_lostHits[iElectron] = '0'; // needed further investigation
			Electron_mass[iElectron] = 0.548579909065; // in MeV taken from Particle Data Group
			Electron_nNano[iElectron] = Emncand;
			Electron_pfRelIso03_all[iElectron] = Emenin[iElectron] / Empt[iElectron];
			Electron_pfRelIso03_chg[iElectron] = Emetneartrk[iElectron][1] / Empt[iElectron];
			Electron_ph[iElectron] = Emtrkph[iElectron];
			Electron_pt[iElectron] = Empt[iElectron];
			Electron_sieie[iElectron] = -999.; // needed furhter investigation
			Electron_sieieR1[iElectron] = -999.; // needed further investigation
			Electron_simId[iElectron] = -1;
			Electron_sip3d[iElectron] = -999.; // needed further investigation
			Electron_tightCharge[iElectron] = -1; // needed further investigation
			Electron_vtxIdx[iElectron] = -1;



			Int_t emTrackID = Emtrknr[iElectron];
			// ********************** Loop over all tracks *************************************************
			// ********************** Check track id and match it with electron track **********************
			for (int nTrk = 0; nTrk < Trk_ntracks; nTrk++) {
				if (Trk_id[nTrk] == emTrackID) {
					Electron_x[iElectron] = Trk_pca[nTrk][0];
					Electron_y[iElectron] = Trk_pca[nTrk][1];
					Electron_z[iElectron] = Trk_pca[nTrk][2];
					Electron_dz[iElectron] = abs(Electron_z[iElectron] - PV_z);
				}
			}

			std::cout << "   ***** Inside electron loop: " << iElectron << std::endl;
		}
		// ********************** Electron loop end ******************


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

	std::cout << "=========================== All particles are converted =====================================" << std::endl;


	// ==================================================================================================================================== 
	// ==================================================================================================================================== 
	// ==================================================================================================================================== 

}




int main() {

	zeustocms_t();

	return 0;
}










































///////////////////////////////////Plot some histograms//////////////////////////////////////////////

/*
	//Create histograms
  	TH1D *hMuon_charge = new TH1D("hMuon_charge", "Muon charge;Charge;Counts", 20, -2, 2);	hMuon_charge->SetDirectory(0);
  	TH1D *hMuon_pt = new TH1D("hMuon_pt", "Muon p_{T};Transverse momentum; Counts", 50, 0, 20);	hMuon_pt->SetDirectory(0);
  	TH1D *hMuon_eta = new TH1D("hMuon_eta", "Muon #eta;#eta;Counts", 100, -7, 7);	hMuon_eta->SetDirectory(0);
	TH1D *hMuon_phi = new TH1D("hMuon_phi", "Muon #phi;#phi;Counts", 100, 0, 7);	hMuon_phi->SetDirectory(0);

	//Loop over the entries of the inputChain
	for(int i=0; i<nEntriesInput; i++) {
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

//inputChain.Close();
//outputFile.Close();