#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"
#include "Math/Vector4Dfwd.h"
#include "Math/Vector4D.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TLatex.h"
#include "TStyle.h"
#include "string.h"
#include "math.h"   
#include "TFile.h"

using namespace ROOT::VecOps;

// Compute the invariant mass of two muon four-vectors
float computeInvariantMass(RVec<float>& pt, RVec<float>& eta, RVec<float>& phi, RVec<float>& mass) {
    ROOT::Math::PtEtaPhiMVector m1(pt[0], eta[0], phi[0], mass[0]);
    ROOT::Math::PtEtaPhiMVector m2(pt[1], eta[1], phi[1], mass[1]);
    return (m1 + m2).mass();
}

// auto getMinProb = [](RVec<float>& var) {
//     return ;
// };


struct Cuts
{
    std::string Emq2da_cut;
    std::string Emyda_cut;
    std::string PV_z_cut;
    std::string Electron_energy_cut;
};


auto applyCuts(Cuts& cuts) {
    
}


void analysis() {

    
    // use electron

    Cuts appliedCuts = {"Emq2da_ZEUS_0 > 125. && Emq2da_ZEUS_0 < 20000.",
                        "Emyda_ZEUS_0 > 0.2 && Emyda_ZEUS_0 < 0.6",
                        "PV_z > -30. && PV_z < 30.",
                        "Electron_energy_0 > 10."};


    gStyle->SetOptStat("nemrous");
    ROOT::EnableImplicitMT(8); 
    ROOT::RDataFrame df("Events", "/home/ihor/cms/zeus/zeustocms.root");
    auto df2 = df
    .Filter("nElectron > 0", "filter1") 
    .Filter("Electron_cutBased[0] > 3", "filter2")
    .Define("Electron_SCeta_0", "Electron_SCeta[0]")
    .Define("Electron_charge_0", "Electron_charge[0]")
    .Define("Electron_convDcot_0", "Electron_convDcot[0]")
    .Define("Electron_convDist_0", "Electron_convDist[0]")
    .Define("Electron_convVeto_0", "Electron_convVeto[0]")
    .Define("Electron_convVetoOld_0", "Electron_convVetoOld[0]")
    .Define("Electron_cutBased_0", "Electron_cutBased[0]")
    .Define("Electron_deltaEtaSC_0", "Electron_deltaEtaSC[0]")
    .Define("Electron_deltaEtaSCtr_0", "Electron_deltaEtaSCtr[0]")
    .Define("Electron_deltaPhiSC_0", "Electron_deltaPhiSC[0]")
    .Define("Electron_deltaPhiSCtr_0", "Electron_deltaPhiSCtr[0]")
    .Define("Electron_dr03EcalRecHitSumEt_0", "Electron_dr03EcalRecHitSumEt[0]")
    .Define("Electron_dr03EcalRecHitSumEtOld_0", "Electron_dr03EcalRecHitSumEtOld[0]")
    .Define("Electron_dr03HcalDepth1TowerSumEt_0", "Electron_dr03HcalDepth1TowerSumEt[0]")
    .Define("Electron_dr03HcalDepth1TowerSumEtOld_0", "Electron_dr03HcalDepth1TowerSumEtOld[0]")
    .Define("Electron_dr03HcalTowerSumEt_0", "Electron_dr03HcalTowerSumEt[0]")
    .Define("Electron_dr03TkSumPt_0", "Electron_dr03TkSumPt[0]")
    .Define("Electron_dr03TkSumPtOld_0", "Electron_dr03TkSumPtOld[0]")
    .Define("Electron_dxy_0", "Electron_dxy[0]")
    .Define("Electron_dxyErr_0", "Electron_dxyErr[0]")
    .Define("Electron_dz_0", "Electron_dz[0]")
    .Define("Electron_dzErr_0", "Electron_dzErr[0]")
    .Define("Electron_eInvMinusPInv_0", "Electron_eInvMinusPInv[0]")
    .Define("Electron_eInvMinusPInvOld_0", "Electron_eInvMinusPInvOld[0]")
    .Define("Electron_eta_0", "Electron_eta[0]")
    .Define("Electron_genPartIdx_0", "Electron_genPartIdx[0]")
    .Define("Electron_hoe_0", "Electron_hoe[0]")
    .Define("Electron_ip3d_0", "Electron_ip3d[0]")
    .Define("Electron_isEB_0", "Electron_isEB[0]")
    .Define("Electron_isEE_0", "Electron_isEE[0]")
    .Define("Electron_isNano_0", "Electron_isNano[0]")
    .Define("Electron_isPFcand_0", "Electron_isPFcand[0]")
    // .Define("Electron_lostHits_0", "Electron_lostHits[0]")
    .Define("Electron_mass_0", "Electron_mass[0]")
    .Define("Electron_nNano_0", "Electron_nNano[0]")
    .Define("Electron_pfRelIso03_all_0", "Electron_pfRelIso03_all[0]")
    .Define("Electron_pfRelIso03_chg_0", "Electron_pfRelIso03_chg[0]")
    .Define("Electron_phi_0", "Electron_phi[0]")
    .Define("Electron_pt_0", "Electron_pt[0]")
    .Define("Electron_sieie_0", "Electron_sieie[0]")
    .Define("Electron_sieieR1_0", "Electron_sieieR1[0]")
    .Define("Electron_simId_0", "Electron_simId[0]")
    .Define("Electron_sip3d_0", "Electron_sip3d[0]")
    .Define("Electron_tightCharge_0", "Electron_tightCharge[0]")
    .Define("Electron_vtxIdx_0", "Electron_vtxIdx[0]")
    .Define("Electron_x_0", "Electron_x[0]")
    .Define("Electron_y_0", "Electron_y[0]")
    .Define("Electron_z_0", "Electron_z[0]")

    .Define("Electron_energy", "Electron_pt * cosh(Electron_eta)")
    .Define("Electron_energy_0", "Electron_energy[0]")
    .Define("Electron_theta", "2 * atan(exp(-Electron_eta))")
    .Define("Electron_theta_0", "Electron_theta[0]")
    .Define("Emxda_ZEUS_0", "Emxda_ZEUS[0]")
    .Define("Emyda_ZEUS_0", "Emyda_ZEUS[0]")
    .Define("Emq2da_ZEUS_0", "Emq2da_ZEUS[0]")
    .Define("log_Emq2da_ZEUS_0", "log(Emq2da_ZEUS_0)")
    .Filter(appliedCuts.Emq2da_cut, "filter Q2")
    .Filter(appliedCuts.Emyda_cut, "filter y")
    .Filter(appliedCuts.PV_z_cut, "filter PV_z")
    .Filter(appliedCuts.Electron_energy_cut, "filter Electron_energy");



    auto h_Electron_SCeta_ptr = df2
    .Histo1D({"Electron_SCeta", "Electron_SCeta", 100u, -2.0, 2.0}, "Electron_SCeta_0");

    auto h_Electron_charge_ptr = df2
    .Histo1D({"Electron_charge", "Electron_charge", 100u, -2, 2}, "Electron_charge_0");

    auto h_Electron_convDcot_ptr = df2
    .Histo1D({"Electron_convDcot", "Electron_convDcot", 100u, -2, 2}, "Electron_convDcot_0");

    auto h_Electron_convDist_ptr = df2
    .Histo1D({"Electron_convDist", "Electron_convDist", 100u, -6., 2.}, "Electron_convDist_0");

    auto h_Electron_convVeto_ptr = df2
    .Histo1D({"Electron_convVeto", "Electron_convVeto", 100u, 1., 1.0006}, "Electron_convVeto_0");

    auto h_Electron_convVetoOld_ptr = df2
    .Histo1D({"Electron_convVetoOld", "Electron_convVetoOld", 100u, 1., 1.0006}, "Electron_convVetoOld_0");

    auto h_Electron_cutBased_ptr = df2
    .Histo1D({"Electron_cutBased", "Electron_cutBased", 100u, 0., 5.}, "Electron_cutBased_0");

    auto h_Electron_deltaEtaSC_ptr = df2
    .Histo1D({"Electron_deltaEtaSC", "Electron_deltaEtaSC", 100u, -0.10, 0.10}, "Electron_deltaEtaSC_0");

    auto h_Electron_deltaEtaSCtr_ptr = df2
    .Histo1D({"Electron_deltaEtaSCtr", "Electron_deltaEtaSCtr", 100u, -0.01, 0.01}, "Electron_deltaEtaSCtr_0");

    auto h_Electron_deltaPhiSC_ptr = df2
    .Histo1D({"Electron_deltaPhiSC", "Electron_deltaPhiSC", 100u, -0.06, 0.06}, "Electron_deltaPhiSC_0");

    auto h_Electron_deltaPhiSCtr_ptr = df2
    .Histo1D({"Electron_deltaPhiSCtr", "Electron_deltaPhiSCtr", 100u, -0.05, 0.06}, "Electron_deltaPhiSCtr_0");

    auto h_Electron_dr03EcalRecHitSumEt_ptr = df2
    .Histo1D({"Electron_dr03EcalRecHitSumEt", "Electron_dr03EcalRecHitSumEt", 100u, 0., 3.}, "Electron_dr03EcalRecHitSumEt_0");

    auto h_Electron_dr03EcalRecHitSumEtOld_ptr = df2
    .Histo1D({"Electron_dr03EcalRecHitSumEtOld", "Electron_dr03EcalRecHitSumEtOld", 100u, 0., 5.}, "Electron_dr03EcalRecHitSumEtOld_0");

    auto h_Electron_dr03HcalDepth1TowerSumEt_ptr = df2
    .Histo1D({"Electron_dr03HcalDepth1TowerSumEt", "Electron_dr03HcalDepth1TowerSumEt", 100u, 0., 1.}, "Electron_dr03HcalDepth1TowerSumEt_0");

    auto h_Electron_dr03HcalDepth1TowerSumEtOld_ptr = df2
    .Histo1D({"Electron_dr03HcalDepth1TowerSumEtOld", "Electron_dr03HcalDepth1TowerSumEtOld", 100u, 0., 1.}, "Electron_dr03HcalDepth1TowerSumEtOld_0");

    auto h_Electron_dr03HcalTowerSumEt_ptr = df2
    .Histo1D({"Electron_dr03HcalTowerSumEt", "Electron_dr03HcalTowerSumEt", 100u, 0., 1.5}, "Electron_dr03HcalTowerSumEt_0");

    auto h_Electron_dr03TkSumPt_ptr = df2
    .Histo1D({"Electron_dr03TkSumPt", "Electron_dr03TkSumPt", 100u, -0.1, 1.5}, "Electron_dr03TkSumPt_0");

    auto h_Electron_dr03TkSumPtOld_ptr = df2
    .Histo1D({"Electron_dr03TkSumPtOld", "Electron_dr03TkSumPtOld", 100u, -0.1, 10.}, "Electron_dr03TkSumPtOld_0"); 

    auto h_Electron_dxy_ptr = df2
    .Histo1D({"Electron_dxy", "Electron_dxy", 100u, -0.03, 2.}, "Electron_dxy_0"); 

    auto h_Electron_dxyErr_ptr = df2
    .Histo1D({"Electron_dxyErr", "Electron_dxyErr", 100u, 0., 0.01}, "Electron_dxyErr_0"); 

    auto h_Electron_dz_ptr = df2
    .Histo1D({"Electron_dz", "Electron_dz", 100u, -0.01, 0.2}, "Electron_dz_0");

    auto h_Electron_dzErr_ptr = df2
    .Histo1D({"Electron_dzErr", "Electron_dzErr", 100u, 0., 5.}, "Electron_dzErr_0"); 

    auto h_Electron_eInvMinusPInv_ptr = df2
    .Histo1D({"Electron_eInvMinusPInv", "Electron_eInvMinusPInv", 100u,  -0.02, 0.02}, "Electron_eInvMinusPInv_0"); 

    auto h_Electron_eInvMinusPInvOld_ptr = df2
    .Histo1D({"Electron_eInvMinusPInvOld", "Electron_eInvMinusPInvOld", 100u, -0.025, 0.025}, "Electron_eInvMinusPInvOld_0"); 

    auto h_Electron_eta_ptr = df2
    .Histo1D({"Electron_eta", "Electron_eta", 100u, -2.5, 3.14}, "Electron_eta_0");

    auto h_Electron_genPartIdx_ptr = df2
    .Histo1D({"Electron_genPartIdx", "Electron_genPartIdx", 100u, -2., 2.}, "Electron_genPartIdx_0");

    auto h_Electron_hoe_ptr = df2
    .Histo1D({"Electron_hoe", "Electron_hoe", 100u, 0., 0.04}, "Electron_hoe_0"); 

    auto h_Electron_ip3d_ptr = df2
    .Histo1D({"Electron_ip3d", "Electron_ip3d", 100u, 0., 0.25}, "Electron_ip3d_0"); 

    auto h_Electron_isEB_ptr = df2
    .Histo1D({"Electron_isEB", "Electron_isEB", 100u, -0.05, 1.04}, "Electron_isEB_0"); 

    auto h_Electron_isEE_ptr = df2
    .Histo1D({"Electron_isEE", "Electron_isEE", 100u, -0.05, 1.04}, "Electron_isEE_0"); 

    auto h_Electron_isNano_ptr = df2
    .Histo1D({"Electron_isNano", "Electron_isNano", 100u, 0., 1.04}, "Electron_isNano_0"); 

    auto h_Electron_isPFcand_ptr = df2
    .Histo1D({"Electron_isPFcand", "Electron_isPFcand", 100u, 0., 1.04}, "Electron_isPFcand_0");

    // auto h_Electron_lostHits_ptr = df2
    // .Histo1D({"h_Electron_lostHits", "Electron_lostHits", 100u, 0., 10.}, "Electron_lostHits_0"); 

    auto h_Electron_mass_ptr = df2
    .Histo1D({"Electron_mass", "Electron_mass", 100u, -0.08, 0.08}, "Electron_mass_0"); 

    auto h_Electron_nNano_ptr = df2
    .Histo1D({"Electron_nNano", "Electron_nNano", 100u, 0., 7.}, "Electron_nNano_0"); 

    auto h_Electron_pfRelIso03_all_ptr = df2
    .Histo1D({"Electron_pfRelIso03_all", "Electron_pfRelIso03_all", 100u, 0., 0.2}, "Electron_pfRelIso03_all_0"); 

    auto h_Electron_pfRelIso03_chg_ptr = df2
    .Histo1D({"Electron_pfRelIso03_chg", "Electron_pfRelIso03_chg", 100u, 0., 0.1}, "Electron_pfRelIso03_chg_0"); 

    auto h_Electron_phi_ptr = df2
    .Histo1D({"Electron_phi", "Electron_phi", 100u, -3.14, 3.14}, "Electron_phi_0"); 

    auto h_Electron_pt_ptr = df2
    .Histo1D({"Electron_pt", "Electron_pt", 100u, 0., 70.}, "Electron_pt_0"); 

    auto h_Electron_sieie_ptr = df2
    .Histo1D({"Electron_sieie", "Electron_sieie", 100u, 0.005, 0.03}, "Electron_sieie_0"); 

    auto h_Electron_sieieR1_ptr = df2
    .Histo1D({"Electron_sieieR1", "Electron_sieieR1", 100u, 0.005, 0.03}, "Electron_sieieR1_0"); 

    auto h_Electron_simId_ptr = df2
    .Histo1D({"Electron_simId", "Electron_simId", 100u, -1.1, -0.9}, "Electron_simId_0"); 

    auto h_Electron_sip3d_ptr = df2
    .Histo1D({"Electron_sip3d", "Electron_sip3d", 100u, 0., 3.5}, "Electron_sip3d_0"); 

    auto h_Electron_tightCharge_ptr = df2
    .Histo1D({"Electron_tightCharge", "Electron_tightCharge", 100u, 0., 2.}, "Electron_tightCharge_0"); 

    auto h_Electron_vtxIdx_ptr = df2
    .Histo1D({"Electron_vtxIdx", "Electron_vtxIdx", 100u, -2.0, 2.5}, "Electron_vtxIdx_0"); 

    auto h_Electron_x_ptr = df2
    .Histo1D({"Electron_x", "Electron_x", 100u, 0.055, 0.09}, "Electron_x_0"); 

    auto h_Electron_y_ptr = df2
    .Histo1D({"Electron_y", "Electron_x", 100u, 0.025, 0.055}, "Electron_y_0"); 

    auto h_Electron_z_ptr = df2
    .Histo1D({"Electron_z", "Electron_z", 100u, -15., 15.}, "Electron_z_0");


    auto h_PV_z_ptr = df2
    .Histo1D({"PV_z", "PV_z", 100u, -50., 50.}, "PV_z");

    auto h_Electron_energy_ptr = df2
    .Histo1D({"Electron_energy", "Electron_energy", 100u, 0., 120.}, "Electron_energy_0");

    auto h_Electron_theta_ptr = df2
    .Histo1D({"Electron_theta", "Electron_theta", 100u, 0.0, 3.0}, "Electron_theta_0");

    auto h_Emxda_ZEUS_ptr = df2
    .Histo1D({"Emxda_ZEUS", "Emxda_ZEUS", 100u, 0., 0.6}, "Emxda_ZEUS_0");

    auto h_Emyda_ZEUS_ptr = df2
    .Histo1D({"Emyda_ZEUS", "Emyda_ZEUS", 100u, 0.2, 0.6}, "Emyda_ZEUS_0");

    auto h_Emq2da_ZEUS_ptr = df2
    .Histo1D({"Emq2da_ZEUS", "Emq2da_ZEUS", 100u, 0., 20000}, "Emq2da_ZEUS_0");

    auto h_log_Emq2da_ZEUS_ptr = df2
    .Histo1D({"log_Emq2da_ZEUS", "log_Emq2da_ZEUS", 100u, 4.828, 10.}, "log_Emq2da_ZEUS_0");




    TH1D h_Electron_SCeta = *h_Electron_SCeta_ptr;
    TH1D h_Electron_charge = *h_Electron_charge_ptr;
    TH1D h_Electron_convDcot = *h_Electron_convDcot_ptr;
    TH1D h_Electron_convDist = *h_Electron_convDist_ptr;
    TH1D h_Electron_convVeto = *h_Electron_convVeto_ptr;
    TH1D h_Electron_convVetoOld = *h_Electron_convVetoOld_ptr;
    TH1D h_Electron_cutBased = *h_Electron_cutBased_ptr;
    TH1D h_Electron_deltaEtaSC = *h_Electron_deltaEtaSC_ptr;
    TH1D h_Electron_deltaEtaSCtr = *h_Electron_deltaEtaSCtr_ptr;
    TH1D h_Electron_deltaPhiSC = *h_Electron_deltaPhiSC_ptr;
    TH1D h_Electron_deltaPhiSCtr = *h_Electron_deltaPhiSCtr_ptr;
    TH1D h_Electron_dr03EcalRecHitSumEt = *h_Electron_dr03EcalRecHitSumEt_ptr;
    TH1D h_Electron_dr03EcalRecHitSumEtOld = *h_Electron_dr03EcalRecHitSumEtOld_ptr;
    TH1D h_Electron_dr03HcalDepth1TowerSumEt = *h_Electron_dr03HcalDepth1TowerSumEt_ptr;
    TH1D h_Electron_dr03HcalDepth1TowerSumEtOld = *h_Electron_dr03HcalDepth1TowerSumEtOld_ptr;
    TH1D h_Electron_dr03HcalTowerSumEt = *h_Electron_dr03HcalTowerSumEt_ptr;
    TH1D h_Electron_dr03TkSumPt = *h_Electron_dr03TkSumPt_ptr;
    TH1D h_Electron_dr03TkSumPtOld = *h_Electron_dr03TkSumPtOld_ptr;
    TH1D h_Electron_dxy = *h_Electron_dxy_ptr;
    TH1D h_Electron_dxyErr = *h_Electron_dxyErr_ptr;
    TH1D h_Electron_dz = *h_Electron_dz_ptr;
    TH1D h_Electron_dzErr = *h_Electron_dzErr_ptr;
    TH1D h_Electron_eInvMinusPInv = *h_Electron_eInvMinusPInv_ptr;
    TH1D h_Electron_eInvMinusPInvOld = *h_Electron_eInvMinusPInvOld_ptr;
    TH1D h_Electron_eta = *h_Electron_eta_ptr;
    TH1D h_Electron_genPartIdx = *h_Electron_genPartIdx_ptr;
    TH1D h_Electron_hoe = *h_Electron_hoe_ptr;
    TH1D h_Electron_ip3d = *h_Electron_ip3d_ptr;
    TH1D h_Electron_isEB = *h_Electron_isEB_ptr;
    TH1D h_Electron_isEE = *h_Electron_isEE_ptr;
    TH1D h_Electron_isNano = *h_Electron_isNano_ptr;
    TH1D h_Electron_isPFcand = *h_Electron_isPFcand_ptr;
    // TH1D h_Electron_lostHits = *h_Electron_lostHits_ptr;
    TH1D h_Electron_mass = *h_Electron_mass_ptr;
    TH1D h_Electron_nNano = *h_Electron_nNano_ptr;
    TH1D h_Electron_pfRelIso03_all = *h_Electron_pfRelIso03_all_ptr;
    TH1D h_Electron_pfRelIso03_chg = *h_Electron_pfRelIso03_chg_ptr;
    TH1D h_Electron_phi = *h_Electron_phi_ptr;
    TH1D h_Electron_pt = *h_Electron_pt_ptr;
    TH1D h_Electron_sieie = *h_Electron_sieie_ptr;
    TH1D h_Electron_sieieR1 = *h_Electron_sieieR1_ptr;
    TH1D h_Electron_simId = *h_Electron_simId_ptr;
    TH1D h_Electron_sip3d = *h_Electron_sip3d_ptr;
    TH1D h_Electron_tightCharge = *h_Electron_tightCharge_ptr;
    TH1D h_Electron_vtxIdx = *h_Electron_vtxIdx_ptr;
    TH1D h_Electron_x = *h_Electron_x_ptr;
    TH1D h_Electron_y = *h_Electron_y_ptr;
    TH1D h_Electron_z = *h_Electron_z_ptr;
    TH1D h_PV_z = *h_PV_z_ptr;
    TH1D h_Electron_energy = *h_Electron_energy_ptr;
    TH1D h_Electron_theta = *h_Electron_theta_ptr;
    TH1D h_Emxda_ZEUS = *h_Emxda_ZEUS_ptr;
    TH1D h_Emyda_ZEUS = *h_Emyda_ZEUS_ptr;
    TH1D h_Emq2da_ZEUS = *h_Emq2da_ZEUS_ptr;
    TH1D h_log_Emq2da_ZEUS = *h_log_Emq2da_ZEUS_ptr;

    TFile outputFile("zeustocms_analysis.root", "RECREATE");


    h_Electron_SCeta.Write();
    h_Electron_charge.Write();
    h_Electron_convDcot.Write();
    h_Electron_convDist.Write();
    h_Electron_convVeto.Write();
    h_Electron_convVetoOld.Write();
    h_Electron_cutBased.Write();
    h_Electron_deltaEtaSC.Write();
    h_Electron_deltaEtaSCtr.Write();
    h_Electron_deltaPhiSC.Write();
    h_Electron_deltaPhiSCtr.Write();
    h_Electron_dr03EcalRecHitSumEt.Write();
    h_Electron_dr03EcalRecHitSumEtOld.Write();
    h_Electron_dr03HcalDepth1TowerSumEt.Write();
    h_Electron_dr03HcalDepth1TowerSumEtOld.Write();
    h_Electron_dr03HcalTowerSumEt.Write();
    h_Electron_dr03TkSumPt.Write();
    h_Electron_dr03TkSumPtOld.Write();
    h_Electron_dxy.Write();
    h_Electron_dxyErr.Write();
    h_Electron_dz.Write();
    h_Electron_dzErr.Write();
    h_Electron_eInvMinusPInv.Write();
    h_Electron_eInvMinusPInvOld.Write();
    h_Electron_eta.Write();
    h_Electron_genPartIdx.Write();
    h_Electron_hoe.Write();
    h_Electron_ip3d.Write();
    h_Electron_isEB.Write();
    h_Electron_isEE.Write();
    h_Electron_isNano.Write();
    h_Electron_isPFcand.Write();
    // h_Electron_lostHits.Write();
    h_Electron_mass.Write();
    h_Electron_nNano.Write();
    h_Electron_pfRelIso03_all.Write();
    h_Electron_pfRelIso03_chg.Write();
    h_Electron_phi.Write();
    h_Electron_pt.Write();
    h_Electron_sieie.Write();
    h_Electron_sieieR1.Write();
    h_Electron_simId.Write();
    h_Electron_sip3d.Write();
    h_Electron_tightCharge.Write();
    h_Electron_vtxIdx.Write();
    h_Electron_x.Write();
    h_Electron_y.Write();
    h_Electron_z.Write();
    h_PV_z.Write();
    h_Electron_energy.Write();
    h_Electron_theta.Write();
    h_Emxda_ZEUS.Write();
    h_Emyda_ZEUS.Write();
    h_Emq2da_ZEUS.Write();
    h_log_Emq2da_ZEUS.Write();

    outputFile.Close();
    gROOT->ProcessLine(".q");

}

int main() {

    analysis();
    return 0;
}