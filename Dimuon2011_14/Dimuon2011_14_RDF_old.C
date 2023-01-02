#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"
#include "Math/Vector4Dfwd.h"
#include "Math/Vector4D.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TLatex.h"
#include "TStyle.h"

using namespace ROOT::VecOps;

void Dimuon2011_14_RDF() {
    // Enable multi-threading
    // The default here is set to a single thread. You can choose the number of threads based on your system.
    ROOT::EnableImplicitMT(12);

    // Run over double muon sample
    ROOT::RDataFrame df_DoubleMu("Events", "/pnfs/desy.de/cms/tier2/store/user/geiser/Data2011/nanoAODplus_v1/DoubleMu/CMS_Run2011A_DoubleMu_AOD_12Oct2013-v1.root");
    auto filter1 = df_DoubleMu.Filter("run<170000", "Run number")
        .Filter("Trig_DoubleMuThresh>12", "Dimuon threshold")
        .Define("Dimu_mass_cut", [](uint nDimu, const RVec<int>& Dimu_charge, const RVec<float>& Dimu_mass, const RVec<int>& Dimu_t1muIdx, const RVec<int>& Dimu_t2muIdx,
            const RVec<float>& Muon_pt, const RVec<bool>& Muon_mediumId) -> RVec<float> {
            RVec<float> Dimu_mass_cut;
            for(size_t i = 0; i < nDimu; i++)
                if(Dimu_charge[i]==0 &&
                    Muon_pt[Dimu_t1muIdx[i]] > 6 &&
                    Muon_pt[Dimu_t2muIdx[i]] > 6 &&
                    Muon_mediumId[Dimu_t1muIdx[i]] &&
                    Muon_mediumId[Dimu_t2muIdx[i]])
                    Dimu_mass_cut.push_back(Dimu_mass[i]);
            return Dimu_mass_cut;
        }, {"nDimu", "Dimu_charge", "Dimu_mass", "Dimu_t1muIdx", "Dimu_t2muIdx", "Muon_pt", "Muon_mediumId"});
    auto h_dimulog1_ptr = filter1.Define("hist1data", "log10(Dimu_mass_cut)")
        .Define("hist1weight", "2./log(10.)/Dimu_mass_cut")
        .Histo1D({"h_dimulog1", "h_dimulog1", 620,-0.4, 2.7}, "hist1data", "hist1weight");
    auto report1 = filter1.Report();

    // Run over MuOnia sample
    ROOT::RDataFrame df_MuOnia("Events", {
        "/nfs/dust/cms/user/yangq2/NanoAODplus_v1/output/data/2011/MuOnia/RunA/00000.root",
        "/nfs/dust/cms/user/yangq2/NanoAODplus_v1/output/data/2011/MuOnia/RunA/00001.root",
        "/nfs/dust/cms/user/yangq2/NanoAODplus_v1/output/data/2011/MuOnia/RunA/20000.root",
        "/nfs/dust/cms/user/yangq2/NanoAODplus_v1/output/data/2011/MuOnia/RunA/210000.root"});
    auto filter4 = df_MuOnia.Filter("run<170000", "Run number")
        .Filter("!(Alsoon_DoubleMu && Trig_DoubleMuThresh>12)", "Dimuon threshold & sample overlap")
        .Filter("Trig_JpsiThresh !=0", "J/Psi threshold")
        .Filter("!HLT_DoubleMu4_LowMass_Displaced && !HLT_DoubleMu4p5_LowMass_Displaced && !HLT_DoubleMu5_LowMass_Displaced && !HLT_Dimuon6p5_LowMass_Displaced && !HLT_Dimuon7_LowMass_Displaced && !HLT_DoubleMu4_Jpsi_Displaced && !HLT_DoubleMu5_Jpsi_Displaced && !HLT_Dimuon6p5_Jpsi_Displaced && !HLT_Dimuon7_Jpsi_Displaced && !HLT_Mu5_L2Mu2", "HLT")
        .Define("Dimu_mass_cut", [](uint nDimu, const RVec<int>& Dimu_charge, const RVec<float>& Dimu_mass, const RVec<int>& Dimu_t1muIdx, const RVec<int>& Dimu_t2muIdx,
            const RVec<float>& Muon_pt, const RVec<bool>& Muon_mediumId) -> RVec<float> {
            RVec<float> Dimu_mass_cut;
            for(size_t i = 0; i < nDimu; i++)
                if(Dimu_charge[i] == 0 &&
                    Muon_pt[Dimu_t1muIdx[i]] > 3 &&
                    Muon_pt[Dimu_t2muIdx[i]] > 3 &&
                    Muon_mediumId[Dimu_t1muIdx[i]] &&
                    Muon_mediumId[Dimu_t2muIdx[i]] &&
                    Dimu_mass[i] > 2)
                    Dimu_mass_cut.push_back(Dimu_mass[i]);
            return Dimu_mass_cut;
        }, {"nDimu", "Dimu_charge", "Dimu_mass", "Dimu_t1muIdx", "Dimu_t2muIdx", "Muon_pt", "Muon_mediumId"});
    auto h_dimulog4_ptr = filter4.Define("hist4data", "log10(Dimu_mass_cut)")
        .Define("hist4weight", "2./log(10.)/Dimu_mass_cut")
        .Histo1D({"h_dimulog4", "h_dimulog4", 620,-0.4, 2.7}, "hist4data", "hist4weight");
    auto report4 = filter4.Report();
    // Add remaining filters here:
    // auto filter2 = ...
    // auto h_dimulog2_ptr = ...
    // auto report2 = filter2.Report();

    // Event loops is run here (once each)
    TH1D& h_dimulog1 = *h_dimulog1_ptr;
    TH1D& h_dimulog4 = *h_dimulog4_ptr;
    //TH1D& h_dimulog2 = *h_dimulog2_ptr;


    TCanvas *c = new TCanvas("c1","c1",1);
    c->SetLogy();
    gStyle->SetOptStat(0); // remove box
    h_dimulog1.SetTitle("Dimuon mass spectrum 2011 7 TeV (1.2 fb-1)");
    h_dimulog1.GetXaxis()->SetTitle("Invariant Log10(Mass) for Nmuon>=2 (in log10(m/GeV/c^2))");
    h_dimulog1.GetYaxis()->SetTitle("Number of Events/10 MeV");
    h_dimulog1.SetMinimum(0.02);
    h_dimulog1.SetMaximum(3.E6);

    h_dimulog1.SetFillColor(8); // Green
    h_dimulog1.Draw("hist");
    h_dimulog4.SetFillColor(9);
    h_dimulog4.Draw("hist same");
    // h_dimulog2.SetFillColor(10);
    // h_dimulog2.Draw("hist same");

    c->Print("Dimuon2011_14_RDF.png");
    c->SaveAs("Dimuon2011_14_RDF.pdf");

    TFile Dimuon2011("Dimuon2011_14_RDF.root", "RECREATE"); 
    h_dimulog1.Write();
    h_dimulog4.Write();

    report1->Print();
    report4->Print();
    // report2->Print();
}


int main() {
    Dimuon2011_14_RDF();
}
