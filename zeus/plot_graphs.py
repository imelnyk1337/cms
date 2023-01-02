import ROOT
import numpy as np

cms_variables = ["Electron_SCeta", "Electron_charge", "Electron_convDcot", "Electron_convDist", "Electron_convVeto", "Electron_convVetoOld",
                "Electron_cutBased", "Electron_deltaEtaSC", "Electron_deltaEtaSCtr", "Electron_deltaPhiSC", "Electron_deltaPhiSCtr", "Electron_dr03EcalRecHitSumEt", 
                "Electron_dr03EcalRecHitSumEtOld", "Electron_dr03HcalDepth1TowerSumEt", "Electron_dr03HcalDepth1TowerSumEtOld", "Electron_dr03HcalTowerSumEt", "Electron_dr03TkSumPt", 
                "Electron_dr03TkSumPtOld", "Electron_dxy", "Electron_dxyErr", "Electron_dz", "Electron_dzErr", "Electron_eInvMinusPInv", "Electron_eInvMinusPInvOld", "Electron_eta",
                "Electron_genPartIdx", "Electron_hoe", "Electron_ip3d", "Electron_isEB", "Electron_isEE", "Electron_isNano", "Electron_isPFcand", "Electron_mass",
                "Electron_nNano", "Electron_pfRelIso03_all", "Electron_pfRelIso03_chg", "Electron_phi", "Electron_pt", "Electron_sieie", "Electron_sieieR1", "Electron_simId",
                "Electron_sip3d", "Electron_tightCharge", "Electron_vtxIdx"] #, "Electron_x", "Electron_y", "Electron_z"

variables_info = {
    "Electron_SCeta": {
        "x1": -np.pi,
        "x2": np.pi,
    },
    "Electron_charge": {
        "x1": 0.,
        "x2": 2.,
    },
    "Electron_convDcot": {
        "x1": -4.,
        "x2": 4.,
    },
    "Electron_convDist": {
        "x1": -20.,
        "x2": 10.,
    },
    "Electron_convVeto": {
        "x1": 1.,
        "x2": 1.5,
    },
    "Electron_convVetoOld": {
        "x1": 0.,
        "x2": 1.,
    },
    "Electron_cutBased": {
        "x1": 0.,
        "x2": 4.,
    },
    "Electron_deltaEtaSC": {
        "x1": -0.15,
        "x2": 0.15,
    },
    "Electron_deltaEtaSCtr": {
        "x1": -0.03,
        "x2": 0.03,
    },  
    "Electron_deltaPhiSC": {
        "x1": -0.10,
        "x2": 0.10,
    },
    "Electron_deltaPhiSCtr": {
        "x1": -0.2,
        "x2": 0.2,
    },
    "Electron_dr03EcalRecHitSumEt": {
        "x1": 0.,
        "x2": 10.,
    },
    "Electron_dr03EcalRecHitSumEtOld": {
        "x1": 0.,
        "x2": 10.,
    },
    "Electron_dr03HcalDepth1TowerSumEt": {
        "x1": 0.,
        "x2": 6.,
    },
    "Electron_dr03HcalDepth1TowerSumEtOld": {
        "x1": 0.,
        "x2": 6.,
    },
    "Electron_dr03HcalTowerSumEt": {
        "x1": 0.,
        "x2": 7.,
    },
    "Electron_dr03TkSumPt": {
        "x1": 0.,
        "x2": 5.,
    },
    "Electron_dr03TkSumPtOld": {
        "x1": 0.,
        "x2": 15.,
    },
    "Electron_dxy": {
        "x1": -0.3,
        "x2": 1.1,
    },
    "Electron_dxyErr": {
        "x1": 0.,
        "x2": 0.25,
    },
    "Electron_dz": {
        "x1": 0.,
        "x2": 25.,
    },
    "Electron_dzErr": {
        "x1": 0.,
        "x2": 0.5,
    },
    "Electron_eInvMinusPInv": {
        "x1": -0.6,
        "x2": 0.2,
    },
    "Electron_eInvMinusPInvOld": {
        "x1": -4.,
        "x2": 1,
    },
    "Electron_eta": {
        "x1": -np.pi,
        "x2": np.pi,
    },
    "Electron_genPartIdx": {
        "x1": -1.5,
        "x2": 0.,
    },
    "Electron_hoe": {
        "x1": 0.,
        "x2": 0.5,
    },
    "Electron_ip3d": {
        "x1": 0.,
        "x2": 50.,
    },
    "Electron_isEB": {
        "x1": 0.,
        "x2": 1.5,
    },
    "Electron_isEE": {
        "x1": 0.,
        "x2": 1.5,
    },
    "Electron_isNano": {
        "x1": 1.,
        "x2": 1.5,
    },
    "Electron_isPFcand": {
        "x1": 1.,
        "x2": 1.05,
    },
    "Electron_mass": {
        "x1": -0.10,
        "x2": 0.10,
    },
    "Electron_nNano": {
        "x1": 0.,
        "x2": 5.,
    },
    "Electron_pfRelIso03_all": {
        "x1": 0.,
        "x2": 1.4,
    },
    "Electron_pfRelIso03_chg": {
        "x1": 0.,
        "x2": 1.5,
    },
    "Electron_phi": {
        "x1": -np.pi,
        "x2": np.pi,
    },
    "Electron_pt": {
        "x1": 0.,
        "x2": 150.,
    },
    "Electron_sieie": {
        "x1": 0.,
        "x2": 0.08,
    },
    "Electron_sieieR1": {
        "x1": 0.,
        "x2": 0.05,
    },
    "Electron_simId": {
        "x1": -1.,
        "x2": 0.,
    },
    "Electron_sip3d": {
        "x1": -20.,
        "x2": 20.,
    },
    "Electron_tightCharge": {
        "x1": 0.,
        "x2": 2.,
    },
    "Electron_vtxIdx": {
        "x1": -1.,
        "x2": 0.,
    }
    # "Electron_x": {
    #     "x1": 0.,
    #     "x2": ,
    # },
    # "Electron_y": {
    #     "x1": ,
    #     "x2": ,
    # },
    # "Electron_z": {
    #     "x1": ,
    #     "x2": ,
    # }
}


cms_file = ROOT.TFile("/home/ihor/cms/zeus/cmsdata_Run2011A_DoubleElectron.root")
zeus_file = ROOT.TFile("/home/ihor/cms/zeus/zeustocms_analysis.root")
# cms_chain = ROOT.TChain("Events")
# cms_chain.Add("/home/ihor/cms/zeus/cmsdata_Run2011A_DoubleElectron.root")

# zeus_chain = ROOT.TChain("Events")
# zeus_chain.Add("/home/ihor/cms/zeus/zeustocms_analysis.root")


def plot_comparison(variables: list, info: dict, cms_file, zeus_file):
    for variable in variables:
        canvas = ROOT.TCanvas("canvas", variable, 1200, 800)
        # canvas.Divide(1, 2)
        # histo_cms = ROOT.TH1D("histo_cms", "CMS Data", 64, info[variable]["x1"], info[variable]["x2"])
        # histo_zeus = ROOT.TH1D("histo_zeus", "ZEUS to CMS Data", 64, info[variable]["x1"], info[variable]["x2"])
        # canvas.cd(1)
        histo_cms = cms_file.Get(f"{variable}")
        xtitle = str(histo_cms.GetXaxis().GetTitle())
        histo_cms.Draw("e1")
        # canvas.cd(2)
        # histo_zeus = zeus_file.Get(f"{variable}")
        # histo_zeus.GetXaxis().SetTitle(xtitle)
        # if variable == "Electron_dr03EcalRecHitSumEt":
        #     histo_zeus.GetXaxis().SetRangeUser(-3.5, 3)
        #     ROOT.gPad.Modified()
        # histo_zeus.Draw("e1")
        canvas.Draw()
        canvas.Print(f"cmsplots/{variable}_cms_only.pdf")
        del histo_cms
        # del histo_zeus



def plot_double_e():
    ROOT.gStyle.SetOptStat("nemrous")
    canvas = ROOT.TCanvas("canvas", "canvas", 1200, 800)
    histo_cms = cms_file.Get("massZto2e")
    histo_cms.Draw("e1")
    canvas.Draw()
    canvas.Print(f"comparison/double_e.pdf")
    


def plot_e_energy():
    ROOT.gStyle.SetOptStat("nemrous")
    canvas = ROOT.TCanvas("canvas", "canvas", 1200, 800)
    histo_zeus = zeus_file.Get("Electron_energy")
    histo_zeus.Draw("e1")
    canvas.Draw()
    canvas.Print(f"comparison/Electron_energy.pdf")

def plot_q2():
    ROOT.gStyle.SetOptStat("nemrous")
    # ROOT.gPad.SetLogx()
    canvas = ROOT.TCanvas("canvas", "canvas", 1200, 800)
    histo_zeus = zeus_file.Get("Emq2da_ZEUS")
    histo_zeus.GetYaxis().SetTitle("log(Number of events)")
    histo_zeus.Draw("e1")
    canvas.SetLogy()
    canvas.Draw()
    canvas.Print(f"comparison/Emq2da_ZEUS.pdf")


def plot_y():
    ROOT.gStyle.SetOptStat("nemrous")
    # ROOT.gPad.SetLogx()
    canvas = ROOT.TCanvas("canvas", "canvas", 1200, 800)
    histo_zeus = zeus_file.Get("Emyda_ZEUS")
    histo_zeus.Draw("e1")
    canvas.Draw()
    canvas.Print(f"comparison/Emyda_ZEUS.pdf")

def plot_x():
    ROOT.gStyle.SetOptStat("nemrous")
    # ROOT.gPad.SetLogx()
    canvas = ROOT.TCanvas("canvas", "canvas", 1200, 800)
    histo_zeus = zeus_file.Get("Emxda_ZEUS")
    histo_zeus.GetYaxis().SetTitle("log(Number of events)")
    histo_zeus.Draw("e1")
    canvas.SetLogy()
    canvas.Draw()
    canvas.Print(f"comparison/Emxda_ZEUS.pdf")

def plot_theta():
    ROOT.gStyle.SetOptStat("nemrous")
    # ROOT.gPad.SetLogx()
    canvas = ROOT.TCanvas("canvas", "canvas", 1200, 800)
    histo_zeus = zeus_file.Get("Electron_theta")
    histo_zeus.Draw("e1")
    canvas.Draw()
    canvas.Print(f"comparison/Electron_theta.pdf")

def plot_xyz():
    ROOT.gStyle.SetOptStat("nemrous")
    canvas_x = ROOT.TCanvas("canvas_x", "canvas", 1200, 800)
    canvas_x.Divide(1, 2)
    canvas_x.cd(2)
    histo_zeus_x = zeus_file.Get("Electron_x")
    histo_zeus_x.Draw("e1")
    canvas_x.cd(1)
    histo_cms_x = cms_file.Get("Electron_x")
    histo_cms_x.Draw("e1")
    canvas_x.Draw()
    canvas_x.Print(f"comparison/Electron_x.pdf")

    canvas_y = ROOT.TCanvas("canvas_y", "canvas", 1200, 800)
    canvas_y.Divide(1, 2)
    canvas_y.cd(2)
    histo_zeus_y = zeus_file.Get("Electron_y")
    histo_zeus_y.Draw("e1")
    canvas_y.cd(1)
    histo_cms_y = cms_file.Get("Electron_y")
    histo_cms_y.Draw("e1")
    canvas_y.Draw()
    canvas_y.Print(f"comparison/Electron_y.pdf")

    canvas_z = ROOT.TCanvas("canvas_z", "canvas", 1200, 800)
    canvas_z.Divide(1, 2)
    canvas_z.cd(2)
    histo_zeus_z = zeus_file.Get("Electron_z")
    histo_zeus_z.Draw("e1")
    canvas_z.cd(1)
    histo_cms_z = cms_file.Get("Electron_z")
    histo_cms_z.Draw("e1")
    canvas_z.Draw()
    canvas_z.Print(f"comparison/Electron_z.pdf")





if __name__== "__main__":
    # ROOT.gPad.SetLogy()
    # print(ROOT.gPad.SetLogy(True))
    # plot_double_e()
    plot_comparison(cms_variables, variables_info, cms_file, zeus_file)
    # plot_e_energy()
    # plot_q2()
    # plot_x()
    # plot_theta()
    # plot_y()
    # plot_xyz()
    ROOT.gROOT.ProcessLine(".q")