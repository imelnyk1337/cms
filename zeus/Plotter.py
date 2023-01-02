import ROOT as r
r.gROOT.SetBatch(1)

# HDA = Higgs Demo Analyzer
# NA = Nano Analyzer
def make_plots(file_old,file_new,ID,name_OLD,name_NEW):

    pathHDA = file_old #"DoubleMu2011/Data_DoubleMu2011_RunA_nanoAOD_OLD.root"
    pathNA = file_new #"DoubleMu2011/Data_DoubleMu2011_RunA_nanoAOD_NEW.root"
    setId = "_"+ID #"_DoubleMu11_RunA"
    rm_letter = "_RunA"
    rm_letter_2 = "_RunB+C"
    for letter in ID:
        ID = ID.replace(rm_letter,'')
        ID = ID.replace(rm_letter_2,'')
    print("The plots will go to:",ID)

# ----------------------------------------------------------------------------------------------------
    NA = r.TFile.Open(pathNA)

    dirNA = NA.GetDirectory("Events")

# ----------------------------------------------------------------------------------------------------
    HDA = r.TFile.Open(pathHDA)

    dirHDA = HDA.GetDirectory("Events")

# ----------------------------------------------------------------------------------------------------
    histoNames = [key.GetName() for key in dirNA.GetListOfKeys()]



# '''
#     # Nombres de los ejes x
    xTitles = ["Number of good reco muons", "Muon p_{T} before the cuts (GeV)", "Muon #eta before the cuts", "Muon p_{T} after the cuts (GeV)", "Muon #eta after the cuts", "m_{#mu#mu} (GeV)", "Relative isolation", "SIP_{3D}", "Transverse impact parameter w.r.t. BS"]
# '''

# Estilo general
    r.gStyle.SetOptTitle(0)
    r.gStyle.SetLineWidth(2)

# Bucle a todos los histogramas

######################################################################################33
    for nn in range(len(histoNames)):
        c = r.TCanvas('c', 'c', 10, 10, 1200, 900)
        c.Divide(1,2)

        low = 0.4
        factor = (1.-low)/low

        c.GetPad(1).SetPad('upperPad', '', 0, low, 1, 1, 0, 0, 0)
        c.GetPad(2).SetPad('lowerPad', '', 0, 0, 1, low, 0, 0, 0)

        c.cd()
        c.cd(1)
        
        # Pad
        r.gPad.SetTopMargin(0.02)
        r.gPad.SetBottomMargin(0.04)
        r.gPad.SetLogy()
        r.gPad.SetTicks(1)

        # Get histo
        hHDA = dirHDA.Get(histoNames[nn])
        hNA = dirNA.Get(histoNames[nn])

        # Settings
        hHDA.SetLineColor(r.kAzure+2)
        hHDA.SetLineWidth(2)

        hNA.SetLineColor(r.kRed)
        hNA.SetLineWidth(2)
        
        max_HDA = 0
        max_NA = 0
        for n in range(0, hHDA.GetXaxis().GetNbins()):
            if hHDA.GetBinContent(n) > max_HDA: max_HDA = hHDA.GetBinContent(n)

        for n in range(0, hNA.GetXaxis().GetNbins()):
            if hNA.GetBinContent(n) > max_NA: max_NA = hNA.GetBinContent(n)

    # ------------------------------------------------------------------
        if max_HDA > max_NA: 

            hHDA.GetXaxis().SetLabelSize(0)
            hHDA.GetYaxis().SetTitle("Events")
            hHDA.GetYaxis().SetLabelSize(0.03)
            hHDA.GetYaxis().SetTitleSize(0.04)
            hHDA.SetTitle("")
            hHDA.SetTitleSize(0)

            hHDA.Draw("axis")
            hHDA.Draw("hist, sames")
            hNA.Draw("hist, sames")   # sames instead of same draws both statistic boxes

            r.gPad.Update()           # necesario para modificar la caja de estadistica       
        
            stHDA = hHDA.FindObject("stats")
            stHDA.SetName(name_OLD)
            stHDA.SetTextColor(r.kAzure+2)

            stNA = hNA.FindObject("stats")
            stNA.SetName(name_NEW)
            stNA.SetTextColor(r.kRed)

            # Posicion de la caja HDA
            X1 = stHDA.GetX1NDC()
            Y1 = stHDA.GetY1NDC()
            X2 = stHDA.GetX2NDC()
            Y2 = stHDA.GetY2NDC()

            # Posicion de la caja NA
            stNA.SetX1NDC(X1)
            stNA.SetX2NDC(X2)
            stNA.SetY1NDC(Y1-(Y2-Y1))
            stNA.SetY2NDC(Y1)

            # Ejes
            hHDA.GetYaxis().SetLabelSize(0.06)
            hHDA.GetYaxis().SetTitleSize(0.08)
            hHDA.GetXaxis().SetNdivisions(6, 5, 0, r.kTRUE)
            hHDA.GetYaxis().SetNdivisions(5, 5, 0, r.kTRUE)
            hHDA.GetXaxis().SetLabelSize(0)
            hHDA.GetYaxis().SetTitleOffset(0.65)

        # ------------------------------------------------------------------
        elif max_HDA <= max_NA:

            hNA.GetXaxis().SetLabelSize(0)
            hNA.GetYaxis().SetTitle("Events")
            hNA.GetYaxis().SetLabelSize(0.03)
            hNA.GetYaxis().SetTitleSize(0.04)

            hNA.SetTitle("")
            hNA.SetTitleSize(0)

            hNA.Draw("axis")
            hHDA.Draw("hist,sames")
            hNA.Draw("hist,sames")

            r.gPad.Update()        
            
            stHDA = hHDA.FindObject("stats")
            stHDA.SetName(name_OLD)
            stHDA.SetTextColor(r.kAzure+2)

            stNA = hNA.FindObject("stats")
            stNA.SetName(name_NEW)
            stNA.SetTextColor(r.kRed)

            # Posicion de la caja HDA
            X1 = stHDA.GetX1NDC()
            Y1 = stHDA.GetY1NDC()
            X2 = stHDA.GetX2NDC()
            Y2 = stHDA.GetY2NDC()

            # Posicion de la caja NA
            stNA.SetX1NDC(X1)
            stNA.SetX2NDC(X2)
            stNA.SetY1NDC(Y1-(Y2-Y1))
            stNA.SetY2NDC(Y1)
            
            # Ejes
            hNA.GetYaxis().SetLabelSize(0.06)
            hNA.GetYaxis().SetTitleSize(0.08)
            hNA.GetXaxis().SetNdivisions(6, 5, 0, r.kTRUE)
            hNA.GetYaxis().SetNdivisions(5, 5, 0, r.kTRUE)
            hNA.GetXaxis().SetLabelSize(0)
            hNA.GetYaxis().SetTitleOffset(0.65)
        # ------------------------------------------------------------------

        l = r.TLegend(X1-0.02-0.18, Y1+0.05, X1-0.02, Y2) # x0, y0, x1, y1

        l.AddEntry(hHDA, name_OLD, "l")
        l.AddEntry(hNA, name_NEW, "l")
        l.SetTextSize(0.04)
        l.SetBorderSize(2)
        l.SetLineColor(r.kBlack)
        l.SetFillColor(r.kWhite)
        l.Draw()

        # -----------------------------------------------------------------
        # cd(2) escala diferente -> parametros x factor
        c.cd(2)

        # Pad
        r.gPad.SetTopMargin(0.04)
        r.gPad.SetBottomMargin(0.28)
        r.gPad.SetTicks(1)
        r.gPad.SetGridy(1)

        ratio = hHDA.Clone('Ratio')
        ratio.Divide(hNA)
        ratio.SetLineColor(r.kGray+3)

        #ratio.GetXaxis().SetTitle(xTitles[nn])
        ratio.GetXaxis().SetLabelSize(0.06*factor)
        ratio.GetXaxis().SetTitleSize(0.08*factor)
        ratio.GetXaxis().SetTitleOffset(1.6/factor)
        ratio.GetYaxis().SetTitle("Ratio")

        ratio.GetYaxis().SetLabelSize(0.06*factor)

        ratio.GetYaxis().SetTitleSize(0.08*factor)
        ratio.GetYaxis().SetTitleOffset(0.6/factor)
        ratio.GetYaxis().SetTitleOffset(0.6/factor)

        ratio.GetXaxis().SetNdivisions(6, 5, 0, r.kTRUE)
        ratio.GetYaxis().SetNdivisions(5, 5, 0, r.kTRUE)

        ratio.SetStats(0)

        ratio.Draw('hist')

        c.Update()
        c.Print(ID+"/plots/"+histoNames[nn] + setId + ".png", "png")
