### This is a short description of the used scripts and files 

#### The directory tree looks like that:

```
[imelnyk@naf-cms12]~/public% pwd
/afs/desy.de/user/i/imelnyk/public

[imelnyk@naf-cms12]~/public% tree
.
├── analysis.cpp
├── cms_data.C
├── Dimuon2011_14_RDF.C
├── plot_graphs.py
├── zeustocms
└── zeustocms.C

0 directories, 6 files
```

There are three logical parts of the work: **zeus**, **cms**, and **common**.

##### 1. ZEUS
The **zeus** part includes the data processing of the ZEUS n-tuples in order to convert the ZEUS variables to the common for CMS _NanoAOD_ format. The agreed variable names for both ZEUS and CMS experiments are different. Some of the ZUES variables are not present in the CMS Documentation, and some of the CMS variables written in the documentation don't have the corresponding analogs in the ZEUS. Therefore, some variables were just renamed, and a lot of variables were calculated or invented (the full list of the variables is shown in the Appendix).

###### 1.1 ZEUS to CMS conversion
The main code used for the electrons and muon (was done on the Summer Program 2021) conversion is **```zeustocms.C```**. The script includes the algorithm of the variables conversion, can be processed on a zeus-naf machine, and doesn't require any special software stack to be sourced from CVMFS. In order to run this program there are two available options: 

1. ```root -l -q zeustocms.C ```

2. ```g++ -o zeustocms `root-config --cflags --libs` zeustocms.C ```

The second command also may be specified by using the ```-g``` flag which means debugging mode. (and also any other GNU flag can be added)
The input ZEUS data is stored in the root files and set by the following command in the conversion script:

```
in_chain.Add("/pnfs/desy.de/dphep/online/zeus/z/ntup/07p/v08b/data/root/data_07p*.root");
```

where ```in_chain``` is a ```ROOT::TChain``` class object. 
The output file with the converted variables is ```zeustocms.root```

###### 1.2 ZEUS analysis
The first program that is used for this part is the **```analysis.cpp```**. This C++ script used the modern **RDataFrame** approach, which has a lot of advantages compared to the old ```ROOT::TTree::Draw``` method of data processing. The core idea of the RDataFrame is a so-called _lazy evaluation_, which is the part of the functional programming paradigm . This concept is implemented by using the _smart pointer_ (the C++ smart pointer, it allows users to make sure that any chunk of memory is not allocated in vain (similar to the links counter and garbage collector in Python), via processing the command only in the time when it is triggered). In the ```analysis.cpp``` some new variables were calculated according to the formulas in the final report. The output of this script is a root file that contains the histograms ```ROOT:TH1D```. The **```analysis.cpp```** also includes some cuts on the variables.


##### 2. CMS
The cms part of the work has the analysis code referred to as **```cms_data.C```**. This script is responsible for the CMS-side analysis. In the file, the invariant mass of _Z -> ee_ decay (Z is produced due to the Drell-Yan process) is calculated. This script reconstructs the invariant mass of the Z boson by combining two electrons whose kinematic distributions were cut off preliminary. Therefore, some cuts were applied to the CMS variables here. The input root file for this script is called
```Run2011A_DoubleElectron_merged.root```, and can be found in the

```root://eospublic.cern.ch//eos/opendata/cms/upload/NanoAODRun1/01-Jul-22/```

directory on the _eospublic_.
The **```cms_data.C```** also produces ```ROOT:TH1D``` objects.

##### 3. COMMON
In the common part of the work it was necessary to compare the converted ZEUS variables with the initial CMS variables obtained on the previous step. The comparison can help to have a look at the quality of the conversion, and adjust something... 

The **```plot_graphs.py```** is the simple python script for histogram drawing.



##### Appendix


###### A.1 Useful links
- [Documentation for CMS Monte Carlo](https://twiki.cern.ch/twiki/pub/CMSPublic/WorkBookNanoAODRun1/doc_DYJetsToLL_M-50_7TeV.html#Electron)
- [ZEUS ROOT Ntuples v08a EM](https://zeusdp.desy.de/ZEUS_ONLY/analysis/comntp/variables/v08/root_variables.html#EM)
- [Introduction to ZEUS analysis](https://www-zeus.desy.de/ZEUS_ONLY/analysis/primer/)
- [Azimuthal correlations in photoproduction and deep inelastic ep scattering at HERA](https://arxiv.org/pdf/2106.12377.pdf)
- [Measurement of e+p Neutral Current Deep Inelastic Scattering with a Longitudinally Polarised Positron Beam and X-ray Radiation Damage for Silicon Sensors](https://www2.physnet.uni-hamburg.de/services/biblio/dissertation/dissfbPhysik/___Volltexte/Friederike___Januschek/Friederike___Januschek.pdf)
- [Jets at High Q2 at HERA and Test Beam Measurements with the EUDET Pixel Telescope](https://zeusdp.desy.de/physics/qcd/thesis/joerg_behr.pdf)
- [```ROOT::RDataFrame``` Class Reference](https://root.cern/doc/master/classROOT_1_1RDataFrame.html)
- [Dataframe tutorials](https://root.cern.ch/doc/master/group__tutorial__dataframe.html)

###### A.2 Full list of the variables
Here is the full list of the variables (some of the variables have not been converted):

| CMS variable | CMS description | Type | Conversion Formula | Status | Detector |
| --- | --- | --- | --- | --- | --- |
| Electron_SCeta | Super Cluster eta | Float_t | ``` Electron_SCeta_x = Emcalpos[iElectron][0] - PV_x;Electron_SCeta_yEmcalpos[iElectron][1] - PV_y;Electron_SCeta_z = Emcalpos[iElectron][2] - PV_z;R =sqrt(pow(Electron_SCeta_x, 2) + pow(Electron_SCeta_y, 2));Theta = atan2(R,Electron_SCeta_z);Electron_SCeta[iElectron] = -log(tan(Theta / 2));```| True | cal |
| Electron_charge | electric charge | Int_t | Emtrkq[EmNcand] (put 0 if e is not from track) | True | trk |
| Electron_convDcot | NanoAODplus extension | Float_t |  | Undef |  |
| Electron_convDist | NanoAODplus extension | Float_t |  | Undef |  |
| Electron_convVeto | pass conversion veto | Bool_t | 1??? | Undef |  |
| Electron_convVetoOld | pass conversion veto | Bool_t |  | Undef |  |
| Electron_cutBased | cut-based ID Fall17 V2 (0:fail, 1:veto, 2:loose, 3:medium, 4:tight) | Int_t | ```if (Emprob[iElectron] < 0.001) { Electron_cutBased[iElectron] = 0;}else { if (Emtrknr[iElectron] == 0) { Electron_cutBased[iElectron] = 2;} else {Electron_cutBased[iElectron] = 4; // cut: Electron_cutBased == 4 -- from track}}``` |  | Undef |  |
| Electron_deltaEtaSC | delta eta (SC, ele) with sign | Float_t | Electron_SCeta - Electron_eta | Undef |  |
| Electron_deltaEtaSCtr | NanoAODplus extension | Float_t |  | Undef |  |
| Electron_deltaPhiSC | delta phi (SC,ele) with sign | Float_t | Calc phi from cluster - Emtrkph | Undef | cal |
| Electron_deltaPhiSCtr | NanoAODplus extension | Float_t |  | Undef |  |
| Electron_dr03EcalRecHitSumEt | Non-PF Ecal isolation within a delta R cone of 0.3 with electron pt > 35 GeV | Float_t | Needed further investigations | Undef |  |
| Electron_dr03EcalRecHitSumEtOld | Non-PF Ecal isolation within a delta R cone of 0.3 | Float_t | Needed further investigations | Undef |  |
| Electron_dr03HcalDepth1TowerSumEt | Non-PF Hcal isolation within a delta R cone of 0.3 with electron pt > 35 GeV | Float_t | Needed further investigations | Undef |  |
| Electron_dr03HcalDepth1TowerSumEtOld | Non-PF Hcal isolation within a delta R cone of 0.3 | Float_t | Needed further investigations | Undef |  |
| Electron_dr03HcalTowerSumEt | NanoAODplus extension | Float_t |  | Undef |  |
| Electron_dr03TkSumPt | Non-PF track isolation within a delta R cone of 0.3 with electron pt > 35 GeV | Float_t | Emetneartrk[EmNcand][1]  (make log scale) | True |  |
| Electron_dr03TkSumPtOld | old: Non-PF track isolation within a delta R cone of 0.3 | Float_t | Emetneartrk[EmNcand][1]  (make log scale) | True |  |
| Electron_dxy | dxy (with sign) wrt first PV, in cm | Float_t | Emdcabeam[EmNcand] hist range up to 2 | True |  |
| Electron_dxyErr | dxy uncertainty, in cm | Float_t | dxy = 0.2 cm | Undef |  |
| Electron_dz | dz (with sign) wrt first PV, in cm | Float_t | sqrt( Emdca[EmNcand]^2 - Emdcabeam[EmNcand]^2 )  | True |  |
| Electron_dzErr | dz uncertainty, in cm | Float_t | dz = 0.4 cm | True |  |
| Electron_eInvMinusPInv | 1/E_SC - 1/p_trk | Float_t | 1/Emcalene[EmNcand] - 1/Emtrkp[EmNcand] | True |  |
| Electron_eInvMinusPInvOld | NanoAODplus Run1 | Float_t |  | Undef |  |
| Electron_eta | eta | Float_t | — ln ( tan ( Emtrkth[ EmNcand ] / 2 ) ) | True | trk |
| Electron_genPartIdx | GenPart index of the associated true electron (-1 if none) | Int_t (index to Genpart) | -1 | True |  |
| Electron_hoe | H over E | Float_t | (1/Emfemc[EmNcand]) - 1 | Undef |  |
| Electron_ip3d | 3D impact parameter wrt first PV, in cm | Float_t | Emdca[EmNcand] | True | trk |
| Electron_isEB | NanoAODplus extension | Bool_t | -120 < Empos[2] < 200 | True |  |
| Electron_isEE | NanoAODplus extension | Bool_t | !(-120 < Empos[2] < 200) | True |  |
| Electron_isNano | NanoAODplus extension | Bool_t | true for all | Undef |  |
| Electron_isPFcand | NanoAODplus extension | Bool_t | true for all | Undef |  |
| Electron_lostHits | NanoAODplus extension | UChar_t |  | Undef |  |
| Electron_mass | mass | Float_t | take it from PDG | True | cal/trk |
| Electron_nNano | electron that follows nanoAOD cut (pT < 5 GeV) | UInt_t | Emncand | True |  |
| Electron_pfRelIso03_all | PF relative isolation dR=0.3 | Float_t | Emenin[EmNcand] / Empt[EmNcand]  | True | cal |
| Electron_pfRelIso03_chg | PF relative isolation dR=0.3, charged component | Float_t | Emetneartrk[EmNcand][1] / Empt[EmNcand]  | True | cal |
| Electron_phi | phi | Float_t | Emtrkph[EmNcand] | True | trk |
| Electron_pt | pt | Float_t | Empt[EmNcand]  | True | cal |
| Electron_sieie | sigma_IetaIeta of the supercluster, calculated with full 5x5 region | Float_t |  | Undef |  |
| Electron_sieieR1 | NanoAODplus Run1 | Float_t |  | Undef |  |
| Electron_simId | GenPart Id of the associated true electron (-1 if none) | Int_t | -l for every cand | True |  |
| Electron_sip3d | 3D impact parameter significance wrt first PV, in cm | Float_t |  | Undef |  |
| Electron_tightCharge | tight charge criteria | Int_t | 2??? | Undef |  |
| Electron_vtxIdx | index of the associated primary vertex (-1 if none) | Int_t (index to Vtx) | -l for each cand | True |  |
| Electron_x | x point of closest approach to beam line, in cm | Float_t | Trk_pca[trk_ntracks][0] | True |  |
| Electron_y | y point of closest approach to beam line, in cm | Float_t | Trk_pca[trk_ntracks][1] | True |  |
| Electron_z | z point of closest approach to beam line, in cm | Float_t | Trk_pca[trk_ntracks][2] | True |  |
| nElectron | number of Electrons | UInt_t | Emncand | True |  |
