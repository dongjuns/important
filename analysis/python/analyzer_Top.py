import ROOT, os, getopt, sys, array, math, glob, json
from ROOT import * 
from array import array


### Rochester ###
ROOT.gROOT.LoadMacro("/cms/scratch/jdj0715/nanoAOD/src/nano/analysis/src/RoccoR.cc+")
roc = ROOT.std.string("/cms/scratch/jdj0715/nanoAOD/src/nano/analysis/data/rcdata.2016.v3/")
rocCor = ROOT.RoccoR(roc)



### Pileup Weight ###
ROOT.gROOT.LoadMacro("/cms/scratch/jdj0715/nanoAOD/src/nano/analysis/scripts/WeightCalculatorFromHistogram.cc+")
pufile_mc="/cms/scratch/jdj0715/nanoAOD/src/nano/analysis/data/pu_root/pileup_profile_Spring16.root"
fmc = ROOT.TFile(pufile_mc)
pufile_data="/cms/scratch/jdj0715/nanoAOD/src/nano/analysis/data/pu_root/PileupData_GoldenJSON_Full2016.root"
fmcrd = ROOT.TFile(pufile_data)


hist_mc = fmc.Get("pu_mc")
hist_mc.SetDirectory(0)

hist_data = fmcrd.Get("pileup")
hist_data.SetDirectory(0)

puWeight = ROOT.WeightCalculatorFromHistogram(hist_mc, hist_data, True, True, False)

### Make TTREE ### 
FileArg = sys.argv
tempdir = FileArg[1]
Dirname = "/cms/scratch/jdj0715/nanoAOD/src/nano/analysis/test/Results/Nano_NewCut/%s/"%tempdir
if not os.path.isdir(Dirname):
    os.makedirs(Dirname)

temp = FileArg[2].split('/').pop()
cattree = Dirname+temp

#print cattree
f = ROOT.TFile(cattree, "recreate")
ALL = ROOT.TTree("nEvent", "nEvent")
"""
Cat1 = ROOT.TTree("Cat1", "Cat1")
Cat2 = ROOT.TTree("Cat2", "Cat2")
Cat3 = ROOT.TTree("Cat3", "Cat3")
Cat4 = ROOT.TTree("Cat4", "Cat4")
Cat5 = ROOT.TTree("Cat5", "Cat5")
Cat6 = ROOT.TTree("Cat6", "Cat6")
Cat7 = ROOT.TTree("Cat7", "Cat7")
Cat8 = ROOT.TTree("Cat8", "Cat8")
Cat9 = ROOT.TTree("Cat9", "Cat9")
Cat10 = ROOT.TTree("Cat10", "Cat10")
"""

### Variables ###
Dilep = ROOT.TLorentzVector()
Mu1 = ROOT.TLorentzVector()
Mu2 = ROOT.TLorentzVector()
GenLep1 = ROOT.TLorentzVector()
GenLep2 = ROOT.TLorentzVector()

Mu_Pt = ROOT.std.vector('float')()
Mu_Eta = ROOT.std.vector('float')()
Mu_Charge = ROOT.std.vector('float')()
Mu_Phi = ROOT.std.vector('float')()
Mu_M = ROOT.std.vector('float')()

El_Pt = ROOT.std.vector('float')()
El_Eta = ROOT.std.vector('float')()
El_Charge = ROOT.std.vector('float')()
El_Phi = ROOT.std.vector('float')()
El_M = ROOT.std.vector('float')()

Jet_Pt = ROOT.std.vector('float')()
Jet_Eta = ROOT.std.vector('float')()
Jet_CSVV2 = ROOT.std.vector('float')()
Jet_M = ROOT.std.vector('float')()
Jet_Phi = ROOT.std.vector('float')()

Event_No = array("i",[0])
Event_Total = array("i",[0])
Nu_Mu = array("i",[0])
Nu_El = array("i",[0])
Nu_Jet = array("i",[0])
Nu_BJet = array("i",[0])
Nu_NonBJet = array("i",[0])
genweight = array("f",[0])
puweight = array("f",[0])
b_weight = array("f",[0])
GJson = array("f",[0])
Step = array("i",[0])
NewStep = array("i",[0])
#Step0 = array("i",[0])
#Step1 = array("i",[0])
#Step2 = array("i",[0])
#Step3 = array("i",[0])
#Step4 = array("i",[0])
#Step5 = array("i",[0])
#Step6 = array("i",[0])

Event_Tot = ROOT.TH1D("Event_total", "Event_total" ,1,0,1)
genweights = ROOT.TH1D("genweight", "genweight" , 1,0,1)
weight = ROOT.TH1D("weight", "weight", 1,0,1)
cutFlow = ROOT.TH1I("cutflow", "cutflow", 11, -0.5, 10.5)

### Branches ###
ALL.Branch("Event_No", Event_No, "Event_No/I")
ALL.Branch("Step", Step, "Step/I")
#ALL.Branch("Step1", Step1, "Step1/O")
#ALL.Branch("Step2", Step2, "Step2/O")
#ALL.Branch("Step3", Step3, "Step3/O")
#ALL.Branch("Step4", Step4, "Step4/O")
#ALL.Branch("Step5", Step5, "Step5/O")
#ALL.Branch("Step6", Step6, "Step6/O")
ALL.Branch("Dilep", "TLorentzVector", Dilep)
ALL.Branch("Mu1", "TLorentzVector", Mu1)
ALL.Branch("Mu2", "TLorentzVector", Mu2)
ALL.Branch("Nu_Mu", Nu_Mu, "Nu_Mu/I")
ALL.Branch("Mu_Pt", Mu_Pt)
ALL.Branch("Mu_Eta", Mu_Eta)
ALL.Branch("Nu_El", Nu_El, "Nu_El/I")
ALL.Branch("El_Pt", El_Pt)
ALL.Branch("El_Eta", El_Eta)
ALL.Branch("Nu_Jet", Nu_Jet, "Nu_Jet/I")
ALL.Branch("Jet_Pt", Jet_Pt)
ALL.Branch("Jet_Eta", Jet_Eta)
ALL.Branch("Nu_BJet", Nu_BJet, "Nu_BJet/I")
ALL.Branch("genweight", genweight, "genweight/F")
ALL.Branch("puweight", puweight, "puweight/F")
"""
Cat10.Branch("Event_No", Event_No, "Event_No/I")
Cat10.Branch("Dilep", "TLorentzVector", Dilep)
Cat10.Branch("Mu1", "TLorentzVector", Mu1)
Cat10.Branch("Mu2", "TLorentzVector", Mu2)
Cat10.Branch("Nu_Mu", Nu_Mu, "Nu_Mu/I")
Cat10.Branch("Mu_Pt", Mu_Pt)
Cat10.Branch("Mu_Eta", Mu_Eta)
Cat10.Branch("Nu_El", Nu_El, "Nu_El/I")
Cat10.Branch("El_Pt", El_Pt)
Cat10.Branch("El_Eta", El_Eta)
Cat10.Branch("Nu_Jet", Nu_Jet, "Nu_Jet/I")
Cat10.Branch("Jet_Pt", Jet_Pt)
Cat10.Branch("Jet_Eta", Jet_Eta)
Cat10.Branch("Nu_BJet", Nu_BJet, "Nu_BJet/I")
Cat10.Branch("genweight", genweight, "genweight/F")
Cat10.Branch("puweight", puweight, "puweight/F")
"""


def LumiCheck(event):
    run = str(event.run)
    if run in Gfile:
        for start, end in Gfile[run]:
            if start <= event.luminosityBlock <= end:
                return True         
        return False
    else:
        return False

def MuScaleFactor (mu_charge, mu_pt, mu_eta, mu_phi, nTrack):
    scaleFactor = 1.0
    u1 = ROOT.gRandom.Rndm()
    u2 = ROOT.gRandom.Rndm()
    if "Run" in FileArg[1]:
        scaleFactor = rocCor.kScaleDT(mu_charge, mu_pt, mu_eta, mu_phi, 0, 0)
    else: 
        if mu_pt == GenLep1.Pt():
            scaleFactor = rocCor.kScaleFromGenMC(mu_charge, mu_pt, mu_eta, mu_phi, nTrack, GenLep1.Pt(), u1, 0, 0);
        if mu_pt == GenLep2.Pt():
            scaleFactor = rocCor.kScaleFromGenMC(mu_charge, mu_pt, mu_eta, mu_phi, nTrack, GenLep2.Pt(), u1, 0, 0);
        else:
            scaleFactor = rocCor.kScaleAndSmearMC(mu_charge, mu_pt, mu_eta, mu_phi, nTrack, u1, u2, 0, 0);   
    
    return scaleFactor 


def MuonSelection (mu_pt , mu_eta, mu_phi, mu_m, mu_iso, mu_charge, mu_id, nTrack, Tracker, Global):
    #if not Tracker: return False 
    #if not Global: return False 
    if not mu_id : return False 

    m = ROOT.TLorentzVector()
    mu = ROOT.TLorentzVector()
    m.SetPtEtaPhiM(mu_pt, mu_eta, mu_phi, mu_m)
    mu = m * MuScaleFactor(mu_charge, mu_pt, mu_eta, mu_phi, nTrack)

    if abs(mu.Pt()) < 20: return False 
    if abs(mu.Eta()) > 2.4: return False 
    if mu_iso > 0.15: return False
    Mu_Pt.push_back(mu.Pt())
    Mu_Eta.push_back(mu.Eta())
    Mu_Charge.push_back(mu_charge)
    Mu_Phi.push_back(mu.Phi())
    Mu_M.push_back(mu.M())
    return True

def ElecSelection (elec_pt, elec_eta, elec_phi, elec_m, elec_iso, elec_id, mu_p):
    ElP4 = ROOT.TLorentzVector()
    ElP4.SetPtEtaPhiM(elec_pt, elec_eta, elec_phi, elec_m)
    #for i, mu in enumerate (mu_p): 
    #    if mu.DeltaR(ElP4) < 0.4:
    #        return False

    if elec_pt < 20: return False   
    if abs(elec_eta) > 2.4: return False 
    if elec_iso > 0.12: return False 
    if elec_id < 3: return False
    
    El_Pt.push_back(elec_pt)
    El_Eta.push_back(elec_eta)
    #El_Charge.push_back(elec_charge)
    #El_Phi.push_back(elec.Phi())
    #El_M.push_back(elec.M())
    return True

def JetSelection (jet_pt, jet_eta, jet_phi, jet_m, jet_id, jet_b, mu_p):
    JET_LOOSE = (1<<0)
    if jet_pt < 30: return False  
    if abs(jet_eta) > 2.4: return False 
    if jet_id < 1: return False 
    JetP4 = ROOT.TLorentzVector()
    JetP4.SetPtEtaPhiM(jet_pt, jet_eta, jet_phi, jet_m)
    
    for i, mu in enumerate (mu_p):
        if mu.DeltaR(JetP4) < 0.4:
            return False      
       
    Jet_Pt.push_back(jet_pt)
    Jet_Eta.push_back(jet_eta)
    Jet_CSVV2.push_back(jet_b)
    Jet_M.push_back(jet_m)
    Jet_Phi.push_back(jet_phi)
    return True 

def BtaggedSelection (Jet_Pt, Jet_Eta, Jet_CSVV2):
    if Jet_Pt < 20: return False 
    if abs(Jet_Eta) > 2.4: return False 
    if Jet_CSVV2 < 0.848: return False 
    return True 

def mainloop (event):
      ### Object Selection ########################################################################################################################
        ### Muon Selection ###       
    Event_Total[0] = 1
    Event_Tot.Fill(0.5, Event_Total[0])
    NuMu = 0
    NuEl = 0
    NuJet = 0
    NuBJet = 0
        
    ### puWeight ###
    if hasattr(event,"Pileup_nTrueInt"):
        nvtx = int(getattr(event,"Pileup_nTrueInt"))
        puweight[0] = puWeight.getWeight(nvtx) if nvtx < hist_mc.GetNbinsX() else 1
    else: puweight[0] = 1

    ### Weights ###
    if "DoubleMuon" not in FileArg[1]:
        genweight[0] = event.genWeight 
        genweights.Fill(0.5, genweight[0])
        b_weight[0] = genweight[0] * puweight[0]
        weight.Fill(0.5, b_weight[0])

    ### Generated Lepton ###
        for i in range(event.nGenDressedLepton): 
            if abs(event.GenDressedLepton_pdgId[i]) != 13 : 
                break    
            
            bosonSample = False
            isfromBoson = False 
            for k in range(event.nGenPart):
                if (event.GenPart_genPartIdxMother[k] == 23 or event.GenPart_genPartIdxMother[k] == 25):
                    bosonSample = True 
                    isfromBoson = True 
                
            if isfromBoson == True:    
                if event.GenDressedLepton_pdgId[i] == 13:
                    GenLep1.SetPtEtaPhiM(event.GenDressedLepton_pt[i], event.GenDressedLepton_eta[i], event.GenDressedLepton_phi[i], event.GenDressedLepton_mass[i])
                else: 
                    GenLep2.SetPtEtaPhiM(event.GenDressedLepton_pt[i], event.GenDressedLepton_eta[i], event.GenDressedLepton_phi[i], event.GenDressedLepton_mass[i])
      
                        ### Muon Selection ###
    if event.nMuon > 0:
        for i in range(event.nMuon):
            if MuonSelection(event.Muon_pt[i], event.Muon_eta[i], event.Muon_phi[i], event.Muon_mass[i], event.Muon_pfRelIso04_all[i], event.Muon_charge[i], event.Muon_tightId[i], event.Muon_nTrackerLayers[i], event.Muon_trackerMu[i], event.Muon_globalMu[i]):
                NuMu += 1
            Nu_Mu[0] = NuMu    
    #if NuMu > 0:        
    #    return True   
        
    ### Muon TLorentzVector ###    
    Mu_P4 = []            
    for i in range(NuMu):
        MuP4 = ROOT.TLorentzVector()
        MuP4.SetPtEtaPhiM(Mu_Pt[i], Mu_Eta[i], Mu_Phi[i], Mu_M[i])
        Mu_P4.append(MuP4)

    ### Election Selection ###  
    if event.nElectron > 0:
        for k in range(event.nElectron):
            if ElecSelection(event.Electron_pt[k], event.Electron_eta[k], event.Electron_phi[k], event.Electron_mass[k], event.Electron_pfRelIso03_all[k], event.Electron_cutBased[k], Mu_P4):
                NuEl += 1 
        Nu_El[0] = NuEl

    ### Jet Selection ###
    if event.nJet > 0:
        for j in range(event.nJet):
            if JetSelection(event.Jet_pt[j], event.Jet_eta[j], event.Jet_phi[j], event.Jet_mass[j], event.Jet_jetId[j], event.Jet_btagCSVV2[j], Mu_P4):
                NuJet += 1
            Nu_Jet[0] = NuJet
    
        ## B-Tagged Jet Selection ###
        for l in xrange(NuJet):
            if BtaggedSelection(Jet_Pt[l], Jet_Eta[l], Jet_CSVV2[l]):
                NuBJet += 1
            Nu_BJet[0] = NuBJet
                        
        
### Muon With Opp Charge ###

    Charge = False
    for i in xrange(NuMu):   
        if ((Mu_Pt[0] > 20) or (Mu_Pt[i] > 20)): 
            if Mu_Charge[0] * Mu_Charge[i] < 0:
                Mu1.SetPtEtaPhiM(Mu_Pt[0],Mu_Eta[0],Mu_Phi[0],Mu_M[0])
                Mu2.SetPtEtaPhiM(Mu_Pt[i],Mu_Eta[i],Mu_Phi[i],Mu_M[i])
                Charge = True
                break
    
    Dilep_ = Mu1 + Mu2
    Dilep.SetPtEtaPhiM(Dilep_.Pt(),Dilep_.Eta(),Dilep_.Phi(),Dilep_.M())
    
    if Charge == True and abs(Dilep.M())>20 :
        return True
    Step[0] = 1
    cutFlow.Fill(1)


    if 76 < Dilep.M() and Dilep.M() < 106 :     
        return False   
    Step[0] = 2        
    cutFlow.Fill(2)
        
                        
    if event.MET_pt > 40 :
        return True
    Step[0] = 3
    cutFlow.Fill(3)

    
    if len(Nu_Jet) > 1 :
        return True
    Step[0] = 4
    cutFlow.Fill(4)

    if len(Nu_BJet) > 0 :
        return True
    Step[0] = 5
    cutFlow.Fill(5)
                
      ### Event Selection ############################################################################################################################
        
        
         ### Muon Triggers ###
        #if not event.HLT_IsoMu24 and not event.HLT_IsoTkMu24:
           # Step6[0] = False        
        #    return True   
        #if Step5[0] == True:  
        #Step[0] = 6    
       
      #  if Step5[0] == True and Step4[0] == True and Step3[0] == True and Step2[0] == True and Step1[0] == True:
      #  NewStep[0] = 6
      
    Event_No[0] = 1    
    return True 

def reset_all():
        Event_No[0] = 0
        Event_Total[0] = 0
        Nu_Mu[0] = 0
        Nu_El[0] = 0
        Nu_Jet[0] = 0
        Nu_BJet[0] = 0
        Nu_NonBJet[0] = 0
        genweight[0] = 0
        puweight[0] = 0
        b_weight[0] = 0
        Step[0] = 0
        NewStep[0] = 0

        Mu_Pt.clear()
        Mu_Eta.clear()
        Mu_Charge.clear()
        Mu_Phi.clear()
        Mu_M.clear()

        El_Pt.clear()
        El_Eta.clear()
        El_Charge.clear()
        El_Phi.clear()
        El_M.clear()


        Jet_Pt.clear()
        Jet_Eta.clear()
        Jet_CSVV2.clear()
        Jet_M.clear()
        Jet_Phi.clear()
        
        Dilep.SetPtEtaPhiM(0.,0.,0.,0.)
        Mu1.SetPtEtaPhiM(0.,0.,0.,0.)
        Mu2.SetPtEtaPhiM(0.,0.,0.,0.)
        GenLep1.SetPtEtaPhiM(0.,0.,0.,0.)
        GenLep2.SetPtEtaPhiM(0.,0.,0.,0.)

GJsonF = open("/cms/scratch/jdj0715/nanoAOD/src/nano/analysis/data/GoldenJson.txt")
Gfile = json.load(GJsonF)

for i,Nfile in enumerate(FileArg[2:]):
    CurFile = TNetXNGFile(Nfile)
    Tree = CurFile.Get("Events")
    
    for ive, event in enumerate(Tree):
       reset_all() 
       Step[0] = 0
       cutFlow.Fill(0)
       if mainloop(event):  
            ALL.Fill()    

f.Write()
f.Close()
GJsonF.close()
