import ROOT, os, getopt, sys, array, math, glob, json
from ROOT import * 
from array import array

#pu weight
ROOT.gROOT.LoadMacro("%s/src/nano/analysis/scripts/WeightCalculatorFromHistogram.cc+" %os.environ['CMSSW_BASE'])
pufile_mc="%s/src/nano/analysis/data/pu_root/pileup_profile_Spring16.root" %os.environ["CMSSW_BASE"]
fmc = ROOT.TFile(pufile_mc)
pufile_data="%s/src/nano/analysis/data/pu_root/PileupData_GoldenJSON_Full2016.root" %os.environ['CMSSW_BASE']
fmcrd = ROOT.TFile(pufile_data)

hist_mc = fmc.Get("pu_mc")
hist_mc.SetDirectory(0)
hist_data = fmcrd.Get("pileup")
hist_data.SetDirectory(0)
puWeight = ROOT.WeightCalculatorFromHistogram(hist_mc, hist_data, True, True, False)

#TTree
FileArg = sys.argv
tempdir = FileArg[1]
Dirname ="%s/src/nano/analysis/test/ttbar_Results/nano/%s/" %(os.environ['CMSSW_BASE'],temdir)
if not os.path.isdir(Dirname):
    os.makedirs(Dirname)
temp = FileArg[2].split('/').pop()
cattree = Dirname +temp



#variables
Mu_pt = ROOT.std.vector('float')()
Mu_eta = ROOT.std.vector('float')()
Mu_charge = ROOT.std.vector('float')()
Mu_phi = ROOT.std.vector('float')()
Mu_mass = ROOT.std.vector('float')()

El_pt = ROOT.std.vector('float')()
El_eta = ROOT.std.vector('float')()
El_phi = ROOT.std.vector('float')()
El_mass = ROOT.std.vector('float')()
El_charge = ROOT.std.vector('float')()

Jet_pt = ROOT.std.vector('float')()
Jet_eta = ROOT.std.vector('float')()
Jet_mass = ROOT.std.vector('float')()
Jet_phi = ROOT.std.vector('float')()


lep1 = ROOT.TLorentzVector()
lep2 = ROOT.TLorentzVector()


#Branch

#Scale Factor
def muScaleFactor (mu_pt, mu_eta, mu_phi, mu_charge, mu_ntrack):
    scaleFactor = 1.0

#Object selection

def selMuon(mu_pt, mu_eta, mu_phi, mu_mass, mu_charge, mu_id, mu_reliso):
    if mu_pt  < 20 : return False
    if abs(mu_eta) > 2.4 : return False
    if not mu_id : return False
    if mu_reliso > 0.15 : return False

    mu = ROOT.TLorentzVector()
    mu.SetPtEtaPhiM(mu_pt, mu_eta, mu_phi, mu_mass)

    Mu_pt.push_back(mu.Pt())
    Mu_eta.push_back(mu.Eta())
    Mu_charge.push_back(mu_charge)
    Mu_phi.push_back(mu.Phi())
    Mu_mass.push_back(mu.M())
    return True

def selElec(el_pt, el_eta, el_phi, el_mass, el_charge, el_id, el_detasc, el_reliso):
    if el_pt < 20 : return False
    if abs(el_eta) > 2.4 : return False
    if el_id < 3: return False
    el_scEta = el_eta + el_detasc
    if abs(el_scEta) > 1.4442 and abs(el_scEta) < 1.566: return False
    #if el_reliso > 0.0571: return False
    elec = ROOT.TLorentzVector()
    elec.SetPtEtaPhiM(el_pt, el_eta, el_phi, el_mass)

    El_pt.push_back(elec.Pt())
    El_eta.push_back(elec.Eta())
    El_charge.push_back(el_charge)
    El_phi.push_back(elec.Phi())
    El_mass.push_back(elec.M())
    return True

def selJet(jet_pt, jet_eta, jet_phi, jet_mass, jet_id):
    if jet_pt < 30 : return False
    if abs(jet_eta) > 2.4 : return False
    if jet_id < 1 : return False
    jet = ROOT.TLorentzVector()
    jet.SetPtEtaPhiM(jet_pt, jet_eta, jet_phi, jet_mass)

    Jet_pt.push_back(jet.Pt())
    Jet_eta.push_back(jet.Eta())
    Jet_phi.push_back(jet.Phi())
    Jet_mass.push_back(jet.M())

    hasOverLap = False
    for lep in recoleps:
        leptlv = ROOT.TLorentzVector()
        leptlv.SetPtEtaPhiM(lep.Pt(), lep.Eta(), lep.Phi(), lep.M())
        if (jet).DeltaR(leptlv) < 0.4 : hasOverLap = True
    if (hasOverLap) : return False

    return True




#for cmeson
'''
def pickD0(cme_pt, cme_eta, cme_phi, cme_mass, cme_pdgid, cme_lxy, cme_l3d):
    for i in range(event.ncmeson):
        if event.cmeson_pdgId[i] != 421: continue
        if event.cmeson_lxy[i] < 0.1: continue
        if event.cmeson_l3D[i] < 0.2: continue

        #d0 = ROOT.TParticle()
        d0 = ROOT.TLorentzVector()
        d0.SetPtEtaPhiM(event.cmeson_pt[i], event.cmeson_eta[i], event.cmeson_phi[i], event.cmeson_mass[i])
        d0s.append(d0)
    return d0s

def pickDstar(event):
    dstars = []
    for i in range(event.ncmeson):
        if event.cmeson_pdgId[i] != 413: continue
        if event.cmeson_lxy[i] < 0.1: continue
        if event.cmeson_l3D[i] < 0.2: continue

        dstar = ROOT.TLorentzVector()
        dstar.SetPtEtaPhiM(event.cmeson_pt[i], event.cmeson_eta[i], event.cmeson_phi[i], event.cmeson_mass[i])
        dstars.append(dstar)
    return dstars

def pickJpsi(event):
    jpsis = []
    for i in range(event.ncmeson):
        if event.cmeson_pdgId[i] != 443: continue

        jpsi = ROOT.TLorentzVector()
        jpsi.SetPtEtaPhiM(event.cmeson_pt[i], event.cmeson_eta[i], event.cmeson_phi[i], event.cmeson_mass[i])
        jpsis.append(jpsi)
    return jpsis

def pickdiffMass(event):
    diffs = []
    for i in range(event.ncmeson):
        diff = ROOT.TLorentzVector()
        diff.SetPtEtaPhiM(event.cmeson_pt[i], event.cmeson_eta[i], event.cmeson_phi[i], event.cmeson_diffMass[i])
        diffs.append(diff)
'''
tstack = ROOT.THStack("dilep mass", "dilep mass")

datalumi = 38*1000

#ttbardir = "/xrootd/store/group/nanoAOD/run2_2016v1/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/171117_140513/0000/"
#dydir10to50 = "/xrootd/store/group/nanoAOD/run2_2016v1/DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/run2_2016MC_NANO_RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/171111_220746/0000/"
#dydir1 = "/xrootd/store/group/nanoAOD/run2_2016v1/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/run2_2016MC_NANO_RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext2-v1/171111_220522/0000/"

ttbardir = "/xrootd/store/group/nanoAOD/run2_2016v2/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/171130_031132/0000/"
dydir10to50 = "/xrootd/store/group/nanoAOD/run2_2016v2/DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/171130_030955/0000/"
dydir1 = "/xrootd/store/group/nanoAOD/run2_2016v2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext2-v1/171130_030907/0000/"

samples = [ttbardir, dydir10to50, dydir1]
#samples = [ttbar]
xsec = [831.76, 18610, 6025.2]
#samples = [ttbardir]
#xsec = [831.76]
#samples = [ttbardir, dydir]



hlist = []
for j, sampledir in enumerate(samples):
    h = copy.deepcopy(ROOT.TH1D("dilep mass", "dilep mass", 60, 20, 320))
    #h = copy.deepcopy(ROOT.TH1D("jpsi mass", "jspi mass", 120, 0, 6))
    #h = copy.deepcopy(ROOT.TH1D("Diffmass", "Diffmass", 120, 0, 6))
    nevents =0
    filelist = [l for l in os.listdir(sampledir) if "root" in l]
    for i, fileName in enumerate(filelist):
        print fileName
        inFile = ROOT.TFile(sampledir+fileName)
        events = inFile.Get("Events")
        nevents += events.GetEntries()

        if i == 3 : break
        recoleps=[]
        for iev, event in enumerate(events):
            Mu_pt.clear()
            Mu_eta.clear()
            Mu_charge.clear()
            Mu_phi.clear()
            Mu_mass.clear()
    
            El_pt.clear()
            El_eta.clear()
            El_charge.clear()
            El_phi.clear()
            El_mass.clear()
    
            Jet_pt.clear()
            Jet_eta.clear()
            Jet_phi.clear()
            Jet_mass.clear()
            N_recolep = 0
            N_muon = 0
            N_elec = 0
            N_jet = 0
            idxs =[]

            ###muon selection
            for i in range(event.nMuon):
                if selMuon(event.Muon_pt[i], event.Muon_eta[i], event.Muon_phi[i], event.Muon_mass[i], event.Muon_charge[i], event.Muon_tightId[i], event.Muon_pfRelIso04_all[i]):
                    N_muon += 1

            ###electron selection
            for j in range(event.nElectron):
                if selElec(event.Electron_pt[j], event.Electron_eta[j], event.Electron_phi[j], event.Electron_mass[j], event.Electron_charge[j], event.Electron_cutBased[j], event.Electron_deltaEtaSC[j], event.Electron_pfRelIso03_all[j]):
                    N_elec += 1

            ###recoleps
            N_recolep = N_muon + N_elec
            if N_recolep != 2 : continue
    
            if N_muon == 2 :
                channel = "MuMu"
                lep1.SetPtEtaPhiM(Mu_pt[0], Mu_eta[0], Mu_phi[0], Mu_mass[0])
                lep2.SetPtEtaPhiM(Mu_pt[1], Mu_eta[1], Mu_phi[1], Mu_mass[1])
                mulcharge = Mu_charge[0] * Mu_charge[1]
            if N_muon == 1 and N_elec == 1 :
                channel = "MuEl"
    
                lep1.SetPtEtaPhiM(Mu_pt[0], Mu_eta[0], Mu_phi[0], Mu_mass[0])
                lep2.SetPtEtaPhiM(El_pt[0], El_eta[0], El_phi[0], El_mass[0])
                mulcharge = Mu_charge[0] * El_charge[0]
    
            if N_elec == 2 :
                channel = "ElEl"
                lep1.SetPtEtaPhiM(El_pt[0], El_eta[0], El_phi[0], El_mass[0])
                lep2.SetPtEtaPhiM(El_pt[1], El_eta[1], El_phi[1], El_mass[1])
                mulcharge = El_charge[0] * El_charge[1]
    
            dilep = lep1 + lep2
            recoleps = [lep1,lep2]

            ###jet selection
            for k in range(event.nJet):
                if selJet(event.Jet_pt[k], event.Jet_eta[k], event.Jet_phi[k], event.Jet_mass[k], event.Jet_jetId[k]):
                    N_jet += 1
                    idxs.append(k)
    
            ### Event Selection ###    
    
            #step1
            if abs(dilep.M()) < 20 : continue
            if mulcharge > 0 : continue
            step1 = True

            #step2        
            if channel != "MuEl" and 76 < abs(dilep.M()) < 106 : continue
            step2 = True
    
            #step3
            if channel != "MuEl" and event.MET_pt < 40 : continue
            step3 = True
    
            #step4
            if N_jet < 2 : continue
            step4 = True

            #step5
            if not any([event.Jet_btagCSVV2[idx]>0.8484 for idx in idxs]): continue
            step5 = True

        inFile.Close()
            
    h.SetTitle(sampledir)
    scale = datalumi*xsec[j]/float(nevents)
    print scale
    print nevents
    h.Scale(scale)
    hlist.append(h)


outFile = ROOT.TFile("dilep.root", "RECREATE")
outFile.cd()

for i, h in enumerate(hlist):
    h.SetLineColor(2+i*2)
    h.SetFillColor(2+i*2)
    if "DY" in h.GetTitle():
        h.SetLineColor(4)
        h.SetFillColor(4)
    tstack.Add(h)
    h.Write()

tstack.Write()

canv = ROOT.TCanvas()
#canv.SetLogy()
tstack.Draw("hist")
#rdhist.Draw("samee1")
canv.Print("dilep.png")

outFile.Write()
outFile.Close()



