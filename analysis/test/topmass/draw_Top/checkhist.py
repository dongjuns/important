import ROOT
rdfile = "/cms/scratch/daniel/nanoAOD/src/nano/analysis/test/Results/Nano_TrigEnd/results_merged/tth2mu_SingleMuon_Run2016.root"
rdFile = ROOT.TFile(rdfile)
tree = rdFile.Get("nEvent")
step_0=0
step_1=0
step_2=0
step_3=0
step_4=0
step_5=0
step_6=0
print "running"
for ive, event in enumerate(tree):
    if event.Dilep.M() > 60: 
        if event.Step == 0 :
            step_0 = step_0 + 1
        if event.Step == 1 :
            step_1 = step_1 + 1
        if event.Step == 2 :
            step_2 = step_2 + 1
        if event.Step == 3 :
            step_3 = step_3 + 1
        if event.Step == 4 :
            step_4 = step_4 + 1
        if event.Step == 5 :
            step_5 = step_5 + 1
        if event.Step == 6 :
            step_6 = step_6 + 1

fout = open("stepchechM60.txt", "w")

print>>fout, "Step 0:  ", step_0      
print>>fout, "Step 1:  ", step_1      
print>>fout, "Step 2:  ", step_2      
print>>fout, "Step 3:  ", step_3      
print>>fout, "Step 4:  ", step_4      
print>>fout, "Step 5:  ", step_5      
print>>fout, "Step 6:  ", step_6      

fout.close()
