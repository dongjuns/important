import ROOT, os, getopt, sys, array, math, glob, json
from ROOT import * 
from array import array


### Make TTREE ### 
FileArg = sys.argv
print FileArg
#tempdir = FileArg[1]
Dirname = "%s/src/nano/analysis/test/Results/Nano_Diff_Trig/"%(os.environ['CMSSW_BASE'])
if not os.path.isdir(Dirname):
    os.makedirs(Dirname)

temp = FileArg[2].split('/').pop()
cattree = Dirname+temp
envName = "%s" %os.environ['CMSSW_BASE']

GJsonF = open("%s/src/nano/analysis/data/Cert_271036-284044_13TeV_PromptReco_Collisions16_JSON.txt"%os.environ['CMSSW_BASE'])
Gfile = json.load(GJsonF)
GJsonF.seek(0)
json_hold = GJsonF.read()

lumiMap = ROOT.std.map('unsigned int', 'std::vector<std::array<unsigned int, 2>>')()

for i in Gfile :
    lumiVector = ROOT.std.vector('std::array<unsigned int, 2>')()
    for j in range(len(Gfile[i])) :
        lumiArray = ROOT.std.array('unsigned int', 2)()
        lumiArray[0] = int(Gfile[i][j][0])
        lumiArray[1] = int(Gfile[i][j][1])
        lumiVector.push_back(lumiArray)
    lumiMap[int(i)] = lumiVector

ROOT.gROOT.LoadMacro("%s/src/nano/analysis/plugins/tth2mu.cc+"%os.environ['CMSSW_BASE'])
analyzer = ROOT.TTH2MuAnalyzer(cattree, envName)
analyzer.LoadLumiMap(lumiMap)

isMC = False
for i,Nfile in enumerate(FileArg[2:]):
    if "Run" not in FileArg[1]:
        isMC = True
    else:
        isMC = False
    analyer.analyze(Nfile, isMC)
analyzer.EndOfAnalyze()
