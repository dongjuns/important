#!/usr/bin/env python 
# catGetDatasetInfo v7-4-4 # to make dataset lists
# sed -i 's/^\/store/root:\/\/cms-xrdr.sdfarm.kr:1094\/\/xrd\/store/g' *

import os, json, array, sys
import numpy as np
from math import ceil       
username = os.environ['USER']

analysis = 'TT'
#analysis = 'TtbarDiLeptonAnalyzer'
pythonCfg = 'analyzer_Top.py'
#analysis=analysis+'Silver'
RunFiles = [
             # 'WMinusH_HToMuMu',
             # 'WPlusH_HToMuMu',
             # 'ZH_HToMuMu',
             # 'VBF_HToMuMu',
             # 'GG_HToMuMu',
             # "WWTo2L2Nu",
             # "WZTo3LNu_amcatnlo",
             # "WZTo2LQQ",
             # "ZZTo2L2Nu",
             # "ZZTo2L2Q",
             # "ZZTo4L",
             # "WWW",
             # "WWZ",
             # "WZZ",
             # "ZZZ",
#              "ttZToLLNuNu",
#              "ttWToLNu",
             # "SingleTop_tW_noHadron",
             # "SingleTbar_tW_noHadron",
             # "SingleTop_tW",
             # "SingleTbar_tW",
#              "TTJets_DiLept",
#              "TTJets_DiLept_Tune4",
              'TT_powheg',
              #'TTJets_aMC', 
              #'DYJets',
	      # 'TT_powheg',
#              'DYJets_MG_10to50',
#              'DYJets_MG2',
#              'DYJets_2J',
#              'DYJets_1J',
#              'DYJets_0J',
#              'DYJets_10to50', 
              'DoubleMuon_Run2016B',
              'DoubleMuon_Run2016C',
              'DoubleMuon_Run2016D',
              'DoubleMuon_Run2016E',
              'DoubleMuon_Run2016F',
              'DoubleMuon_Run2016G',
              'DoubleMuon_Run2016H',
             # 'SingleMuon_Run2016H_v3',
              ]
datadir = '/cms/scratch/daniel/nanoAOD/src/nano/analysis/data/dataset/' 
#ersion = os.environ["CMSSW_VERSION"]


for i in RunFiles:
    datasetName = i
    fileList = datadir + 'dataset_' + datasetName + '.txt'
    jobName = analysis+'_'+datasetName 

    Dirname = "/cms/scratch/jdj0715/nanoAOD/src/nano/analysis/test/topmass/batch_Top/%s/"%jobName
    DirnameJDS = "/cms/scratch/jdj0715/nanoAOD/src/nano/analysis/test/topmass/batch_Top/"
    if os.path.isdir(Dirname):
        print "ERROR: output directory already existing."
        sys.exit()
    else: os.makedirs(Dirname)

    files = np.array([])
    for f in open(fileList).readlines():
        f = f.strip()
        f = f.strip('\',"')
        if len(f) < 5: continue
        if '#' == f[0] or '.root' != f[-5:]: continue
        files = np.append(files,[f])
    nFiles = len(files)     
    maxFiles = 10
    nSection = int(ceil(1.0*nFiles/maxFiles))
    count = 0
    for section in range(nSection):
        begin = section*maxFiles
        end = min(begin + maxFiles, nFiles)
        FileNames = files[begin:end]
        FileNamesStr = " ".join(str(i) for i in FileNames)

        print "@@ Writing run script..."
        jds = "%ssubmit.jds" %Dirname 
        fout = open(jds, "w")
        print>>fout, "# Job description file for condor job"
        print>>fout, """executable = /cms/scratch/jdj0715/nanoAOD/src/nano/analysis/test/topmass/batch_Top/ttbar.sh
universe   = vanilla

log = condor.log

getenv     = True
should_transfer_files = YES
when_to_transfer_output = ON_EXIT
output = %s/job_%d.log
error = %s/job_%d.err
transfer_input_files = NanoAOD
queue""" % (Dirname, count, Dirname, count)
        fout.close()
        count += 1 
        #jobName = analysis+'_'+datasetName
        subBatch = "condor_submit -batch-name %s -append 'arguments=%s %s' %s" %(datasetName ,datasetName,FileNamesStr, jds)
        #print createbatch
        print subBatch 
            
        os.system(subBatch)
