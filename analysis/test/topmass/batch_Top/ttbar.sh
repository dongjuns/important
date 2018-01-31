#!/bin/bash 

RESULT_DIR="NanoAOD"
DIR=`hostname`

/bin/mkdir -p $RESULT_DIR/$DIR

echo "hostname" >> $RESULT_DIR/$DIR/`hostname`.txt
echo "date" >> $RESULT_DIR/$DIR/`hostname`.txt
source /cvmfs/cms.cern.ch/cmsset_default.sh >> $RESULT_DIR/$DIR/`hostname`.txt
python /cms/scratch/jdj0715/nanoAOD/src/nano/analysis/python/analyzer_Top.py $@
echo "Done." >> $RESULT_DIR/$DIR/`hostname`.txt

