#!/bin/sh
#$ -S /bin/bash 
source /cvmfs/cms.cern.ch/cmsset_default.sh 
cd /cms/ldap_home/jongho/CMSSW/CMSSW_9_3_7/src 
eval `scramv1 runtime -sh` 
##cd /cms/ldap_home/jongho/L1PixelTrigger/L1PixEle-cmssw9/Signal/L1PixEle-Eff/Results/SE_PU200_719445/Job_1/ 
cd /cms/ldap_home/jongho/L1PixelTrigger/TrkIso_test
rm -rf test_C.so 
rm -rf test_C.d 
rm -rf result.log 
root -l -b < x_test.C >& process.log 
rm -rf test_C.so 
rm -rf test_C.d 
