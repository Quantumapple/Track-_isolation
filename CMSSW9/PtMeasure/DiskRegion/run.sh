#!/bin/bash

##source /cvmfs/cms.cern.ch/cmsset_default.sh 
##
##cd /cms/ldap_home/jongho/CMSSW/CMSSW_9_3_7/src/ 
##eval `scramv1 runtime -sh`
##
##cd /cms/ldap_home/jongho/L1PixelTrigger/Parameters/Pt 
##root -l -b < x_test.C >& result.log  

dir=`pwd`
cd $dir
root -l -b < x_test.C >& result.log & 
