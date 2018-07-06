def makeBatchConfigFile( job_dir):

    config='#!/bin/sh'
    config+='\n'
    config+='#$ -S /bin/bash \n'
    ##config+='cd /share/apps/root_v5-34-32/root/ \n'
    ##config+='cd /u/user/moon/cmssw/ \n'
    ##config+='cd /cms/ldap_home/jongho \n'
    ##config+='. bin/thisroot.sh \n'
    ##config+='./root.sh \n'
    config+='source /cvmfs/cms.cern.ch/cmsset_default.sh \n'
    config+='cd /cms/ldap_home/jongho/CMSSW/CMSSW_9_3_7/src \n'
    config+='eval `scramv1 runtime -sh` \n'
    ##config+='cd - \n'
    config+='cd ' + job_dir +  ' \n'
    config+='rm -rf test_C.so \n'
    config+='rm -rf test_C.d \n'
    config+='rm -rf result.log \n'
    config+='root -l -b < x_test.C \n'
    config+='rm -rf test_C.so \n'
    config+='rm -rf test_C.d \n'

    return config
