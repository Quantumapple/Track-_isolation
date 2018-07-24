#!/bin/bash

###path=/u/user/quantumapple/TrkIso/Testing
path=/u/user/quantumapple/L1Pixel_KNU/TrkEleIso/Pt_phi_fit/Testing
for i in 0 1 2 3 4 5 
    do
	if [ $i = 0 ]; 
	then directory="Bsl1_l1l2"
	   if [ ! -d $directory ]; 
	   then mkdir $directory
	   fi
	   cd $directory
           mkdir Plots
	   cp $path/test.C . 
	   cp $path/test.h .
	   exec 1>x_test.C
           echo ".L test.C()" 
	   echo "test a"
	   echo "a.Loop($i)"

           root -l -b < x_test.C >& result.log &
	   cd ../
	fi
	
	if [ $i = 1 ]; 
	then directory="Bsl1_l1l3"
	   if [ ! -d $directory ]; 
	   then mkdir $directory
	   fi
	   cd $directory
           mkdir Plots
	   cp $path/test.C . 
	   cp $path/test.h .
	   exec 1>x_test.C
           echo ".L test.C()" 
	   echo "test a"
	   echo "a.Loop($i)"

           root -l -b < x_test.C >& result.log &
	   cd ../
	fi
	
	if [ $i = 2 ]; 
	then directory="Bsl1_l1l4"
	   if [ ! -d $directory ]; 
	   then mkdir $directory
	   fi
	   cd $directory
           mkdir Plots
	   cp $path/test.C . 
	   cp $path/test.h .
	   exec 1>x_test.C
           echo ".L test.C()" 
	   echo "test a"
	   echo "a.Loop($i)"

           root -l -b < x_test.C >& result.log &
	   cd ../
	fi
	
	if [ $i = 3 ]; 
	then directory="Bsl2_l2l3"
	   if [ ! -d $directory ]; 
	   then mkdir $directory
	   fi
	   cd $directory
           mkdir Plots
	   cp $path/test.C . 
	   cp $path/test.h .
	   exec 1>x_test.C
           echo ".L test.C()" 
	   echo "test a"
	   echo "a.Loop($i)"

           root -l -b < x_test.C >& result.log &
	   cd ../
	fi
	
	if [ $i = 4 ]; 
	then directory="Bsl2_l2l4"
	   if [ ! -d $directory ]; 
	   then mkdir $directory
	   fi
	   cd $directory
           mkdir Plots
	   cp $path/test.C . 
	   cp $path/test.h .
	   exec 1>x_test.C
           echo ".L test.C()" 
	   echo "test a"
	   echo "a.Loop($i)"

           root -l -b < x_test.C >& result.log &
	   cd ../
	fi
	
	if [ $i = 5 ]; 
	then directory="Bsl3_l3l4"
	   if [ ! -d $directory ]; 
	   then mkdir $directory
	   fi
	   cd $directory
           mkdir Plots
	   cp $path/test.C . 
	   cp $path/test.h .
	   exec 1>x_test.C
           echo ".L test.C()" 
	   echo "test a"
	   echo "a.Loop($i)"

           root -l -b < x_test.C >& result.log &
	   cd ../
	fi
    done

