#!/bin/bash

dir=`pwd`
cd $dir
#root -l -b < x_test.C >& result.log &
root -l -b < x_para.C >& result.log &
