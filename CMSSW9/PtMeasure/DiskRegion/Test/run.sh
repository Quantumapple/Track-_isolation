#!/bin/bash

dir=`pwd`
cd $dir

mkdir Plots

root -l -b < x_run.C >& result.log &
