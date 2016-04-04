#!/bin/bash

## This script has to be executed from the
## sequence run output directory
## (there were the Data/Intensities dir)

echo "==============";
echo "## currentDirectory"
echo $PWD
echo "## conversionLog";
cat ./conversionLog.txt; 
echo "## conversionError tail"
cat ./conversionError.txt | tail -1; 
echo "==============";

