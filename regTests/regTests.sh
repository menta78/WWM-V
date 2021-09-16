#!/bin/bash
# author: Lorenzo Mentaschi



# setting the environment variables
. set_environ.sh



# checking python
if [ -z $SCWWM_PY ]; then
  #SCWWM_PY is empty. Defaulting to python3
  SCWWM_PY=python3
fi

$SCWWM_PY -c '
import sys
if sys.version_info < (3,5):
  raise Exception()
'
if [ $? != 0 ]; then
  echo 'python should be at least 3.5'
  exit
fi
  
$($SCWWM_PY -c "import numpy, netCDF4")
if [ $? != 0 ]; then 
  echo "numpy and netCDF4 need to be installed in python"
  exit
fi



# launching the tests
if [[ $1 == '-a' ]]; then
  echo "launching all the tests"
  $SCWWM_PY -m unittest discover
elif [[ $1 == '-k' ]]; then
  echo "launching the tests corresponding to the given name patter"
  ptrn=$(echo $2 | sed "s:/: :g")
  $SCWWM_PY -m unittest -k $ptrn
else
  echo "examples:"
  echo "./regTests.sh -a #runs all the tests"
  echo "./regTests.sh -k myTest* #runs the tests with name corresponding to the given pattern"
fi
  

