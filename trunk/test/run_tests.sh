#!/bin/sh
#
# Run all the tests using different modules
#
# Note: you can add an arbitrary argument to stop testing test_20_rpy
#
echo ===================================================================
echo 
echo   RUNNING TESTS FOR
echo
echo         binary libraries
echo
echo 
SIMUALLELETYPE=binary
export SIMUALLELETYPE
if python run_tests.py $1
then
  echo SUCCEED
else
  echo BINARY-ALLELE MODULES FAILED
  exit 1
fi

sleep 10
echo ===================================================================
echo 
echo   RUNNING TESTS FOR
echo
echo         short libraries
echo
echo 
SIMUALLELETYPE=short
export SIMUALLELETYPE
if python run_tests.py $1
then
  echo SUCCEED
else
  echo SHORT-ALLELS MODULES FAILED
  exit 1
fi
sleep 10

echo ===================================================================
echo 
echo   RUNNING TESTS FOR
echo
echo         long libraries
echo
echo 
SIMUALLELETYPE=long
export SIMUALLELETYPE
if python run_tests.py $1
then
  echo SUCCEED
else
  echo LONG-ALLELE MODULES FAILED
  exit 1
fi

