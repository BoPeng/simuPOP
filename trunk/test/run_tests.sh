#!/bin/sh
echo ===================================================================
echo 
echo   RUNNING TESTS FOR
echo
echo         binary libraries
echo
echo 
SIMUALLELETYPE=binary
export SIMUALLELETYPE
if python run_tests.py
then
  echo SUCCEED
else
  echo BINARY-ALLELE MODULES FAILED
  exit 1
fi

echo ===================================================================
echo 
echo   RUNNING TESTS FOR
echo
echo         short libraries
echo
echo 
SIMUALLELETYPE=short
export SIMUALLELETYPE
if python run_tests.py
then
  echo SUCCEED
else
  echo SHORT-ALLELS MODULES FAILED
  exit 1
fi

echo ===================================================================
echo 
echo   RUNNING TESTS FOR
echo
echo         long libraries
echo
echo 
SIMUALLELETYPE=long
export SIMUALLELETYPE
if python run_tests.py
then
  echo SUCCEED
else
  echo LONG-ALLELE MODULES FAILED
  exit 1
fi

