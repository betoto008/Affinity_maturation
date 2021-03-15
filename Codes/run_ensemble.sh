#!/bin/sh

#  run_ensemble.sh
#  
#
#  Created by Roberto Moran Tovar on 15.03.21.
#  

c++ Dynamics_ensemble.cpp -o ensemble.x

if [ $? -eq '0' ]       ## check the compiler return, run only on success
then
    ./ensemble.x 12 10000 19 14 0.001 1000
else
    printf "\nAn error occurred, executable not called\n\n" >&2
fi
