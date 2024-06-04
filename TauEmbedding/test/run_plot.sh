#!/bin/bash

for i in `seq 1 1 500`
do
    #python3 plotTausFromMINIAOD.py $i
    python3 plotTausBoosted.py $i
done
