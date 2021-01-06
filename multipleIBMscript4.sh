#!/bin/bash

RUNLIST=IBMTest_runlist
PROGNAME=RunIBM
MAILADDRESS="hanna.tenbrink@eawag.ch"
WALLTIME=00:05
NCORES=1

#How many runs are there in the runlist file#
NR=$(wc -l < "$RUNLIST")

echo $WALLTIME

for (( i=1; i<=$NR; i++ ))
do
   input=$(awk -v jindex=$i 'NR==jindex' $RUNLIST)
   Echo $input
   ./$PROGNAME "$input" &
done






