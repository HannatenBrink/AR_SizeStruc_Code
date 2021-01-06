#!/bin/bash

input=$(awk -v jindex=1 'NR==jindex' IBMTest_runlist)

PROGNAME = RUNIBM
Echo $input
./$PROGNAME $input