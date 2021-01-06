#!/bin/bash

input=TestRun
COPYING=0
CLEANING=0




function cleanup() {
        echo "cleanup"
}


function copyback() {
        if [[ $COPYING -eq 0 ]]; then
                COPYING=1
                echo "Tar everything"
                tar cfz $input.tgz ${input}*
                echo "done with tar"
        fi
}



EV=0

trap 'copyback; cleanup' SIGUSR2 SIGINT SIGTERM


eval ./RunIBM $input


EV=$?

copyback





