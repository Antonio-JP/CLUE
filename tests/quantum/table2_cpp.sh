#!/bin/bash
#####################################################################################################################
###
### TESTS FOR DATA ON TABLE 1
###
#####################################################################################################################
MY_SCRIPT=./CLUE_main_script

## This are the different types of experiments
NEXPERIMENTS=17
EXPERIMENTS=("ae" "dj" "ghz" "graphstate" "hhl" "portfolioqaoa" "portfoliovqe" "pricingcall" "pricingput" "qft" "qnn" "qpeexact" "qpeinexact" "qwalk" "tsp" "vqe" "wstate")
OBSERVABLE=("0")
## This are the different types of approaches
TESTS="DDSIM_ALONE DDSIM CLUE"
REPEATS=5
MIN=3
MAXSIZES=(10 30 16 15 11 10 10 15 15 27 10 11 10 11 9 10 10)
MAX=27
TIMEOUT=500s

cd ../../cpp/clue/build/apps

# AE: Stops naturally at 10 qubits
# DJ: CLUE stops at 16 due to memory issues (surprising) -- We may increase the size for DDSIM
# GHZ: Stops naturally at 16 qubits
# GRAPHSTATE: Stops naturally at 15 qubits
# HHL: We only have some circuits up to 11 qubits
# PORTFOLIOQAOA: Stops naturally at 10 qubits
# PORTFOLIOVQE: Stops naturally at 10 qubits
# PRICINGCALL: we need to run these cases from 11 to 15
# PRIVINGPUT: we need to run these cases from 11 to 15
# QFT: CLUE stops at 13 qubits. Then DDSIM+CLUE stops at 28
# QNN: Stops naturally at 10 qubits
# QPEEXACT: Stops naturally at 11 qubits
# QPEINEXACT: Stops naturally at 10 qubits
# QWALK: Stops naturally at 11 qubits
# TSP: We only have some circuits up to 9 qubits
# VQE: Stops naturally at 10 qubits
# WSTATE: Stops naturally at 10 qubits

for size in $(seq $MIN 1 $MAX)
do
    for repeat in $(seq 1 1 $REPEATS)
    do
        n=0;
        while [ $n -lt $NEXPERIMENTS ]
        do
            if [ ${MAXSIZES[n]} -ge $size ]
            then
                experiment=${EXPERIMENTS[n]};
                for observable in $OBSERVABLE
                do
                    for test in $TESTS
                    do
                        echo timeout $TIMEOUT $MY_SCRIPT $experiment -m $size -M $size -r 1 -t $test -obs $observable
                    done
                done
            fi
            n=$(( $n + 1 ))
        done
    done
done
