REPS=50
m=5
M=15

#####################################################################################################################
###
### TESTS FOR DATA ON TABLE 1
###
#####################################################################################################################
for s in $(seq $m $M); do # we change for each size
    # Executing tests for DDSIM column
    ## Grover
    for n in $(seq 5); do # this test is not random -> only 5 repetitions
        python3 q_search.py -to 500 -t full_ddsim -m $s -M $s -r 1; 
    done
    ## SAT
    for n in $(seq $REPS); do
        python3 q_sat.py -to 500 -t full_ddsim -m $s -M $s -r 1;
    done
    ## MaxCut
    for n in $(seq $REPS); do
        python3 q_maxcut.py -to 500 -t full_ddsim -m $s -M $s -r 1;
    done

    # Executing tests for Quokka# column
    ## Grover
    for n in $(seq 5); do # this test is not random -> only 5 repetitions
        python3 q_search.py -to 500 -t full_quokka# -m $s -M $s -r 1; 
    done
    ## SAT
    for n in $(seq $REPS); do
        python3 q_sat.py -to 500 -t full_quokka# -m $s -M $s -r 1;
    done
    ## MaxCut
    for n in $(seq $REPS); do
        python3 q_maxcut.py -to 500 -t full_quokka# -m $s -M $s -r 1;
    done

    # Executing tests for CLUE column
    ## Grover
    for n in $(seq 5); do # this test is not random -> only 5 repetitions
        python3 q_search.py -to 500 -t full_clue -m $s -M $s -r 1; 
    done
    ## SAT
    for n in $(seq $REPS); do
        python3 q_sat.py -to 500 -t full_direct -m $s -M $s -r 1;
    done
    ## MaxCut
    for n in $(seq $REPS); do
        python3 q_maxcut.py -to 500 -t full_direct -m $s -M $s -r 1;
    done

    # Executing tests for d column
    ## SAT
    for n in $(seq $REPS); do
        python3 q_sat.py -to 500 -t direct -m $s -M $s -r 1;
    done
    ## MaxCut
    for n in $(seq $REPS); do
        python3 q_maxcut.py -to 500 -t direct -m $s -M $s -r 1;
    done
done
