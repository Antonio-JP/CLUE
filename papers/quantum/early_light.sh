REPS=5
#####################################################################################################################
###
### TESTS FOR DATA ON TABLE 1
###
#####################################################################################################################
# Executing tests for DDSIM column
## Grover
for n in $(seq $REPS); do
    python3 q_search.py -to 500 -t full_ddsim -m 5 -M 6 -r 1; # this test is not random -> only 5 repetitions
done
## SAT
for n in $(seq $REPS); do
    python3 q_sat.py -to 500 -t full_ddsim -m 5 -M 6 -r 1;
done
## MaxCut
for n in $(seq $REPS); do
    python3 q_maxcut.py -to 500 -t full_ddsim -m 5 -M 6 -r 1;
done

# Executing tests for Quokka# column
## Grover
for n in $(seq $REPS); do
    python3 q_search.py -to 500 -t full_quokka# -m 5 -M 6 -r 1; # this test is not random -> only 5 repetitions
done
## SAT
for n in $(seq $REPS); do
    python3 q_sat.py -to 500 -t full_quokka# -m 5 -M 6 -r 1;
done
## MaxCut
for n in $(seq $REPS); do
    python3 q_maxcut.py -to 500 -t full_quokka# -m 5 -M 6 -r 1;
done

# Executing tests for CLUE column
## Grover
for n in $(seq $REPS); do
    python3 q_search.py -to 500 -t full_clue -m 5 -M 6 -r 1; # this test is not random -> only 5 repetitions
done
## SAT
for n in $(seq $REPS); do
    python3 q_sat.py -to 500 -t full_direct -m 5 -M 6 -r 1;
done
## MaxCut
for n in $(seq $REPS); do
    python3 q_maxcut.py -to 500 -t full_direct -m 5 -M 6 -r 1;
done

# Executing tests for d column
## SAT
for n in $(seq $REPS); do
    python3 q_sat.py -to 500 -t direct -m 5 -M 6 -r 1;
done
## MaxCut
for n in $(seq $REPS); do
    python3 q_maxcut.py -to 500 -t direct -m 5 -M 6 -r 1;
done

#####################################################################################################################
###
### TESTS FOR DATA ON TABLE 2
###
#####################################################################################################################
MQT_BENCH_NAMES=$(python -c "from mqt.bench.benchmarks import *; print(str(get_available_benchmark_names()).replace('\'','').replace('[','').replace(',','').replace(']',''))")

# Tests for the column "d/N wrt S_{\ket{0}}"
for bench in $MQT_BENCH_NAMES; do
    for n in $(seq $REPS); do
        python3 q_benchmark.py -to 500 -n $bench -m 3 -M 4 -r 1;
    done
done

## Test for the column "DDSIM time"
for bench in $MQT_BENCH_NAMES; do
    for n in $(seq $REPS); do
        python3 q_benchmark.py -to 500 -t full_ddsim -obs 0 -obs H -n $bench -m 3 -M 4 -r 1;
    done
done

## Test for the column "Quokka# time"
for bench in $MQT_BENCH_NAMES; do
    for n in $(seq $REPS); do
        python3 q_benchmark.py -to 500 -t full_quokka# -obs 0 -obs H -n $bench -m 3 -M 4 -r 1;
    done
done
