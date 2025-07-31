#####################################################################################################################
###
### TESTS FOR DATA ON TABLE 1
###
#####################################################################################################################
# Executing tests for DDSIM column
## Grover
python3 q_search.py -to 500 -t full_ddsim -m 5 -M 6 -r 5; # this test is not random -> only 5 repetitions
## SAT
python3 q_sat.py -to 500 -t full_ddsim -m 5 -M 6 -r 5;
## MaxCut
python3 q_maxcut.py -to 500 -t full_ddsim -m 5 -M 6 -r 5;

# Executing tests for Quokka# column
## Grover
python3 q_search.py -to 500 -t full_quokka# -m 5 -M 6 -r 5; # this test is not random -> only 5 repetitions
## SAT
python3 q_sat.py -to 500 -t full_quokka# -m 5 -M 6 -r 5;
## MaxCut
python3 q_maxcut.py -to 500 -t full_quokka# -m 5 -M 6 -r 5;

# Executing tests for CLUE column
## Grover
python3 q_search.py -to 500 -t full_clue -m 5 -M 6 -r 5; # this test is not random -> only 5 repetitions
## SAT
python3 q_sat.py -to 500 -t full_direct -m 5 -M 6 -r 5;
## MaxCut
python3 q_maxcut.py -to 500 -t full_direct -m 5 -M 6 -r 5;

# Executing tests for d column
## SAT
python3 q_sat.py -to 500 -t direct -m 5 -M 6 -r 5;
## MaxCut
python3 q_maxcut.py -to 500 -t direct -m 5 -M 6 -r 5;

#####################################################################################################################
###
### TESTS FOR DATA ON TABLE 2
###
#####################################################################################################################
MQT_BENCH_NAMES=$(python -c "from mqt.bench.benchmarks import *; print(str(get_available_benchmark_names()).replace('\'','').replace('[','').replace(',','').replace(']',''))")

# Tests for the column "d/N wrt S_{\ket{0}}"
for bench in $MQT_BENCH_NAMES; do
    python3 q_benchmark.py -to 500 -n $bench -m 3 -M 4 -r 5;
done

## Test for the column "DDSIM time"
for bench in $MQT_BENCH_NAMES; do
    python3 q_benchmark.py -to 500 -t full_ddsim -obs 0 -obs H -n $bench -m 3 -M 4 -r 5;
done

## Test for the column "Quokka# time"
for bench in $MQT_BENCH_NAMES; do
    python3 q_benchmark.py -to 500 -t full_quokka# -obs 0 -obs H -n $bench -m 3 -M 4 -r 5;
done
