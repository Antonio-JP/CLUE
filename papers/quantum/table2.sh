#####################################################################################################################
###
### TESTS FOR DATA ON TABLE 2
###
#####################################################################################################################
MQT_BENCH_NAMES=$(python -c "from mqt.bench.benchmarks import *; print(str(get_available_benchmark_names()).replace('\'','').replace('[','').replace(',','').replace(']',''))")
REPS=5
m=3
M=7

for s in ($seq $m $M); do # we change for each size
    # Tests for the column "d/N wrt S_{\ket{0}}"
    for bench in $MQT_BENCH_NAMES; do
        for n in $(seq $REPS); do
            python3 q_benchmark.py -to 500 -n $bench -m $s -M $s -r 1;
        done
    done

    ## Test for the column "DDSIM time"
    for bench in $MQT_BENCH_NAMES; do
        for n in $(seq $REPS); do
            python3 q_benchmark.py -to 500 -t full_ddsim -obs 0 -obs H -n $bench -m $s -M $s -r 1;
        done
    done

    ## Test for the column "Quokka# time"
    for bench in $MQT_BENCH_NAMES; do
        for n in $(seq $REPS); do
            python3 q_benchmark.py -to 500 -t full_quokka# -obs 0 -obs H -n $bench -m $s -M $s -r 1;
        done
    done
done