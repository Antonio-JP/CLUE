## Tests with several types for SAT
python3 q_sat.py -m 5 -M 5 -r 1 -t clue;
python3 q_sat.py -m 5 -M 5 -r 1 -t ddsim;
python3 q_sat.py -m 5 -M 5 -r 1 -t direct;
python3 q_sat.py -m 5 -M 5 -r 1 -t full_clue;
python3 q_sat.py -m 5 -M 5 -r 1 -t full_ddsim;
python3 q_sat.py -m 5 -M 5 -r 1 -t full_quokka#;
python3 q_sat.py -m 5 -M 5 -r 1 -t full_direct;
## Tests with several types for MAXCUT
python3 q_maxcut.py -m 5 -M 5 -r 1 -t clue;
python3 q_maxcut.py -m 5 -M 5 -r 1 -t ddsim;
python3 q_maxcut.py -m 5 -M 5 -r 1 -t direct;
python3 q_maxcut.py -m 5 -M 5 -r 1 -t full_clue;
python3 q_maxcut.py -m 5 -M 5 -r 1 -t full_ddsim;
python3 q_maxcut.py -m 5 -M 5 -r 1 -t full_quokka#;
python3 q_maxcut.py -m 5 -M 5 -r 1 -t full_direct;
## Tests with several types for SEARCH
python3 q_search.py -m 5 -M 5 -r 1 -t clue;
python3 q_search.py -m 5 -M 5 -r 1 -t ddsim;
# ---- no direct case because it is not implemented for Grover
python3 q_search.py -m 5 -M 5 -r 1 -t full_clue;
python3 q_search.py -m 5 -M 5 -r 1 -t full_ddsim;
python3 q_search.py -m 5 -M 5 -r 1 -t full_quokka#;
# ---- no direct case because it is not implemented for Grover
## Tests with several types for ORDER
python3 q_order.py -m 5 -M 5 -r 1 -t clue;
python3 q_order.py -m 5 -M 5 -r 1 -t full_clue;
## Tests for ech benchmark family in the smallest case with observable 0 (so it is fast)

MQT_BENCH_NAMES=$(python -c "from mqt.bench.benchmarks import *; print(str(get_available_benchmark_names()).replace('\'','').replace('[','').replace(',','').replace(']',''))")
for bench in $MQT_BENCH_NAMES; do
    python3 q_benchmark.py -t clue -n $bench -m 3 -M 3 -r 1 -obs 0;
    python3 q_benchmark.py -t full_ddsim -n $bench -m 3 -M 3 -r 1 -obs 0;
    python3 q_benchmark.py -t full_quokka# -n $bench -m 3 -M 3 -r 1 -obs 0;
done