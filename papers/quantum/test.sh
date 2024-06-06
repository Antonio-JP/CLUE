## Tests with several types for SAT
python3 q_sat.py -m 5 -M 5 -r 1 -t clue;
python3 q_sat.py -m 5 -M 5 -r 1 -t ddsim;
python3 q_sat.py -m 5 -M 5 -r 1 -t direct;
python3 q_sat.py -m 5 -M 5 -r 1 -t full_clue;
python3 q_sat.py -m 5 -M 5 -r 1 -t full_ddsim;
python3 q_sat.py -m 5 -M 5 -r 1 -t full_direct;
## Tests with several types for MAXCUT
python3 q_maxcut.py -m 5 -M 5 -r 1 -t clue;
python3 q_maxcut.py -m 5 -M 5 -r 1 -t ddsim;
python3 q_maxcut.py -m 5 -M 5 -r 1 -t direct;
python3 q_maxcut.py -m 5 -M 5 -r 1 -t full_clue;
python3 q_maxcut.py -m 5 -M 5 -r 1 -t full_ddsim;
python3 q_maxcut.py -m 5 -M 5 -r 1 -t full_direct;
## Tests with several types for SEARCH
python3 q_search.py -m 5 -M 5 -r 1 -t clue;
python3 q_search.py -m 5 -M 5 -r 1 -t ddsim;
# ---- no direct case because it is not implemented for Grover
python3 q_search.py -m 5 -M 5 -r 1 -t full_clue;
python3 q_search.py -m 5 -M 5 -r 1 -t full_ddsim;
# ---- no direct case because it is not implemented for Grover
## Tests with several types for ORDER
python3 q_order.py -m 5 -M 5 -r 1 -t clue;
python3 q_order.py -m 5 -M 5 -r 1 -t full_clue;
## Tests for ech benchmark family in the smallest case with observable 0 (so it is fast)
python3 q_benchmark.py -t clue -n ae -m 3 -M 3 -r 1 -obs 0;
python3 q_benchmark.py -t clue -n dj -m 3 -M 3 -r 1 -obs 0;
python3 q_benchmark.py -t clue -n ghz -m 3 -M 3 -r 1 -obs 0;
python3 q_benchmark.py -t clue -n graphstate -m 3 -M 3 -r 1 -obs 0;
python3 q_benchmark.py -t clue -n hhl -m 2 -M 2 -r 1 -obs 0;
python3 q_benchmark.py -t clue -n pricingput -m 2 -M 2 -r 1 -obs 0;
python3 q_benchmark.py -t clue -n pricingcall -m 2 -M 2 -r 1 -obs 0;
python3 q_benchmark.py -t clue -n portfolioqaoa -m 3 -M 3 -r 1 -obs 0;
python3 q_benchmark.py -t clue -n portfoliovqe -m 3 -M 3 -r 1 -obs 0;
python3 q_benchmark.py -t clue -n qft -m 3 -M 3 -r 1 -obs 0;
python3 q_benchmark.py -t clue -n qpeexact -m 3 -M 3 -r 1 -obs 0;
python3 q_benchmark.py -t clue -n qpeinexact -m 3 -M 3 -r 1 -obs 0;
python3 q_benchmark.py -t clue -n qwalk -m 3 -M 3 -r 1 -obs 0;
python3 q_benchmark.py -t clue -n tsp -m 2 -M 2 -r 1 -obs 0;
python3 q_benchmark.py -t clue -n qnn -m 3 -M 3 -r 1 -obs 0;
python3 q_benchmark.py -t clue -n vqe -m 3 -M 3 -r 1 -obs 0;
python3 q_benchmark.py -t clue -n wstate -m 3 -M 3 -r 1 -obs 0;

python3 q_benchmark.py -t full_ddsim -n ae -m 3 -M 3 -r 1 -obs 0;
python3 q_benchmark.py -t full_ddsim -n dj -m 3 -M 3 -r 1 -obs 0;
python3 q_benchmark.py -t full_ddsim -n ghz -m 3 -M 3 -r 1 -obs 0;
python3 q_benchmark.py -t full_ddsim -n graphstate -m 3 -M 3 -r 1 -obs 0;
python3 q_benchmark.py -t full_ddsim -n hhl -m 2 -M 2 -r 1 -obs 0;
python3 q_benchmark.py -t full_ddsim -n pricingput -m 2 -M 2 -r 1 -obs 0;
python3 q_benchmark.py -t full_ddsim -n pricingcall -m 2 -M 2 -r 1 -obs 0;
python3 q_benchmark.py -t full_ddsim -n portfolioqaoa -m 3 -M 3 -r 1 -obs 0;
python3 q_benchmark.py -t full_ddsim -n portfoliovqe -m 3 -M 3 -r 1 -obs 0;
python3 q_benchmark.py -t full_ddsim -n qft -m 3 -M 3 -r 1 -obs 0;
python3 q_benchmark.py -t full_ddsim -n qpeexact -m 3 -M 3 -r 1 -obs 0;
python3 q_benchmark.py -t full_ddsim -n qpeinexact -m 3 -M 3 -r 1 -obs 0;
python3 q_benchmark.py -t full_ddsim -n qwalk -m 3 -M 3 -r 1 -obs 0;
python3 q_benchmark.py -t full_ddsim -n tsp -m 2 -M 2 -r 1 -obs 0;
python3 q_benchmark.py -t full_ddsim -n qnn -m 3 -M 3 -r 1 -obs 0;
python3 q_benchmark.py -t full_ddsim -n vqe -m 3 -M 3 -r 1 -obs 0;
python3 q_benchmark.py -t full_ddsim -n wstate -m 3 -M 3 -r 1 -obs 0;