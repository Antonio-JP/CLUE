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
# Tests for the columns "d/N wrt S_{\ket{0}}" and "Avg. ..."
python3 q_benchmark.py -to 500 -n ae -m 3 -M 4 -r 5;
python3 q_benchmark.py -to 500 -n dj -m 3 -M 4 -r 5;
python3 q_benchmark.py -to 500 -n ghz -m 3 -M 4 -r 5;
python3 q_benchmark.py -to 500 -n graphstate -m 3 -M 4 -r 5;
# python3 q_benchmark.py -to 500 -n hhl -m 2 -M 5 -r 5;
# python3 q_benchmark.py -to 500 -n pricingcall -m 2 -M 4 -r 5;
# python3 q_benchmark.py -to 500 -n pricingput -m 2 -M 4 -r 5;
python3 q_benchmark.py -to 500 -n qft -m 3 -M 4 -r 5;
python3 q_benchmark.py -to 500 -n qwalk -m 3 -M 4 -r 5;
python3 q_benchmark.py -to 500 -n tsp -m 2 -M 2 -r 5;
python3 q_benchmark.py -to 500 -n portfolioqaoa -m 3 -M 4 -r 5;
python3 q_benchmark.py -to 500 -n portfoliovqe -m 3 -M 4 -r 5;
python3 q_benchmark.py -to 500 -n qnn -m 3 -M 4 -r 5;
python3 q_benchmark.py -to 500 -n qpeexact -m 3 -M 4 -r 5;
python3 q_benchmark.py -to 500 -n qpeinexact -m 3 -M 4 -r 5;
python3 q_benchmark.py -to 500 -n vqe -m 3 -M 4 -r 5;
python3 q_benchmark.py -to 500 -n wstate -m 3 -M 4 -r 5;

## Test for the column "DDSIM time"
python3 q_benchmark.py -to 500 -t full_ddsim -obs 0 -obs H -n ae -m 3 -M 4 -r 5;
python3 q_benchmark.py -to 500 -t full_ddsim -obs 0 -obs H -n dj -m 3 -M 4 -r 5;
python3 q_benchmark.py -to 500 -t full_ddsim -obs 0 -obs H -n ghz -m 3 -M 4 -r 5;
python3 q_benchmark.py -to 500 -t full_ddsim -obs 0 -obs H -n graphstate -m 3 -M 4 -r 5;
# python3 q_benchmark.py -to 500 -t full_ddsim -obs 0 -obs H -n hhl -m 2 -M 5 -r 5;
# python3 q_benchmark.py -to 500 -t full_ddsim -obs 0 -obs H -n pricingcall -m 2 -M 4 -r 5;
# python3 q_benchmark.py -to 500 -t full_ddsim -obs 0 -obs H -n pricingput -m 2 -M 4 -r 5;
python3 q_benchmark.py -to 500 -t full_ddsim -obs 0 -obs H -n qft -m 3 -M 4 -r 5;
python3 q_benchmark.py -to 500 -t full_ddsim -obs 0 -obs H -n qwalk -m 3 -M 4 -r 5;
python3 q_benchmark.py -to 500 -t full_ddsim -obs 0 -obs H -n tsp -m 2 -M 2 -r 5;
python3 q_benchmark.py -to 500 -t full_ddsim -obs 0 -obs H -n portfolioqaoa -m 3 -M 4 -r 5;
python3 q_benchmark.py -to 500 -t full_ddsim -obs 0 -obs H -n portfoliovqe -m 3 -M 4 -r 5;
python3 q_benchmark.py -to 500 -t full_ddsim -obs 0 -obs H -n qnn -m 3 -M 4 -r 5;
python3 q_benchmark.py -to 500 -t full_ddsim -obs 0 -obs H -n qpeexact -m 3 -M 4 -r 5;
python3 q_benchmark.py -to 500 -t full_ddsim -obs 0 -obs H -n qpeinexact -m 3 -M 4 -r 5;
python3 q_benchmark.py -to 500 -t full_ddsim -obs 0 -obs H -n vqe -m 3 -M 4 -r 5;
python3 q_benchmark.py -to 500 -t full_ddsim -obs 0 -obs H -n wstate -m 3 -M 4 -r 5;