#####################################################################################################################
###
### TESTS FOR DATA ON TABLE 2
###
#####################################################################################################################
# Tests for the column "d/N wrt S_{\ket{0}}"
python3 q_benchmark.py -to 500 -n ae -m 3 -M 7 -r 5;
python3 q_benchmark.py -to 500 -n dj -m 3 -M 7 -r 5;
python3 q_benchmark.py -to 500 -n ghz -m 3 -M 7 -r 5;
python3 q_benchmark.py -to 500 -n graphstate -m 3 -M 7 -r 5;
python3 q_benchmark.py -to 500 -n hhl -m 2 -M 4 -r 5;
python3 q_benchmark.py -to 500 -n pricingcall -m 2 -M 4 -r 5;
python3 q_benchmark.py -to 500 -n pricingput -m 2 -M 4 -r 5;
python3 q_benchmark.py -to 500 -n qft -m 3 -M 7 -r 5;
python3 q_benchmark.py -to 500 -n qwalk -m 3 -M 7 -r 5;
python3 q_benchmark.py -to 500 -n tsp -m 2 -M 3 -r 5;
python3 q_benchmark.py -to 500 -n portfolioqaoa -m 3 -M 7 -r 5;
python3 q_benchmark.py -to 500 -n portfoliovqe -m 3 -M 7 -r 5;
python3 q_benchmark.py -to 500 -n qnn -m 3 -M 7 -r 5;
python3 q_benchmark.py -to 500 -n qpeexact -m 3 -M 7 -r 5;
python3 q_benchmark.py -to 500 -n qpeinexact -m 3 -M 7 -r 5;
python3 q_benchmark.py -to 500 -n vqe -m 3 -M 7 -r 5;
python3 q_benchmark.py -to 500 -n wstate -m 3 -M 7 -r 5;

## Test for the column "DDSIM time"
python3 q_benchmark.py -to 500 -t full_ddsim -obs 0 -obs H -n ae -m 3 -M 7 -r 5;
python3 q_benchmark.py -to 500 -t full_ddsim -obs 0 -obs H -n dj -m 3 -M 7 -r 5;
python3 q_benchmark.py -to 500 -t full_ddsim -obs 0 -obs H -n ghz -m 3 -M 7 -r 5;
python3 q_benchmark.py -to 500 -t full_ddsim -obs 0 -obs H -n graphstate -m 3 -M 7 -r 5;
python3 q_benchmark.py -to 500 -t full_ddsim -obs 0 -obs H -n hhl -m 2 -M 4 -r 5;
python3 q_benchmark.py -to 500 -t full_ddsim -obs 0 -obs H -n pricingcall -m 2 -M 4 -r 5;
python3 q_benchmark.py -to 500 -t full_ddsim -obs 0 -obs H -n pricingput -m 2 -M 4 -r 5;
python3 q_benchmark.py -to 500 -t full_ddsim -obs 0 -obs H -n qft -m 3 -M 7 -r 5;
python3 q_benchmark.py -to 500 -t full_ddsim -obs 0 -obs H -n qwalk -m 3 -M 7 -r 5;
python3 q_benchmark.py -to 500 -t full_ddsim -obs 0 -obs H -n tsp -m 2 -M 3 -r 5;
python3 q_benchmark.py -to 500 -t full_ddsim -obs 0 -obs H -n portfolioqaoa -m 3 -M 7 -r 5;
python3 q_benchmark.py -to 500 -t full_ddsim -obs 0 -obs H -n portfoliovqe -m 3 -M 7 -r 5;
python3 q_benchmark.py -to 500 -t full_ddsim -obs 0 -obs H -n qnn -m 3 -M 7 -r 5;
python3 q_benchmark.py -to 500 -t full_ddsim -obs 0 -obs H -n qpeexact -m 3 -M 7 -r 5;
python3 q_benchmark.py -to 500 -t full_ddsim -obs 0 -obs H -n qpeinexact -m 3 -M 7 -r 5;
python3 q_benchmark.py -to 500 -t full_ddsim -obs 0 -obs H -n vqe -m 3 -M 7 -r 5;
python3 q_benchmark.py -to 500 -t full_ddsim -obs 0 -obs H -n wstate -m 3 -M 7 -r 5;
