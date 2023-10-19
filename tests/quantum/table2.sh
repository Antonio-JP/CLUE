# Running al the case studies in the repository
## FIRST: We run only the 0 and entangled state for the test on the paper
python3 q_benchmark.py -to 500 -n dj -m 3 -M 4 -r 5 -obs 0 -obs H;
python3 q_benchmark.py -to 500 -n ghz -m 3 -M 4 -r 5 -obs 0 -obs H;
python3 q_benchmark.py -to 500 -n graphstate -m 3 -M 4 -r 5 -obs 0 -obs H;
python3 q_benchmark.py -to 500 -n pricingcall -m 2 -M 4 -r 5 -obs 0 -obs H;
python3 q_benchmark.py -to 500 -n pricingput -m 2 -M 4 -r 5 -obs 0 -obs H;
python3 q_benchmark.py -to 500 -n qft -m 3 -M 4 -r 5 -obs 0 -obs H;
python3 q_benchmark.py -to 500 -n qwalk -m 3 -M 4 -r 5 -obs 0 -obs H;
python3 q_benchmark.py -to 500 -n tsp -m 2 -M 3 -r 5 -obs 0 -obs H;
## SECOND: We run only the remaining states for the test on the paper
python3 q_benchmark.py -to 500 -n dj -m 3 -M 4 -r 5;
python3 q_benchmark.py -to 500 -n ghz -m 3 -M 4 -r 5;
python3 q_benchmark.py -to 500 -n graphstate -m 3 -M 4 -r 5;
python3 q_benchmark.py -to 500 -n pricingcall -m 2 -M 4 -r 5;
python3 q_benchmark.py -to 500 -n pricingput -m 2 -M 4 -r 5;
python3 q_benchmark.py -to 500 -n qft -m 3 -M 4 -r 5;
python3 q_benchmark.py -to 500 -n qwalk -m 3 -M 4 -r 5;
python3 q_benchmark.py -to 500 -n tsp -m 2 -M 3 -r 5;
## THIRD: We run the DDSIM case for the table cases and input 0 and H
python3 q_benchmark.py -to 500 -t full_ddsim -obs 0 -obs H -n dj -m 3 -M 8 -r 5;
python3 q_benchmark.py -to 500 -t full_ddsim -obs 0 -obs H -n ghz -m 3 -M 8 -r 5;
python3 q_benchmark.py -to 500 -t full_ddsim -obs 0 -obs H -n graphstate -m 3 -M 8 -r 5;
python3 q_benchmark.py -to 500 -t full_ddsim -obs 0 -obs H -n pricingcall -m 2 -M 4 -r 5;
python3 q_benchmark.py -to 500 -t full_ddsim -obs 0 -obs H -n pricingput -m 2 -M 4 -r 5;
python3 q_benchmark.py -to 500 -t full_ddsim -obs 0 -obs H -n qft -m 3 -M 8 -r 5;
python3 q_benchmark.py -to 500 -t full_ddsim -obs 0 -obs H -n qwalk -m 3 -M 8 -r 5;
python3 q_benchmark.py -to 500 -t full_ddsim -obs 0 -obs H -n tsp -m 2 -M 3 -r 5;
## FOURTH: We run the remaining case studies for the states 0 and entangled
python3 q_benchmark.py -to 500 -n portfolioqaoa -m 5 -M 8 -r 5 -obs 0 -obs H;
python3 q_benchmark.py -to 500 -n qnn -m 5 -M 8 -r 5 -obs 0 -obs H;
python3 q_benchmark.py -to 500 -n qpeexact -m 5 -M 8 -r 5 -obs 0 -obs H;
python3 q_benchmark.py -to 500 -n qpeinexact -m 5 -M 8 -r 5 -obs 0 -obs H;
python3 q_benchmark.py -to 500 -n qwalk -m 5 -M 8 -r 5 -obs 0 -obs H;
python3 q_benchmark.py -to 500 -n vqe -m 5 -M 8 -r 5 -obs 0 -obs H;
python3 q_benchmark.py -to 500 -n wstate -m 5 -M 8 -r 5 -obs 0 -obs H;
# ## FIFTH: We run the remaining case studies for the remaining states
python3 q_benchmark.py -to 500 -n portfolioqaoa -m 5 -M 8 -r 5;
python3 q_benchmark.py -to 500 -n qnn -m 5 -M 8 -r 5;
python3 q_benchmark.py -to 500 -n qpeexact -m 5 -M 8 -r 5;
python3 q_benchmark.py -to 500 -n qpeinexact -m 5 -M 8 -r 5;
python3 q_benchmark.py -to 500 -n qwalk -m 5 -M 8 -r 5;
python3 q_benchmark.py -to 500 -n vqe -m 5 -M 8 -r 5;
python3 q_benchmark.py -to 500 -n wstate -m 5 -M 8 -r 5;


