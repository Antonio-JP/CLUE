# Running al the case studies in the repository
## FIRST: We run only the 0 and entangled state for the test on the paper
python3 q_benchmark.py -to 1000 -n dj -m 5 -M 8 -r 50 -obs 0 -obs H;
git add [result]q_benchmark_clue.csv; git commit -m "DJ (0 and H) finished (up to 8 q-bits)"; git push;
python3 q_benchmark.py -to 1000 -n ghz -m 5 -M 8 -r 50 -obs 0 -obs H;
git add [result]q_benchmark_clue.csv; git commit -m "GHZ (0 and H) finished (up to 8 q-bits)"; git push;
python3 q_benchmark.py -to 1000 -n graphstate -m 5 -M 8 -r 50 -obs 0 -obs H;
git add [result]q_benchmark_clue.csv; git commit -m "graphstate (0 and H) finished (up to 8 q-bits)"; git push;
python3 q_benchmark.py -to 1000 -n pricingcall -m 2 -M 4 -r 50 -obs 0 -obs H;
git add [result]q_benchmark_clue.csv; git commit -m "pricingcall (0 and H) finished (up to 8 q-bits)"; git push;
python3 q_benchmark.py -to 1000 -n pricingput -m 2 -M 4 -r 50 -obs 0 -obs H;
git add [result]q_benchmark_clue.csv; git commit -m "pricingput (0 and H) finished (up to 8 q-bits)"; git push;
python3 q_benchmark.py -to 1000 -n qft -m 5 -M 8 -r 50 -obs 0 -obs H;
git add [result]q_benchmark_clue.csv; git commit -m "QFT (0 and H) finished (up to 8 q-bits)"; git push;
python3 q_benchmark.py -to 1000 -n qwalk -m 5 -M 8 -r 50 -obs 0 -obs H;
git add [result]q_benchmark_clue.csv; git commit -m "Q-Walk (0 and H) finished (up to 8 q-bits)"; git push;
python3 q_benchmark.py -to 1000 -n tsp -m 2 -M 3 -r 50 -obs 0 -obs H;
git add [result]q_benchmark_clue.csv; git commit -m "TSP (0 and H) finished (up to 8 q-bits)"; git push;
## SECOND: We run only the remaining states for the test on the paper
python3 q_benchmark.py -to 1000 -n dj -m 5 -M 8 -r 50;
git add [result]q_benchmark_clue.csv; git commit -m "DJ (all) finished (up to 8 q-bits)"; git push;
python3 q_benchmark.py -to 1000 -n ghz -m 5 -M 8 -r 50;
git add [result]q_benchmark_clue.csv; git commit -m "GHZ (all) finished (up to 8 q-bits)"; git push;
python3 q_benchmark.py -to 1000 -n gstate -m 5 -M 8 -r 50;
git add [result]q_benchmark_clue.csv; git commit -m "graphstate (all) finished (up to 8 q-bits)"; git push;
python3 q_benchmark.py -to 1000 -n pricingcall -m 2 -M 4 -r 50;
git add [result]q_benchmark_clue.csv; git commit -m "pricingcall (all) finished (up to 8 q-bits)"; git push;
python3 q_benchmark.py -to 1000 -n pricingput -m 2 -M 4 -r 50;
git add [result]q_benchmark_clue.csv; git commit -m "pricingput (all) finished (up to 8 q-bits)"; git push;
python3 q_benchmark.py -to 1000 -n qft -m 5 -M 8 -r 50;
git add [result]q_benchmark_clue.csv; git commit -m "QFT (all) finished (up to 8 q-bits)"; git push;
python3 q_benchmark.py -to 1000 -n qwalk -m 5 -M 8 -r 50;
git add [result]q_benchmark_clue.csv; git commit -m "Q-Walk (all) finished (up to 8 q-bits)"; git push;
python3 q_benchmark.py -to 1000 -n tsp -m 2 -M 3 -r 50;
git add [result]q_benchmark_clue.csv; git commit -m "TSP (all) finished (up to 8 q-bits)"; git push;
## THIRD: We run the remaining case studies for the states 0 and entangled
python3 q_benchmark.py -to 1000 -n portfolioqaoa -m 5 -M 8 -r 50 -obs 0 -obs H;
git add [result]q_benchmark_clue.csv; git commit -m "Portfolio (0 and H) finished (up to 8 q-bits)"; git push;
python3 q_benchmark.py -to 1000 -n qnn -m 5 -M 8 -r 50 -obs 0 -obs H;
git add [result]q_benchmark_clue.csv; git commit -m "Neural Network (0 and H) finished (up to 8 q-bits)"; git push;
python3 q_benchmark.py -to 1000 -n qpeexact -m 5 -M 8 -r 50 -obs 0 -obs H;
git add [result]q_benchmark_clue.csv; git commit -m "QPE (0 and H) finished (up to 8 q-bits)"; git push;
python3 q_benchmark.py -to 1000 -n qpeinexact -m 5 -M 8 -r 50 -obs 0 -obs H;
git add [result]q_benchmark_clue.csv; git commit -m "QPEI (0 and H) finished (up to 8 q-bits)"; git push;
python3 q_benchmark.py -to 1000 -n qwalk -m 5 -M 8 -r 50 -obs 0 -obs H;
git add [result]q_benchmark_clue.csv; git commit -m "Q-Walk (0 and H) finished (up to 8 q-bits)"; git push;
python3 q_benchmark.py -to 1000 -n vqe -m 5 -M 8 -r 50 -obs 0 -obs H;
git add [result]q_benchmark_clue.csv; git commit -m "VQE (0 and H) finished (up to 8 q-bits)"; git push;
python3 q_benchmark.py -to 1000 -n wstate -m 5 -M 8 -r 50 -obs 0 -obs H;
git add [result]q_benchmark_clue.csv; git commit -m "W-State (0 and H) finished (up to 8 q-bits)"; git push;
## FOURTH: We run the remaining case studies for the remaining states
python3 q_benchmark.py -to 1000 -n portfolioqaoa -m 5 -M 8 -r 50;
git add [result]q_benchmark_clue.csv; git commit -m "Portfolio (all) finished (up to 8 q-bits)"; git push;
python3 q_benchmark.py -to 1000 -n qnn -m 5 -M 8 -r 50;
git add [result]q_benchmark_clue.csv; git commit -m "Neural Network (all) finished (up to 8 q-bits)"; git push;
python3 q_benchmark.py -to 1000 -n qpeexact -m 5 -M 8 -r 50;
git add [result]q_benchmark_clue.csv; git commit -m "QPE (all) finished (up to 8 q-bits)"; git push;
python3 q_benchmark.py -to 1000 -n qpeinexact -m 5 -M 8 -r 50;
git add [result]q_benchmark_clue.csv; git commit -m "QPEI (all) finished (up to 8 q-bits)"; git push;
python3 q_benchmark.py -to 1000 -n qwalk -m 5 -M 8 -r 50;
git add [result]q_benchmark_clue.csv; git commit -m "Q-Walk (all) finished (up to 8 q-bits)"; git push;
python3 q_benchmark.py -to 1000 -n vqe -m 5 -M 8 -r 50;
git add [result]q_benchmark_clue.csv; git commit -m "VQE (all) finished (up to 8 q-bits)"; git push;
python3 q_benchmark.py -to 1000 -n wstate -m 5 -M 8 -r 50;
git add [result]q_benchmark_clue.csv; git commit -m "W-State (all) finished (up to 8 q-bits)"; git push;


