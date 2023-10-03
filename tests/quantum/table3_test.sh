tsp -S 1; #number of parallel processes
tsp python3 q_maxcut.py -t full_ddsim -m 5 -M 9 -r 50;
tsp python3 q_sat.py -t full_ddsim -m 5 -M 9 -r 50;
tsp python3 q_maxcut.py -t full_clue -m 5 -M 9 -r 50;
tsp python3 q_sat.py -t full_clue -m 5 -M 9 -r 50;
tsp python3 q_search.py -t full_ddsim -m 5 -M 9 -r 50;
tsp python3 q_search.py -t full_clue -m 5 -M 9 -r 50;
tsp python3 q_maxcut.py -t ddsim -m 5 -M 9 -r 50;
tsp python3 q_sat.py -t ddsim -m 5 -M 9 -r 50;
tsp python3 q_search.py -t ddsim -m 5 -M 9 -r 50;

