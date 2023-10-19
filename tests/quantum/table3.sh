# Running the FULL-DDSIM case from 5 to 10 qubits
python3 q_search.py -to 500 -t full_ddsim -m 5 -M 10 -r 5; # no randomness -> only 5 repetitions
python3 q_maxcut.py -to 500 -t full_ddsim -m 5 -M 10 -r 50;
python3 q_sat.py -to 500 -t full_ddsim -m 5 -M 10 -r 50;
# Running the FULL-CLUE case from 5 to 10 qubits
python3 q_search.py -to 500 -t full_clue -m 5 -M 10 -r 5; # no randomness -> only 5 repetitions
python3 q_maxcut.py -to 500 -t full_direct -m 5 -M 10 -r 50;
python3 q_sat.py -to 500 -t full_direct -m 5 -M 10 -r 50;

# Running the FULL-DDSIM case from 11 to 15 qubits
python3 q_search.py -to 500 -t full_ddsim -m 11 -M 15 -r 5; # no randomness -> only 5 repetitions
python3 q_maxcut.py -to 500 -t full_ddsim -m 11 -M 15 -r 50;
python3 q_sat.py -to 500 -t full_ddsim -m 11 -M 15 -r 50;
# Running the FULL-CLUE case from 11 to 15 qubits
python3 q_search.py -to 500 -t full_clue -m 11 -M 15 -r 5; # no randomness -> only 5 repetitions
python3 q_maxcut.py -to 500 -t full_direct -m 11 -M 15 -r 50;
python3 q_sat.py -to 500 -t full_direct -m 11 -M 15 -r 50;
