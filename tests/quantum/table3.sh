# Running the FULL-DDSIM case from 5 to 10 qubits
python3 q_search.py -to 1000 -t full_ddsim -m 5 -M 10 -r 5; # no randomness -> only 5 repetitions
git add results/[result]q_search_full_ddsim.csv; git commit -m "Grover Full DDSIM (with 5-10 qubits)"; git push;
python3 q_maxcut.py -to 1000 -t full_ddsim -m 5 -M 10 -r 50;
git add results/[result]q_maxcut_full_ddsim.csv; git commit -m "SAT Full DDSIM (with 5-10 qubits)"; git push;
python3 q_sat.py -to 1000 -t full_ddsim -m 5 -M 10 -r 50;
git add results/[result]q_sat_full_ddsim.csv; git commit -m "MAXCUT Full DDSIM (with 5-10 qubits)"; git push;
# Running the FULL-CLUE case from 5 to 10 qubits
python3 q_search.py -to 1000 -t full_clue -m 5 -M 10 -r 5; # no randomness -> only 5 repetitions
git add results/[result]q_search_full_clue.csv; git commit -m "Grover Full CLUE (with 5-10 qubits)"; git push;
python3 q_maxcut.py -to 1000 -t full_clue -m 5 -M 10 -r 50;
git add results/[result]q_maxcut_full_clue.csv; git commit -m "SAT Full CLUE (with 5-10 qubits)"; git push;
python3 q_sat.py -to 1000 -t full_clue -m 5 -M 10 -r 50;
git add results/[result]q_sat_full_clue.csv; git commit -m "MAXCUT Full CLUE (with 5-10 qubits)"; git push;

# Running the FULL-DDSIM case from 11 to 15 qubits
python3 q_search.py -to 1000 -t full_ddsim -m 11 -M 15 -r 5; # no randomness -> only 5 repetitions
git add results/[result]q_search_full_ddsim.csv; git commit -m "Grover Full DDSIM (with 11-15 qubits)"; git push;
python3 q_maxcut.py -to 1000 -t full_ddsim -m 11 -M 15 -r 50;
git add results/[result]q_maxcut_full_ddsim.csv; git commit -m "SAT Full DDSIM (with 11-15 qubits)"; git push;
python3 q_sat.py -to 1000 -t full_ddsim -m 11 -M 15 -r 50;
git add results/[result]q_sat_full_ddsim.csv; git commit -m "MAXCUT Full DDSIM (with 11-15 qubits)"; git push;
# Running the FULL-CLUE case from 11 to 15 qubits
python3 q_search.py -to 1000 -t full_clue -m 11 -M 15 -r 5; # no randomness -> only 5 repetitions
git add results/[result]q_search_full_clue.csv; git commit -m "Grover Full CLUE (with 11-15 qubits)"; git push;
python3 q_maxcut.py -to 1000 -t full_clue -m 11 -M 15 -r 50;
git add results/[result]q_maxcut_full_clue.csv; git commit -m "SAT Full CLUE (with 11-15 qubits)"; git push;
python3 q_sat.py -to 1000 -t full_clue -m 11 -M 15 -r 50;
git add results/[result]q_sat_full_clue.csv; git commit -m "MAXCUT Full CLUE (with 11-15 qubits)"; git push;

# Running the FULL-DDSIM case from 16 to 20 qubits
python3 q_search.py -to 1000 -t full_ddsim -m 16 -M 20 -r 5; # no randomness -> only 5 repetitions
git add results/[result]q_search_full_ddsim.csv; git commit -m "Grover Full DDSIM (with 16-20 qubits)"; git push;
python3 q_maxcut.py -to 1000 -t full_ddsim -m 16 -M 20 -r 50;
git add results/[result]q_maxcut_full_ddsim.csv; git commit -m "SAT Full DDSIM (with 16-20 qubits)"; git push;
python3 q_sat.py -to 1000 -t full_ddsim -m 16 -M 20 -r 50;
git add results/[result]q_sat_full_ddsim.csv; git commit -m "MAXCUT Full DDSIM (with 16-20 qubits)"; git push;
# Running the FULL-CLUE case from 16 to 20 qubits
python3 q_search.py -to 1000 -t full_clue -m 16 -M 20 -r 5; # no randomness -> only 5 repetitions
git add results/[result]q_search_full_clue.csv; git commit -m "Grover Full CLUE (with 16-20 qubits)"; git push;
python3 q_maxcut.py -to 1000 -t full_clue -m 16 -M 20 -r 50;
git add results/[result]q_maxcut_full_clue.csv; git commit -m "SAT Full CLUE (with 16-20 qubits)"; git push;
python3 q_sat.py -to 1000 -t full_clue -m 16 -M 20 -r 50;
git add results/[result]q_sat_full_clue.csv; git commit -m "MAXCUT Full CLUE (with 16-20 qubits)"; git push;
