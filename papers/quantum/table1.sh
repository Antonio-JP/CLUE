#####################################################################################################################
###
### TESTS FOR DATA ON TABLE 1
###
#####################################################################################################################
# Executing tests for DDSIM column
## Grover
python3 q_search.py -to 500 -t full_ddsim -m 5 -M 15 -r 5; # this test is not random -> only 5 repetitions
## SAT
python3 q_sat.py -to 500 -t full_ddsim -m 5 -M 15 -r 50;
## MaxCut
python3 q_maxcut.py -to 500 -t full_ddsim -m 5 -M 15 -r 50;

# Executing tests for CLUE column
## Grover
python3 q_search.py -to 500 -t full_clue -m 5 -M 15 -r 5; # this test is not random -> only 5 repetitions
## SAT
python3 q_sat.py -to 500 -t full_direct -m 5 -M 15 -r 50;
## MaxCut
python3 q_maxcut.py -to 500 -t full_direct -m 5 -M 15 -r 50;

# Executing tests for d column
## SAT
python3 q_sat.py -to 500 -t direct -m 5 -M 15 -r 50;
## MaxCut
python3 q_maxcut.py -to 500 -t direct -m 5 -M 15 -r 50;
