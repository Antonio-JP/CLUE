#####################################################################
### TESTS FOR TABLE 1
#####################################################################
## Faster tests
python3 q_sat.py -M 8 -r 20;
python3 q_maxcut.py -M 8 -r 20;
python3 q_sat.py -M 10 -r 20 -t ddsim;
python3 q_maxcut.py -M 10 -r 20 -t ddsim;
## Test for Q_SAT
python3 q_sat.py -m 11 -M 20 -r 10 -t ddsim;
python3 q_sat.py -m 21 -M 25 -r 5 -t ddsim;
python3 q_sat.py -r 20;
## Test for MAXCUT
python3 q_maxcut.py -m 11 -M 20 -r 10 -t ddsim;
python3 q_maxcut.py -m 21 -M 25 -r 5 -t ddsim;
python3 q_maxcut.py -r 20;
## Test in case of finishing
python3 q_sat.py -M 20 -r 5;
python3 q_maxcut.py -M 20 -r 5;

#####################################################################
### TESTS FOR TABLE 2
#####################################################################
python3 q_sat.py -m 8 -M 9 -r 100 -t full_clue
python3 q_sat.py -m 8 -M 9 -r 100 -t full_ddsim
python3 q_maxcut.py -m 8 -M 9 -r 100 -t full_clue
python3 q_maxcut.py -m 8 -M 9 -r 100 -t full_ddsim