## Test for Q_SAT
python3 q_sat.py -M 10 -r 20 -t ddsim;
python3 q_sat.py -m 11 -M 20 -r 10 -t ddsim;
python3 q_sat.py -m 21 -M 25 -r 5 -t ddsim;
python3 q_sat.py -r 20;
## Test for Q_SAT
python3 q_maxcut.py -M 10 -r 20 -t ddsim;
python3 q_maxcut.py -m 11 -M 20 -r 10 -t ddsim;
python3 q_maxcut.py -m 21 -M 25 -r 5 -t ddsim;
python3 q_maxcut.py -r 20;
## Test in case of finishing
python3 q_sat.py -M 20 -r 5;
python3 q_maxcut.py -M 20 -r 5;

