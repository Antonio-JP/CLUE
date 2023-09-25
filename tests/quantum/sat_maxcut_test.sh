#####################################################################
### TESTS FOR TABLE 1
#####################################################################
git pull
## Faster tests
python3 q_sat.py -M 8 -r 20;
git add [result]q_sat_clue.csv;
git commit -m "Updated results of Q-SAT (clue)";
git push;
python3 q_maxcut.py -M 8 -r 20;
git add [result]q_maxcut_clue.csv;
git commit -m "Updated results of Q-MAXCUT (clue)";
git push;
python3 q_sat.py -M 10 -r 20 -t ddsim;
git add [result]q_sat_ddsim.csv;
git commit -m "Updated results of Q-SAT (ddsim)";
git push;
python3 q_maxcut.py -M 10 -r 20 -t ddsim;
git add [result]q_maxcut_ddsim.csv;
git commit -m "Updated results of Q-MAXCUT (ddsim)";
git push;
## Test for Q_SAT
python3 q_sat.py -m 11 -M 20 -r 10 -t ddsim;
git add [result]q_sat_ddsim.csv;
git commit -m "Updated results of Q-SAT (ddsim) from 11 to 20";
git push;
python3 q_sat.py -m 21 -M 25 -r 5 -t ddsim;
git add [result]q_sat_ddsim.csv;
git commit -m "Updated results of Q-SAT (ddsim) from 21 to 25";
git push;
python3 q_sat.py -r 20 -to 600;
git add [result]q_sat_clue.csv;
git commit -m "Updated results of Q-SAT (clue)";
git push;
## Test for MAXCUT
python3 q_maxcut.py -m 11 -M 20 -r 10 -t ddsim;
git add [result]q_maxcut_ddsim.csv;
git commit -m "Updated results of Q-MAXCUT (ddsim) from 11 to 20";
git push;
python3 q_maxcut.py -m 21 -M 25 -r 5 -t ddsim;
git add [result]q_maxcut_ddsim.csv;
git commit -m "Updated results of Q-MAXCUT (ddsim) from 21 to 25";
git push;
python3 q_maxcut.py -r 20 -to 600;
git add [result]q_maxcut_clue.csv;
git commit -m "Updated results of Q-MAXCUT (clue)";
git push;

#####################################################################
### TESTS FOR TABLE 2
#####################################################################
python3 q_sat.py -m 8 -M 9 -r 100 -t full_clue -to 600;
git add [result]q_sat_full_clue.csv;
git commit -m "Updated results of Q-SAT (full-clue)";
git push;
python3 q_sat.py -m 8 -M 9 -r 100 -t full_ddsim -to 600;
git add [result]q_sat_full_ddsim.csv;
git commit -m "Updated results of Q-SAT (full-ddsim)";
git push;
python3 q_maxcut.py -m 8 -M 9 -r 100 -t full_clue -to 600;
git add [result]q_maxcut_full_clue.csv;
git commit -m "Updated results of Q-MAXCUT (full-clue)";
git push;
python3 q_maxcut.py -m 8 -M 9 -r 100 -t full_ddsim -to 600;
git add [result]q_maxcut_full_ddsim.csv;
git commit -m "Updated results of Q-MAXCUT (full-ddsim)";
git push;


## Test in case of finishing
python3 q_sat.py -M 20 -r 5 -to 600;
python3 q_maxcut.py -M 20 -r 5 -to 600;

