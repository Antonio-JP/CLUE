# Executing examples Up to order 15
python3 q_sat.py -to 1000 -t direct -m 5 -M 15 -r 50;
git add [result]q_sat_direct.csv; git commit -m "Table 1 (SAT) up to order 15"; git push;
python3 q_maxcut.py -to 1000 -t direct -m 5 -M 15 -r 50;
git add [result]q_maxcut_direct.csv; git commit -m "Table 1 (MAXCUT) up to order 15"; git push;
python3 q_order.py -to 1000 -m 5 -M 12 -r 50;
git add [result]q_order_clue.csv; git commit -m "Table 1 (ORDER) up to order 12"; git push;
python3 q_search.py -to 1000 -m 5 -M 10 -r 50;
git add [result]q_search_clue.csv; git commit -m "Table 1 (GROVER) up to order 10"; git push;
# Executing examples Up to order 20
python3 q_sat.py -to 1000 -t direct -m 16 -M 20 -r 50;
git add [result]q_sat_direct.csv; git commit -m "Table 1 (SAT) up to order 20"; git push;
python3 q_maxcut.py -to 1000 -t direct -m 16 -M 20 -r 50;
git add [result]q_maxcut_direct.csv; git commit -m "Table 1 (MAXCUT) up to order 20"; git push;
python3 q_order.py -to 1000 -m 13 -M 15 -r 50;
git add [result]q_order_clue.csv; git commit -m "Table 1 (ORDER) up to order 20"; git push;
python3 q_search.py -to 1000 -m 11 -M 15 -r 50;
git add [result]q_search_clue.csv; git commit -m "Table 1 (GROVER) up to order 20"; git push;