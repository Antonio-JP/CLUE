#####################################################################################################################
###
### TESTS FOR DATA ON TABLE 1
###
#####################################################################################################################
MY_SCRIPT=./CLUE_main_script

cd ../../cpp/clue/build/apps

# Executing tests for DDSIM column
## Grover
$MY_SCRIPT search -m 5 -M 15 -repeats 50 -t DDSIM_ALONE
$MY_SCRIPT search -m 5 -M 15 -repeats 50 -t DDSIM
$MY_SCRIPT search -m 5 -M 15 -repeats 50 -t CLUE

## SAT
$MY_SCRIPT sat -m 5 -M 15 -repeats 50 -t DDSIM_ALONE
$MY_SCRIPT sat -m 5 -M 15 -repeats 50 -t DDSIM
$MY_SCRIPT sat -m 5 -M 15 -repeats 50 -t DIRECT
$MY_SCRIPT sat -m 5 -M 15 -repeats 50 -t CLUE

## MAXCUT
$MY_SCRIPT maxcut -m 5 -M 15 -repeats 50 -t DDSIM_ALONE
$MY_SCRIPT maxcut -m 5 -M 15 -repeats 50 -t DDSIM
$MY_SCRIPT maxcut -m 5 -M 15 -repeats 50 -t DIRECT
$MY_SCRIPT maxcut -m 5 -M 15 -repeats 50 -t CLUE

