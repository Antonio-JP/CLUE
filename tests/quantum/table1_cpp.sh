#####################################################################################################################
###
### TESTS FOR DATA ON TABLE 1
###
#####################################################################################################################
MY_SCRIPT=./CLUE_main_script

## This are the different types of experiments
EXPERIMENTS="search sat maxcut"
## This are the different types of approaches
TESTS="DDSIM_ALONE DDSIM CLUE"
REPEATS=50
MIN=5
MAX=50
TIMEOUT=500s

cd ../../cpp/clue/build/apps


for size in $(seq $MIN 1 $MAX)
do
    for repeat in $(seq 1 1 $REPEATS)
    do
        for experiment in $EXPERIMENTS
        do
            for test in $TESTS
            do
                echo "timeout $TIMEOUT $MY_SCRIPT $experiment -m $size -M $size -r 1 -t $test"
            done
        done
    done
done
