#####################################################################################################################
###
### TESTS FOR DATA ON TABLE 1
###
#####################################################################################################################
MY_SCRIPT=./CLUE_main_script

## This are the different types of experiments
EXPERIMENTS="ae dj ghz graphstate hhl portfolioqaoa portfoliovqe pricingcall pricingput \
             qft qnn qpeexact qpeinexact qwalk tsp vqe wstate"
OBSERVABLE="0"
## This are the different types of approaches
TESTS="DDSIM_ALONE DDSIM CLUE"
REPEATS=5
MIN=3
MAX=25
TIMEOUT=500s

cd ../../cpp/clue/build/apps


for size in $(seq $MIN 1 $MAX)
do
    for repeat in $(seq 1 1 $REPEATS)
    do
        for experiment in $EXPERIMENTS
        do
            for observable in $OBSERVABLE
            do
                for test in $TESTS
                do
                    timeout $TIMEOUT $MY_SCRIPT $experiment -m $size -M $size -r 1 -t $test -obs $observable
                done
            done
        done
    done
done