#!/bin/bash
rm outputs/*.csv

for ((j = 2; j <= 32; j *= 2)); do
    for i in 3.51 4.51 5.51
    do
        #echo "$j$i"
        ./ising3d_wolff_task3 $i $j &   #remove parallelization if less than 15 cores (to measure time)
        ./ising3d_mr2t2_task3 $i $j &   #I used a R9 5950x so no problem in my case
    done
    #wait 
done

wait 

python3 plot_task3.py
