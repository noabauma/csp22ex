#!/bin/bash
rm data_task1.csv

for ((j = 2; j <= 16; j *= 2)); do
    for i in $(seq 0.1 0.1 5.0)
    do
        #echo "$j$i"
        ./ising3d_mr2t2 $i $j > tmp_data/$j$i.txt &
    done
    #wait 
done

wait 

cat tmp_data/*.txt > data_task1.csv

rm tmp_data/*.txt

python3 plot_task1_2.py data_task1.csv
