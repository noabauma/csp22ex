#!/bin/bash
rm data_task2.csv

for ((j = 2; j <= 16; j += 2)); do
    for i in 0.1 1.0 2.0 3.0 3.5 3.75 4.0 4.125 4.25 4.375 4.45 4.5 4.55 4.625 4.75 4.875 5.0 5.25 5.5 6.0 7.0 8.0 9.0 10.0
    do
        #echo "$j$i"
        ./ising3d_wolff $i $j > tmp_data/$j$i.txt &
        #./ising3d_mr2t2 $i $j > tmp_data/$j$i.txt &
    done
    #wait 
done

wait 

cat tmp_data/*.txt > data_task2.csv

rm tmp_data/*.txt

python3 plot_task2.py data_task2.csv
