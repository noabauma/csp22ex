#!/bin/bash
rm task3_data.csv

for ((j = 2; j <= 32; j *= 2)); do
    for ((L = j, N = L*L*L, k = 63; k >= 0; --k )); do
        i=$((-3*N*k/63 ))
        #echo "$j$i"
        ./ising3D_creutz $i $j > tmp_data/$j$i.txt &
    done
    wait 
done

wait 

cat tmp_data/*.txt > task3_data.csv

rm tmp_data/*.txt

python3 plot_task3.py task3_data.csv
