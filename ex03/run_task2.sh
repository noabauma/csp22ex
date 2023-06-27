#!/bin/bash
rm task2_data.csv

for ((L = 16, N = L*L*L, k = 63; k >= 0; --k )); do
    i=$((-3*N*k/63 ))
    #echo "$i"
    ./ising3D_creutz $i > tmp_data/$i.txt &
done

wait 

cat tmp_data/*.txt > task2_data.csv

rm tmp_data/*.txt

python3 plot.py task2_data.csv
