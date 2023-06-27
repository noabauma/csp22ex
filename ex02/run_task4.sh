#!/bin/bash
rm task4_data.csv

for i in $(seq 0.1 0.4 10.0)
do
    for j in {3..10..1}
    do 
        #echo $i $j
        ./ising3d_task4 $i $j > tmp_data/$i$j.txt &
    done
done

wait 

cat tmp_data/*.txt > task4_data.csv

rm tmp_data/*.txt

python3 plot_task4.py task4_data.csv
