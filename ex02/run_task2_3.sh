#!/bin/bash
rm task2_3_data.csv

#for i in $(seq 0.1 0.4 10.0)
for i in 0.1 1.0 2.0 3.0 3.5 3.75 4.0 4.125 4.25 4.375 4.45 4.5 4.55 4.625 4.75 4.875 5.0 5.25 5.5 6.0 7.0 8.0 9.0 10.0
do
	#echo $i
    ./ising3d_task2_3 $i > tmp_data/$i.txt &
done

wait 

cat tmp_data/*.txt > task2_3_data.csv

rm tmp_data/*.txt

python3 plot_task2_3.py task2_3_data.csv
