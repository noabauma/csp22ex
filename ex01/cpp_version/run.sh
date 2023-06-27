#!/bin/bash
rm data.csv

for i in $(seq 0.2 0.4 8.0)
do
	#echo $i
    ./main2 $i > $i.txt &
done

wait 

cat *.txt > data.csv

rm *.txt

python3 plot.py data.csv
