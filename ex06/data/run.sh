#!/bin/bash
rm train_data/*
rm test_data/*

for i in $(seq 0.2 0.2 8.0)
do
    for j in 1
    do
	    #echo $i $j
        ./ising2d $i $j &
    done
done

wait 
