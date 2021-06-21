#!/bin/bash

make -j 

N=(4 8 16 32)
selection=2
k1=8
k2=1
k3=1
file=./output.txt

echo "Running code: " > "$file"
 

for N in "${N[@]}";
do
    printf "N=$N selection=$selection k1=$k1 k2=$k2 k3=$k3 \n" >> "$file"
    ./../bin/dolphin $N $selection $k1 $k2 $k3 >> "$file" 
done
