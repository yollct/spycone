#!/bin/bash

for file in $(ls | grep "txt");
do
    sed -i $'s/\t/  /g' $file

done