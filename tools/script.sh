#!/bin/bash
for fileName in `ls run75*.root`
do
    echo $fileName
    root -l -q "sort.cpp+O(\"$fileName\")"
done
