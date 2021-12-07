#!/bin/bash

runnb=$1
volume1=$2
volume2=$3

echo "Put Parameters: run_nbr; volume_from; volume_to"
echo "Now starting run the selector $runnb"	

volnb=$volume1

#echo "file: run$runnb""_""$volnb"_"ssgant1.root" 
#echo "$vol1 $vol2"
#exit

while test $volnb -le $volume2
do 
echo "file: run$runnb""_""$volnb"_"ssgant1.root" 
  if test -f "run$runnb""_""$volnb"_"ssgant1.root"
  then
   echo "file: run$runnb""_""$volnb"_"ssgant1.root FOUND"
   rootcommand=sort.cpp"($runnb,$volnb)"  
   echo "$rootcommand"
   root -l -b -q $rootcommand  
  fi  
 volnb=$(($volnb + 1))  
done


