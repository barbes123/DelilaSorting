#!/bin/bash


FIRSTrun=$1
volume1=$2
volume2=$3

echo "Will run the selector for runs $FIRSTrun upto $LASTrun "
echo "Will run the selector for volumes $volume1 upto $volume2 "		


space=" "

runnb=$FIRSTrun
volnb=$volume1
command="hadd selected_run_$runnb"."root "


if test -f "selected_run_$runnb"."root"
   then
   rm "selected_run_$runnb"."root"
fi

while test $volnb -le $volume2
 do 
 echo "selected_run_$runnb""_$volnb.root"
   if test -f "selected_run_$runnb""_$volnb.root"
   then
     command="$command""selected_run_$runnb""_$volnb.root""$space"
     echo "$command"
     volnb=$(($volnb + 1))
     else 
     break  
   fi  
done

echo "$command"
$command
