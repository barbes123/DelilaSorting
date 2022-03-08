#!/bin/bash

FIRSTrun=$1
#LASTrun=$3
volume1=$2
volume2=$3
numberofevents=$4

#echo "Will run the selector for runs $FIRSTrun upto $LASTrun "
#echo "Will run the selector for volumes $volume1 upto $volume2 "		


echo "Put Parameters: run_nbr; volume_from; volume_to; number of events"

runnb=$FIRSTrun

volnb=$volume1
while test $volnb -le $volume2
do
     rootcommand=Delila_selector_elifant.C+"($runnb,$runnb,$volnb,$volnb,$numberofevents)"  
     root -l -b -q $rootcommand
     volnb=$(($volnb + 1))
done 
