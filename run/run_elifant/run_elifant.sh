#!/bin/bash

FIRSTrun=$1
#LASTrun=$3
volume1=${2:-0}
volume2=${3:-$volume1}
numberofevents=${4:-0}

#echo "Will run the selector for runs $FIRSTrun upto $LASTrun "
#echo "Will run the selector for volumes $volume1 upto $volume2 "		


echo "Put Parameters: run_nbr; volume_from; volume_to; number of events"

lut_path="$HOME/onlineAnalysis/LookUpTables/2022_w10_paul/"
lut_link="$HOME/DelilaSorting/"


lut_file="LUT_DELILA_Co_Eu_1191_1195_1196.dat"
lut_conf="coinc_gates_test_new.dat"
lut_ta="TimeCalib_run_1195.dat"

echo "$lut_path$lut_file" "$lut_link""LUT_DELILA.dat"

unlink "$lut_link""LUT_DELILA.dat"
ln -s "$lut_path$lut_file" "$lut_link""LUT_DELILA.dat"

unlink "$lut_link""LUT_CONF.dat"
ln -s "$lut_path$lut_conf" "$lut_link""LUT_CONF.dat"

unlink "$lut_link""LUT_TA.dat"
ln -s "$lut_path$lut_ta" "$lut_link""LUT_TA.dat"

runnb=$FIRSTrun

volnb=$volume1
while test $volnb -le $volume2
do
     rootcommand=Delila_selector_elifant.C+"($runnb,$runnb,$volnb,$volnb,$numberofevents)"  
     root -l -b -q $rootcommand
     volnb=$(($volnb + 1))
done 
