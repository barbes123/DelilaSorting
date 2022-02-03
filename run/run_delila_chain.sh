#!/bin/bash

FIRSTrun=$1
#LASTrun=$3
volume1=$2
volume2=$3

#echo "Will run the selector for runs $FIRSTrun upto $LASTrun "
#echo "Will run the selector for volumes $volume1 upto $volume2 "		

echo "Put Parameters: run_nbr; volume_from; volume_to"

runnb=$FIRSTrun
LASTrun=$FIRSTrun

while test $runnb -le $LASTrun
do
 echo "Now starting run the selector $runnb"	
 volnb=$volume1
 while test $volnb -le $volume2
 do 
   if test -f "run$runnb_$volnb.root" 
   then
     volnb=$(($volnb + 1))
   else 
     break  
   fi  
 done

 rootcommand=Delila_selector_chain.C+"($runnb,$runnb,$volume1,$volume2)"  

 root -l -b -q $rootcommand
 
 runnb=$(($runnb + 1))
done 
