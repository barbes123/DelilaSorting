#!/bin/bash

#echo "MASS first last"

# DIST=$1
FIRSTdomain=$1
LASTdomain=$2
activity152Eu=$3
#activityCs137=$3

CUR_DIR=$(pwd)
cd $CUR_DIR

echo "CUR_DIR $CUR_DIR"


domnb=$FIRSTdomain
echo "mDelila_raw_py_$domnb.spe" 
 if test -f "mDelila_raw_py_$domnb.spe" 
   then
   echo "yessssssssssss"
   else
   echo "nooooooooooooooooo"
 fi

echo "Will run the calib for domains $FIRSTdomain upto $LASTdomain "

##################################   
 if test -f "eliade.calib" 
   then 
   	rm "eliade.calib"   
   fi
##################################   
    if test -f "res.temp" 
   then 
   	rm "res.temp"   
   fi
##################################      
    if test -f "152Eu.temp" 
   then 
   	rm "152Eu.temp"   
   fi
##################################      
    if test -f "temp.calib" 
   then 
   	rm "temp.calib"   
   fi
##################################
    if test -f "resolution.txt" 
   then 
   	rm "resolution.txt"   
   fi
##################################
    if test -f "fulldata.calib" 
   then 
   	rm "fulldata.calib"   
   fi
##################################      
   
   
domnb=$FIRSTdomain
while test $domnb -le $LASTdomain
do
 echo "Now starting calib for domain $domnb"	
 
 if [[ "$domnb" == "35" ]];
 then
  lim1=50
  lim2=800
  fwhm=7
  ampl=40
  echo "changed for segment $domnb"
 else
  lim1=50
  lim2=800
  fwhm=7
  ampl=40
  echo "changed for segment $domnb"
 fi
 
 if test -f "mDelila_raw_py_$domnb.spe" 
   then
      #get full data on calibration
      #/data/live/IT/tools/RecalEnergy -spe mDelila_raw_py_$domnb.spe -fmt A 16384 -lim $lim1 $lim2 -dwa $fwhm $ampl -60Co -v 2 -poly1 -ener 1460.82
      #/data/live/IT/tools/RecalEnergy -spe mDelila_raw_py_$domnb.spe -fmt A 16384 -lim $lim1 $lim2 -dwa $fwhm $ampl -60Co -v 2 -poly1 -ener 1460.82 > fulldata.calib
      ~/EliadeSorting/EliadeTools/RecalEnergy -spe mDelila_raw_py_$domnb.spe -fmt A 16384 -lim $lim1 $lim2 -dwa $fwhm $ampl -152Eu -40K -v 2 -poly2
      ~/EliadeSorting/EliadeTools/RecalEnergy -spe mDelila_raw_py_$domnb.spe -fmt A 16384 -lim $lim1 $lim2 -dwa $fwhm $ampl -152Eu -40K -v 2 -poly2 > fulldata.calib
      #delete first 14 lines
      awk 'NR > 14 { print }' fulldata.calib >  temp.calib
      
      echo -n "domain $domnb " >> eliade.calib      
      if grep 'Cal1' temp.calib
	then
 	   grep 'Cal1'  temp.calib >> eliade.calib      
 	   echo "Found Cal1"
 	else echo " ">> eliade.calib #if domain is emtry add a new line
	fi 	
	##################################      
	if grep '121.779' temp.calib
	then
           #echo -n "domain $domnb " >> 152Eu.calib    
 	  #grep '121.779'  temp.calib >> 152Eu.calib
 	   grep '121.779'  temp.calib >> res.temp
 	   grep '121.779'  temp.calib >> 152Eu.temp      
 	   echo "Found 121.779"
	fi
	##################################     
	if grep '244.693' temp.calib	    
	then
           #echo -n "domain $domnb " >> 152Eu.calib 
 	   grep '244.693'  temp.calib >> res.temp
 	   grep '244.693'  temp.calib >> 152Eu.temp      
 	   echo "Found"
	fi
	##################################
	if grep '344.272' temp.calib	    
	then
           #echo -n "domain $domnb " >> 152Eu.calib 
 	   grep '344.272'  temp.calib >> res.temp
 	   grep '344.272'  temp.calib >> 152Eu.temp      
 	   echo "Found"
	fi
	##################################     
	if grep '411.111' temp.calib	    
	then
           #echo -n "domain $domnb " >> 152Eu.calib 
 	   grep '411.111'  temp.calib >> res.temp
 	   grep '411.111'  temp.calib >> 152Eu.temp      
 	   echo "Found"
	fi
	##################################     
	if grep '443.979' temp.calib	    
	then
           #echo -n "domain $domnb " >> 152Eu.calib 
 	   grep '443.979'  temp.calib >> res.temp
 	   grep '443.979'  temp.calib >> 152Eu.temp      
 	   echo "Found"
	fi
	##################################     
	if grep '778.890' temp.calib	    
	then
           #echo -n "domain $domnb " >> 152Eu.calib 
 	   grep '778.890'  temp.calib >> res.temp
 	   grep '778.890'  temp.calib >> 152Eu.temp      
 	   echo "Found"
	fi
	##################################     
	if grep '964.014' temp.calib	    
	then
           #echo -n "domain $domnb " >> 152Eu.calib 
 	   grep '964.014'  temp.calib >> res.temp
 	   grep '964.014'  temp.calib >> 152Eu.temp      
 	   echo "Found"
	fi
	##################################     
	if grep '1085.793' temp.calib	    
	then
           #echo -n "domain $domnb " >> 152Eu.calib 
 	   grep '1085.793'  temp.calib >> res.temp
 	   grep '1085.793'  temp.calib >> 152Eu.temp      
 	   echo "Found"
	fi
	##################################     
	if grep '1112.070' temp.calib	    
	then
           #echo -n "domain $domnb " >> 152Eu.calib 
 	   grep '1112.070'  temp.calib >> res.temp
 	   grep '1112.070'  temp.calib >> 152Eu.temp      
 	   echo "Found"
	fi
	##################################     
	if grep '1407.993' temp.calib	    
	then
           #echo -n "domain $domnb " >> 152Eu.calib 
 	   grep '1407.993'  temp.calib >> res.temp
 	   grep '1407.993'  temp.calib >> 152Eu.temp      
 	   echo "Found"
	fi
	##################################          	

	if grep '#2' res.temp
	then
 	   grep '#2'  res.temp >> data.calib      
 	   echo "Found"
	fi            
   else 
     echo no file "mDelila_raw_py_$domnb.spe" 
     echo "domain $domnb #2     0  Slope = 0.000000    Cal1=[ -0.0000  1.000000 ]" >> eliade.calib       
     #echo "#2" >> eliade.calib 
     #sed -i "s/#2/$domnb/" eliade.calib
   fi  

sed -i "s/#2/$domnb/" data.calib
sed -i "s/#2/$domnb/" 152Eu.temp

 rm res.temp
 rm temp.calib	
 domnb=$(($domnb + 1))
done 
 
cat eliade.calib | sed 's/|/ /' | awk '{print $2, $9, $10}' >  eliade.coeff

#awk -F " " '{ print $1 " " $5 " " $9 " " $6 " " $7 " " $7/$6*$9 " keV"}' 152Eu.calib > resolution.txt

#awk -F " " '{ print $1 " " $5 " " $9 " " $6 " " $7 " " $7/$6*$9 " keV eff: " $5/'$activity'}' 152Eu.calib > resolution.txt

awk -F " " '{ print $1 " " $5 " " $9 " " $6 " " $7 " " $7/$6*$9 " keV" }' data.calib > resolution.txt 


awk '!/#/{print $0}' 152Eu.temp > 152Eu_1.temp # delete lines with "#"

awk -F " " '{ print $1 " " $q5 " " $9 " " $6 " " $7 " " $7/$6*$9 " [keV]; eff [%]: " $5/'$activity152Eu'*100}' 152Eu_1.temp > eff_152Eu.txt

rm 152Eu.temp
rm 152Eu_1.temp
rm data.calib

# case  $DIST in
# ##########################################################################################################
# "20")	if [ "$FIRST" -eq "0" ]; then 
# 	   FIRST=2000 LAST=2013
# 	fi
# 	
# 	echo "Will run the selector for run mass $MASS for runs from $FIRST up to $LAST"	
# 	
# 	domnb=$FIRST
# 	while test $domnb -le $LAST
# 	do
# 	
# 	eu_lut_file="LUT_EUCLIDES_w46_Testov_zero_time.dat";	
# 	eu_cut_file="euclides_cuts_w46_testov_300um_all.root";
# 	ga_lut_file="LUT_GALILEO_w46_Testov_run_2000.dat";
# 	
# 	
# 	if [ "$domnb" -eq "2000" ]; then 
# 	   eu_lut_file="LUT_EUCLIDES_w46_Testov_run_2000.dat";
# 	   eu_cut_file="euclides_cuts_run_2000_20um.root";
# 	fi
# 	
# 	if [ "$domnb" -ge "2001" ] && [ "$domnb" -le "2004" ]; then
# 	   eu_lut_file="LUT_EUCLIDES_w46_Testov_run_2001.dat";
# 	   eu_cut_file="euclides_cuts_run_2000_20um.root";
# 	fi
# 	
# 	if [ "$domnb" -eq "2005" ]; then
# 	   eu_lut_file="LUT_EUCLIDES_w46_Testov_run_2005.dat";
# 	   eu_cut_file="euclides_cuts_run_2013_20um.root";
# 	fi
# 	
# 	if [ "$domnb" -eq "2006" ]; then
# 	   eu_lut_file="LUT_EUCLIDES_w46_Testov_run_2006.dat";
# 	   eu_cut_file="euclides_cuts_run_2013_20um.root";
# 	fi
# 	
# 	if [ "$domnb" -eq "2007" ]; then
# 	   eu_lut_file="LUT_EUCLIDES_w46_Testov_run_2007.dat";
# 	   eu_cut_file="euclides_cuts_run_2013_20um.root";
# 	fi
# 	
# 	if [ "$domnb" -ge "2009" ]; then
# 	   eu_lut_file="LUT_EUCLIDES_w46_Testov_run_2013.dat";
# 	   eu_cut_file="euclides_cuts_run_2013_20um.root";
# 	fi
# 	
# 	echo "--------------------------------------------------------"
# 	echo "Distance is set to $DIST um file [$file_prefix$domnb.root"]
# 	echo "--------------------------------------------------------"	
# 	echo -ne "LUT file:: "; tput setaf 2; echo " $ga_lut_file";tput sgr0;
# 	echo -ne "LUT file:: "; tput setaf 2; echo " $eu_lut_file";tput sgr0;
# 	echo -ne "CUT file:: "; tput setaf 2; echo " $eu_cut_file";tput sgr0;
# 	echo -ne "DATA file  "; tput setaf 6; echo " $data_path$file_prefix$domnb.root";tput sgr0;
# 	echo "LUT dir::   $lut_path"
# 	echo "CUT dir::   $cut_path"
# 	if [ "$TEST" -eq "1" ]; then
# 	tput setaf 1; echo "Testing: 1000000 events...."; tput sgr0;    
# 	fi
# 	echo "--------------------------------------------------------"
# 
# 
# 	#echo "$cut_path""cut_charged_particles.root"
# 	#echo "$cut_path""$eu_cut_file" "$cut_path""cut_charged_particles.root"
# 	unlink "$cut_path""cut_charged_particles.root"
# 	ln -s "$cut_path$eu_cut_file" "$cut_path""cut_charged_particles.root"
# 
# 	unlink "$lut_path""LUT_EUCLIDES.dat"
# 	ln -s "$lut_path$eu_lut_file" "$lut_path""LUT_EUCLIDES.dat"
# 
# 	unlink "$lut_path""LUT_GALILEO.dat"
# 	ln -s "$lut_path$ga_lut_file" "$lut_path""LUT_GALILEO.dat"	
# 	
#         if test -f "$data_path$file_prefix$domnb.root"
# 	    then	
# 		echo $rootcommand
# 		echo starting "$inputpath/$prefix$domnb.root"        
# 		rootcommand=automated_selector.C"($domnb,$domnb,$DIST,$TEST)" 
# 		echo "rootcommand $rootcommand"
# 		if [ $FLAG -eq 1 ]; then	    
# 		    root -l -b -q $rootcommand
# 		    else 
# 		    echo "<<<bash-script test>>"
# 		fi
# 	    fi
# 	    domnb=$(($domnb + 1))
# 	done
# 	FINISHED=1;
# 	;;
# ##########################################################################################################	
# "100")	if [ "$FIRST" -eq "0" ]; then #100 um from November 2016; december 2016 is below
# 	   FIRST=2100 LAST=2102
# 	fi
# 	
# 	echo "Will run the selector for run mass $MASS for runs from $FIRST up to $LAST"	
# 	
# 	domnb=$FIRST
# 	while test $domnb -le $LAST
# 	do
# 	
# 	eu_lut_file="LUT_EUCLIDES_w46_Testov_100um.dat";
# 	eu_cut_file="euclides_cuts_w46_testov_300um_all.root";
# 	ga_lut_file="LUT_GALILEO_w46_Testov_100um.dat";	
# 	
# 	if [ "$domnb" -eq "2101" ]; then 
# 	   ga_lut_file="LUT_GALILEO_w46_Testov_run_2101.dat";
# 	fi	
# 	
# 	echo "--------------------------------------------------------"
# 	echo "Distance is set to $DIST um file [$file_prefix$domnb.root"]
# 	echo "--------------------------------------------------------"	
# 	echo -ne "LUT file:: "; tput setaf 2; echo " $ga_lut_file";tput sgr0;
# 	echo -ne "LUT file:: "; tput setaf 2; echo " $eu_lut_file";tput sgr0;
# 	echo -ne "CUT file:: "; tput setaf 2; echo " $eu_cut_file";tput sgr0;
# 	echo -ne "DATA file  "; tput setaf 6; echo " $data_path$file_prefix$domnb.root";tput sgr0;
# 	echo "LUT dir::   $lut_path"
# 	echo "CUT dir::   $cut_path"
# 	if [ "$TEST" -eq "1" ]; then 
# 
# 	tput setaf 1; echo "Testing: 1000000 events...."; tput sgr0;    
# 	fi
# 	echo "--------------------------------------------------------"
# 
# 
# 	#echo "$cut_path""cut_charged_particles.root"
# 	#echo "$cut_path""$eu_cut_file" "$cut_path""cut_charged_particles.root"
# 	unlink "$cut_path""cut_charged_particles.root"
# 	ln -s "$cut_path$eu_cut_file" "$cut_path""cut_charged_particles.root"
# 
# 	unlink "$lut_path""LUT_EUCLIDES.dat"
# 	ln -s "$lut_path$eu_lut_file" "$lut_path""LUT_EUCLIDES.dat"
# 
# 	unlink "$lut_path""LUT_GALILEO.dat"
# 	ln -s "$lut_path$ga_lut_file" "$lut_path""LUT_GALILEO.dat"	
# 	
#         if test -f "$data_path$file_prefix$domnb.root"
# 	    then	
# 		echo $rootcommand
# 		echo starting "$inputpath/$prefix$domnb.root"        
# 		rootcommand=automated_selector.C"($domnb,$domnb,$DIST,$TEST)" 
# 		echo "rootcommand $rootcommand"
# 		if [ $FLAG -eq 1 ]; then	    
# 		    root -l -b -q $rootcommand
# 		    else 
# 		    echo "<<<bash-script test>>"
# 		fi
# 	    fi
# 	    domnb=$(($domnb + 1))
# 	done
# 	FINISHED=1;
# 	;;
# ##########################################################################################################	
# "200")	if [ "$FIRST" -eq "0" ]; then 
# 	   FIRST=2200 LAST=2214 #run 2203 2205 2215 is activation 
# 	fi
# 	eu_lut_file="LUT_EUCLIDES_w46_Testov.dat";
# 	eu_cut_file="euclides_cuts_w46_testov_300um_all.root";
# 	ga_lut_file="LUT_GALILEO_w46_Testov.dat";
# 	
# 	
# 	if [ "$domnb" -ge "2200" ]; then % -ge
# 	   eu_lut_file="LUT_EUCLIDES_w46_Testov_run_2200.dat";
# 	   eu_cut_file="euclides_cuts_run_2013_20um.root";
# 	fi
# 	
# 	if [ "$domnb" -ge "2201" ]; then % -ge
# 	   eu_lut_file="LUT_EUCLIDES_w46_Testov_run_2201.dat";
# 	   eu_cut_file="euclides_cuts_run_2013_20um.root";
# 	fi
# 	
# 	if [ "$domnb" -ge "2202" ]; then % -ge
# 	   eu_lut_file="LUT_EUCLIDES_w46_Testov_run_2202.dat";
# 	   eu_cut_file="euclides_cuts_run_2013_20um.root";
# 	fi
# 	
# 	if [ "$domnb" -ge "2209" ]; then % -ge
# 	   eu_lut_file="LUT_EUCLIDES_w46_Testov_run_2210.dat";
# 	   eu_cut_file="euclides_cuts_run_2209_200um.root";
# 	fi
# 	
# 	if [ "$domnb" -ge "2210" ]; then % -ge
# 	   eu_lut_file="LUT_EUCLIDES_w46_Testov_run_2210.dat";
# 	   eu_cut_file="euclides_cuts_run_2210_200um.root";
# 	fi
# 	
# 	if [ "$domnb" -ge "2213" ]; then % -ge
# 	   eu_lut_file="LUT_EUCLIDES_w46_Testov_run_2213.dat";	  
# 	   eu_cut_file="euclides_cuts_run_2210_200um.root";#to be checked
# 	fi
# 	
# 	if [ "$domnb" -ge "2214" ]; then % -ge
# 	   eu_lut_file="LUT_EUCLIDES_w46_Testov_run_2214.dat";
# 	   eu_cut_file="euclides_cuts_run_2210_200um.root";#to be checked
# 	fi
# 	
# 	
# 	
# 	FINISHED=0;
# 	;;
# ##########################################################################################################	
# "300")	if [ "$FIRST" -eq "0" ]; then 
# 	   FIRST=2300 LAST=2314 #run 2302 2307 2308 2310 2311 are bad
# 	fi
# 	
# 	echo "Will run the selector for run mass $MASS for runs from $FIRST up to $LAST"	
# 	
# 	domnb=$FIRST
# 	while test $domnb -le $LAST
# 	do
# 	
# 	eu_lut_file="LUT_EUCLIDES_w46_Testov.dat";
# 	eu_cut_file="euclides_cuts_w46_testov_300um_all.root";
# 	ga_lut_file="LUT_GALILEO_w46_Testov.dat";
# 	
# 	if [ "$domnb" -eq "2300" ]; then 
# 	   eu_lut_file="LUT_EUCLIDES_w46_Testov_run_2300.dat";
# 	   eu_cut_file="eu_cuts_w46_testov_300um_run_2300.root";
#            #eu_cut_file="eu_cuts_w46_testov_300um_run_2300_1a_2a_1p.root";
# 	fi
# 	
# 	if [ "$domnb" -eq "2301" ]; then 
# 	   eu_lut_file="LUT_EUCLIDES_w46_Testov_run_2300.dat";
# 	   eu_cut_file="eu_cuts_w46_testov_300um_run_2300.root";
#            #eu_cut_file="eu_cuts_w46_testov_300um_run_2300_1a_2a_1p.root";
# 	fi
# 	
# 	if [ "$domnb" -eq "2303" ]; then 
# 	   eu_lut_file="LUT_EUCLIDES_w46_Testov_run_2300.dat";
# 	   eu_cut_file="eu_cuts_w46_testov_300um_run_2300.root";
#            #eu_cut_file="eu_cuts_w46_testov_300um_run_2300_1a_2a_1p.root";
# 	fi
# 	
# 	if [ "$domnb" -eq "2304" ]; then 
# 	   eu_lut_file="LUT_EUCLIDES_w46_Testov_run_2304.dat";
# 	   eu_cut_file="eu_cuts_w46_testov_300um_run_2300.root";
#            #eu_cut_file="eu_cuts_w46_testov_300um_run_2300_1a_2a_1p.root";
# 	fi
# 	
# 	if [ "$domnb" -eq "2305" ]; then 
# 	   eu_lut_file="LUT_EUCLIDES_w46_Testov_run_2305.dat";
#     	   eu_cut_file="eu_cuts_w46_testov_300um_run_2300.root";
#            #eu_cut_file="eu_cuts_w46_testov_300um_run_2300_1a_2a_1p.root";
# 	fi
# 	
# 	if [ "$domnb" -eq "2306" ]; then 
# 	   eu_lut_file="LUT_EUCLIDES_w46_Testov_run_2306.dat";
# 	   #eu_cut_file="eu_cuts_w46_testov_300um_run_2300_1a_2a_1p.root";
#   	   eu_cut_file="eu_cuts_w46_testov_300um_run_2300.root";
# #  	   eu_cut_file="eu_cuts_w46_testov_300um_run_2300_alpha_proton.root";
# 	fi
# 	
# 	if [ "$domnb" -eq "2309" ]; then 
# 	   eu_lut_file="LUT_EUCLIDES_w46_Testov_run_2306.dat";
# 	   eu_cut_file="eu_cuts_w46_testov_300um_run_2300.root";
#            #eu_cut_file="eu_cuts_w46_testov_300um_run_2300_1a_2a_1p.root";
# 	fi
# 	
# 	if [ "$domnb" -eq "2312" ]; then 
# 	   eu_lut_file="LUT_EUCLIDES_w46_Testov_run_2312.dat";
#  	   eu_cut_file="eu_cuts_w46_testov_300um_run_2300.root";
#            #eu_cut_file="eu_cuts_w46_testov_300um_run_2300_1a_2a_1p.root";
# 	fi
# 	
# 	if [ "$domnb" -eq "2313" ]; then 
# 	   eu_lut_file="LUT_EUCLIDES_w46_Testov_run_2313.dat";
# 	   #eu_cut_file="eu_cuts_w46_testov_300um_run_2313_1a_2a_1p.root";
#  	   eu_cut_file="eu_cuts_w46_testov_300um_run_2313.root";
# 	fi
# 	
# 	
# 	if [ "$domnb" -eq "2314" ]; then 
# 	   eu_lut_file="LUT_EUCLIDES_w46_Testov_run_2314.dat";
#  	   eu_cut_file="eu_cuts_w46_testov_300um_run_2314.root";
# 	   #eu_cut_file="eu_cuts_w46_testov_300um_run_2314_alpha_proton.root";
# 	   #eu_cut_file="eu_cuts_w46_testov_300um_run_2314_1a_2a_1p.root";
# 	fi
# 	
# 	echo "--------------------------------------------------------"
# 	echo "Distance is set to $DIST um file [$file_prefix$domnb.root"]
# 	echo "--------------------------------------------------------"	
# 	echo -ne "LUT file:: "; tput setaf 2; echo " $ga_lut_file";tput sgr0;
# 	echo -ne "LUT file:: "; tput setaf 2; echo " $eu_lut_file";tput sgr0;
# 	echo -ne "CUT file:: "; tput setaf 2; echo " $eu_cut_file";tput sgr0;
# 	echo -ne "DATA file  "; tput setaf 6; echo " $data_path$file_prefix$domnb.root";tput sgr0;
# 	echo "LUT dir::   $lut_path"
# 	echo "CUT dir::   $cut_path"
# 	if [ "$TEST" -eq "1" ]; then
# 	tput setaf 1; echo "Testing: 1000000 events...."; tput sgr0;    
# 	fi
# 	echo "--------------------------------------------------------"
# 
# 
# 	#echo "$cut_path""cut_charged_particles.root"
# 	#echo "$cut_path""$eu_cut_file" "$cut_path""cut_charged_particles.root"
# 	unlink "$cut_path""cut_charged_particles.root"
# 	ln -s "$cut_path$eu_cut_file" "$cut_path""cut_charged_particles.root"
# 
# 	unlink "$lut_path""LUT_EUCLIDES.dat"
# 	ln -s "$lut_path$eu_lut_file" "$lut_path""LUT_EUCLIDES.dat"
# 
# 	unlink "$lut_path""LUT_GALILEO.dat"
# 	ln -s "$lut_path$ga_lut_file" "$lut_path""LUT_GALILEO.dat"	
# 	
#         if test -f "$data_path$file_prefix$domnb.root"
# 	    then	
# 		echo $rootcommand
# 		echo starting "$inputpath/$prefix$domnb.root"        
# 		rootcommand=automated_selector.C"($domnb,$domnb,$DIST,$TEST)" 
# 		echo "rootcommand $rootcommand"
# 		if [ $FLAG -eq 1 ]; then	    
# 		    root -l -b -q $rootcommand
# 		    else 
# 		    echo "<<<bash-script test>>"
# 		fi
# 	    fi
# 	    domnb=$(($domnb + 1))
# 	done
# 	FINISHED=1;
# 	;;
# ##########################################################################################################	
# "1000")	if [ "$FIRST" -eq "0" ]; then 
# 	   FIRST=2400 LAST=2406 
# 	fi
# 	eu_lut_file="LUT_EUCLIDES_w46_Testov.dat";
# # 	eu_cut_file="euclides_cuts_w46_testov_300um_all.root";
# 	eu_cut_file="eu_cuts_w46_testov_1000um_run_2400.root";
# 	ga_lut_file="LUT_GALILEO_w46_Testov.dat";
# 	
# 	if [ "$domnb" -eq "2400" ]; then 	   
# 	   eu_cut_file="eu_cuts_w46_testov_1000um_run_2400.root";
#            eu_lut_file="LUT_EUCLIDES_w46_run_2400.dat";
# 	fi
# 	
# 	if [ "$domnb" -eq "2401" ]; then 	   
# 	   eu_cut_file="eu_cuts_w46_testov_1000um_run_2400.root";
#            eu_lut_file="LUT_EUCLIDES_w46_run_2400.dat";
# 	fi
# 	
# 	if [ "$domnb" -eq "2402" ]; then 	   
# 	   eu_cut_file="eu_cuts_w46_testov_1000um_run_2402.root";
#            eu_lut_file="LUT_EUCLIDES_w46_run_2402.dat";
# 	fi
# 	#runs 2403 2404 bad
# 	if [ "$domnb" -eq "2405" ]; then 	   
# 	   eu_cut_file="eu_cuts_w46_testov_1000um_run_2402.root";
#            eu_lut_file="LUT_EUCLIDES_w46_run_2405.dat";
# 	fi
# 	
# 	if [ "$domnb" -eq "2406" ]; then 	   
# 	   eu_cut_file="eu_cuts_w46_testov_1000um_run_2402.root";
#            eu_lut_file="LUT_EUCLIDES_w46_run_2406.dat";
# 	fi	
# 	FINISHED=0;
# 	;;
# ##########################################################################################################		
# "1500")	if [ "$FIRST" -eq "0" ]; then 
# 	   FIRST=2500 LAST=2503 #run 2501 is bad
# 	fi
# 	
# 	echo "Will run the selector for run mass $MASS for runs from $FIRST up to $LAST"	
# 	
# 	domnb=$FIRST
# 	while test $domnb -le $LAST
# 	do
# 	
# # 	eu_lut_file="LUT_EUCLIDES_w46_Testov.dat";
# # 	eu_cut_file="euclides_cuts_w46_testov_300um_all.root";
#         eu_lut_file="LUT_EUCLIDES_w46_run_2500.dat";
#         eu_cut_file="eu_cuts_w46_testov_1000um_run_2400.root";
# 	ga_lut_file="LUT_GALILEO_w46_Testov.dat";
# 	
# 	if [ "$domnb" -eq "2500" ]; then 	   
#            eu_lut_file="LUT_EUCLIDES_w46_run_2500_Testov.dat";
#            eu_cut_file="eu_cuts_w46_testov_1500um_run_2500.root";
# 	fi
# 	
# 	if [ "$domnb" -eq "2502" ]; then 
#  	   eu_lut_file="LUT_EUCLIDES_w46_run_2502_Testov.dat";
#  	   eu_cut_file="eu_cuts_w46_testov_1000um_run_2400.root";
# #	   eu_lut_file="LUT_EUCLIDES_w46_Testov_zero_time.dat";	
# 	fi	
# 	
# 	
# 	if [ "$domnb" -eq "2503" ]; then 
#  	   eu_lut_file="LUT_EUCLIDES_w46_run_2503_Testov.dat";
#  	   eu_cut_file="eu_cuts_w46_testov_1000um_run_2400.root";
# #	   eu_lut_file="LUT_EUCLIDES_w46_Testov_zero_time.dat";	
# 	fi	
# 	
# 	echo "--------------------------------------------------------"
# 	echo "Distance is set to $DIST um file [$file_prefix$domnb.root"]
# 	echo "--------------------------------------------------------"	
# 	echo -ne "LUT file:: "; tput setaf 2; echo " $ga_lut_file";tput sgr0;
# 	echo -ne "LUT file:: "; tput setaf 2; echo " $eu_lut_file";tput sgr0;
# 	echo -ne "CUT file:: "; tput setaf 2; echo " $eu_cut_file";tput sgr0;
# 	echo -ne "DATA file  "; tput setaf 6; echo " $data_path$file_prefix$domnb.root";tput sgr0;
# 	echo "LUT dir::   $lut_path"
# 	echo "CUT dir::   $cut_path"
# 	if [ "$TEST" -eq "1" ]; then 
# 
# 	tput setaf 1; echo "Testing: 1000000 events...."; tput sgr0;    
# 	fi
# 	echo "--------------------------------------------------------"
# 
# 
# 	#echo "$cut_path""cut_charged_particles.root"
# 	#echo "$cut_path""$eu_cut_file" "$cut_path""cut_charged_particles.root"
# 	unlink "$cut_path""cut_charged_particles.root"
# 	ln -s "$cut_path$eu_cut_file" "$cut_path""cut_charged_particles.root"
# 
# 	unlink "$lut_path""LUT_EUCLIDES.dat"
# 	ln -s "$lut_path$eu_lut_file" "$lut_path""LUT_EUCLIDES.dat"
# 
# 	unlink "$lut_path""LUT_GALILEO.dat"
# 	ln -s "$lut_path$ga_lut_file" "$lut_path""LUT_GALILEO.dat"	
# 	
#         if test -f "$data_path$file_prefix$domnb.root"
# 	    then	
# 		echo $rootcommand
# 		echo starting "$inputpath/$prefix$domnb.root"        
# 		rootcommand=automated_selector.C"($domnb,$domnb,$DIST,$TEST)" 
# 		echo "rootcommand $rootcommand"
# 		if [ $FLAG -eq 1 ]; then	    
# 		    root -l -b -q $rootcommand
# 		    else 
# 		    echo "<<<bash-script test>>"
# 		fi
# 	    fi
# 	    domnb=$(($domnb + 1))
# 	done
# 	FINISHED=1;
# 	;;
# ##########################################################################################################	
# #Decembre 2016
# ##########################################################################################################	
# "101")	if [ "$FIRST" -eq "0" ]; then  #100 um Decembre 2016
# 	   FIRST=3100 LAST=3118
# 	fi
# 	
# 	echo "Will run the selector for run mass $MASS for runs from $FIRST up to $LAST"	
# 	
# 	domnb=$FIRST
# 	while test $domnb -le $LAST
# 	do
# 	BADFILE=0;
#         eu_lut_file="LUT_EUCLIDES_w49_Testov_run_3103.dat";
# 	#eu_cut_file="euclides_cuts_w46_testov_300um_all.root";
# 	eu_cut_file="euclides_cuts_run_3200_700um.root";
# 	ga_lut_file="LUT_GALILEO_w49_Testov_new.dat";
# 	
# 	if [ "$domnb" -eq "3102" ]; then
# 	   BADFILE=1;
#  	   eu_lut_file="LUT_EUCLIDES_w49_Testov_run_3102.dat";   
# 	   #eu_lut_file="LUT_EUCLIDES_w49_Testov_zero_time.dat";	
# 	fi
# 	
# 	if [ "$domnb" -eq "3105" ]; then 
#  	   eu_lut_file="LUT_EUCLIDES_w49_Testov_run_3105.dat";   
# 	   #eu_lut_file="LUT_EUCLIDES_w49_Testov_zero_time.dat";	
# 	fi	
# 	
# 	if [ "$domnb" -eq "3106" ]; then 
#  	   eu_lut_file="LUT_EUCLIDES_w49_Testov_run_3106.dat";   
# 	   #eu_lut_file="LUT_EUCLIDES_w49_Testov_zero_time.dat";	
# 	fi
# 	
# 	if [ "$domnb" -eq "3109" ]; then 
#  	   eu_lut_file="LUT_EUCLIDES_w49_Testov_run_3109.dat";   
# 	   #eu_lut_file="LUT_EUCLIDES_w49_Testov_zero_time.dat";	
# 	fi
# 	
# 	if [ "$domnb" -eq "3110" ]; then 
#  	   eu_lut_file="LUT_EUCLIDES_w49_Testov_run_3109.dat";   
# 	   #eu_lut_file="LUT_EUCLIDES_w49_Testov_zero_time.dat";	
# 	fi
# 	
# 	if [ "$domnb" -eq "3111" ]; then 
#  	   eu_lut_file="LUT_EUCLIDES_w49_Testov_run_3109.dat";   
# 	   #eu_lut_file="LUT_EUCLIDES_w49_Testov_zero_time.dat";	
# 	fi
# 	
# 	if [ "$domnb" -eq "3112" ]; then 
#  	   eu_lut_file="LUT_EUCLIDES_w49_Testov_run_3109.dat";   
# 	   #eu_lut_file="LUT_EUCLIDES_w49_Testov_zero_time.dat";	
# 	fi
# 	
# 	if [ "$domnb" -eq "3113" ]; then 
#  	   eu_lut_file="LUT_EUCLIDES_w49_Testov_run_3113.dat";   
# 	   #eu_lut_file="LUT_EUCLIDES_w49_Testov_zero_time.dat";	
# 	fi
# 	
# 	if [ "$domnb" -eq "3114" ]; then 
#  	   eu_lut_file="LUT_EUCLIDES_w49_Testov_run_3114.dat";   
# 	   #eu_lut_file="LUT_EUCLIDES_w49_Testov_zero_time.dat";	
# 	fi
# 	
# 	if [ "$domnb" -ge "3115" ]; then 
#  	   eu_lut_file="LUT_EUCLIDES_w49_Testov_run_3117.dat";   
# 	   #eu_lut_file="LUT_EUCLIDES_w49_Testov_zero_time.dat";	
# 	fi	
# 	
# 	echo "--------------------------------------------------------"
# 	echo "Distance is set to $DIST um file [$file_prefix$domnb.root"]
# 	echo "--------------------------------------------------------"	
# 	echo -ne "LUT file:: "; tput setaf 2; echo " $ga_lut_file";tput sgr0;
# 	echo -ne "LUT file:: "; tput setaf 2; echo " $eu_lut_file";tput sgr0;
# 	echo -ne "CUT file:: "; tput setaf 2; echo " $eu_cut_file";tput sgr0;
# 	echo -ne "DATA file  "; tput setaf 6; echo " $data_path$file_prefix$domnb.root";tput sgr0;
# 	echo "LUT dir::   $lut_path"
# 	echo "CUT dir::   $cut_path"
# 	if [ "$TEST" -eq "1" ]; then 
# 
# 	tput setaf 1; echo "Testing: 1000000 events...."; tput sgr0;    
# 	fi
# 	echo "--------------------------------------------------------"
# 	
# 	#if [ "$BADFILE" -eq "1" ]; then 
# 	#    continue
# 	#fi
# 
# 	#echo "$cut_path""cut_charged_particles.root"
# 	#echo "$cut_path""$eu_cut_file" "$cut_path""cut_charged_particles.root"
# 	unlink "$cut_path""cut_charged_particles.root"
# 	ln -s "$cut_path$eu_cut_file" "$cut_path""cut_charged_particles.root"
# 
# 	unlink "$lut_path""LUT_EUCLIDES.dat"
# 	ln -s "$lut_path$eu_lut_file" "$lut_path""LUT_EUCLIDES.dat"
# 
# 	unlink "$lut_path""LUT_GALILEO.dat"
# 	ln -s "$lut_path$ga_lut_file" "$lut_path""LUT_GALILEO.dat"	
# 	
#         if test -f "$data_path$file_prefix$domnb.root"
# 	    then	
# 		echo $rootcommand
# 		echo starting "$inputpath/$prefix$domnb.root"        
# 		rootcommand=automated_selector.C"($domnb,$domnb,$DIST,$TEST)" 
# 		echo "rootcommand $rootcommand"
# 		if [ $FLAG -eq 1 ]; then	    
# 		    root -l -b -q $rootcommand
# 		    else 
# 		    echo "<<<bash-script test>>"
# 		fi
# 	    fi
# 	    domnb=$(($domnb + 1))
# 	done
# 	FINISHED=1;
# 	;;
# ##########################################################################################################		
# "700")	if [ "$FIRST" -eq "0" ]; then  
# 	   FIRST=3201 LAST=3203
# 	fi
# 	eu_lut_file="LUT_EUCLIDES_w49_Testov_run_3106.dat";
# 	eu_cut_file="euclides_cuts_run_3200_700um.root";
# 	ga_lut_file="LUT_GALILEO_w49_Testov_new.dat";
# 	
# 	if [ "$domnb" -eq "3202" ]; then 
#  	   eu_lut_file="LUT_EUCLIDES_w49_Testov_run_3203.dat";   	   
# 	fi	
# 	
# 	if [ "$domnb" -eq "3203" ]; then 
#  	   eu_lut_file="LUT_EUCLIDES_w49_Testov_run_3203.dat";   	   
# 	fi
# 	
# 	FINISHED=0;
# 	;;
# ##########################################################################################################	
# #Blackout
# ##########################################################################################################
# "10")	if [ "$FIRST" -eq "0" ]; then  
# 	   FIRST=3400 LAST=3408
# 	fi
# 	#eu_lut_file="LUT_EUCLIDES_w49_Testov_run_3400.dat";
# 	eu_lut_file="LUT_EUCLIDES_w49_Testov_run_3900.dat";
# 	eu_cut_file="euclides_cuts_w46_testov_300um_all.root";
# 	ga_lut_file="LUT_GALILEO_w49_Testov_run_3400.dat";
# 	FINISHED=0;
# 	;;
# ##########################################################################################################
# "60")	if [ "$FIRST" -eq "0" ]; then  
# 	   FIRST=3500 LAST=3505
# 	fi
# 	#eu_lut_file="LUT_EUCLIDES_w49_Testov_run_3400.dat";
# 	eu_lut_file="LUT_EUCLIDES_w49_Testov_run_3900.dat";
# 	eu_cut_file="euclides_cuts_w46_testov_300um_all.root";
# 	ga_lut_file="LUT_GALILEO_w49_Testov_run_3400.dat";
# 	FINISHED=0;
# 	;;
# ##########################################################################################################
# "150")	if [ "$FIRST" -eq "0" ]; then  
# 	   FIRST=3600 LAST=3605
# 	fi
# 	#eu_lut_file="LUT_EUCLIDES_w49_Testov_run_3400.dat";
# 	eu_lut_file="LUT_EUCLIDES_w49_Testov_run_3900.dat";
# 	eu_cut_file="euclides_cuts_w46_testov_300um_all.root";
# 	ga_lut_file="LUT_GALILEO_w49_Testov_run_3400.dat";
# 	FINISHED=0;
# 	;;
# ##########################################################################################################
# "18")	if [ "$FIRST" -eq "0" ]; then  
# 	   FIRST=3700 LAST=3702
# 	fi
# 	#eu_lut_file="LUT_EUCLIDES_w49_Testov_run_3400.dat";
# 	eu_lut_file="LUT_EUCLIDES_w49_Testov_run_3400.dat";
# 	eu_cut_file="euclides_cuts_w46_testov_300um_all.root";
# 	ga_lut_file="LUT_GALILEO_w49_Testov_run_3400.dat";
# 	FINISHED=0;
# 	;;
# ##########################################################################################################
# "2000")	if [ "$FIRST" -eq "0" ]; then  
# 	   FIRST=3800 LAST=3803
# 	fi
# 	#eu_lut_file="LUT_EUCLIDES_w49_Testov_run_3400.dat";
# 	eu_lut_file="LUT_EUCLIDES_w49_Testov_run_3900.dat";   	   
# 	eu_cut_file="euclides_cuts_run_3900_130um.root";
# 	ga_lut_file="LUT_GALILEO_w49_Testov_run_3400.dat";
# 	FINISHED=0;
# 	;;
# ##########################################################################################################
# "130")	if [ "$FIRST" -eq "0" ]; then  
# 	   FIRST=3900 LAST=3904
# 	fi		
# 	
# 	eu_lut_file="LUT_EUCLIDES_w49_Testov_run_3900.dat";
# 	
# # 	eu_cut_file="euclides_cuts_w46_testov_300um_all.root";
#  	eu_cut_file="euclides_cuts_run_3900_130um.root";
#  	
# 	ga_lut_file="LUT_GALILEO_w49_Testov_run_3400.dat";
# 	
# 	
# 	if [ "$domnb" -ge "3900" ]; then 
#  	   eu_lut_file="LUT_EUCLIDES_w49_Testov_run_3900.dat";   	   
# 	fi	
# 	
# 	
# 	FINISHED=0;
# 	;;
# esac
# 
# if [ "$FINISHED" -eq "0" ]; then 
# 
# echo "Will run the selector for run mass $MASS for runs from $FIRST up to $LAST"	
# 	
# 	domnb=$FIRST
# 	while test $domnb -le $LAST
# 	do	
# 	
# 	echo "--------------------------------------------------------"
# 	echo "Distance is set to $DIST um file [$file_prefix$domnb.root"]
# 	echo "--------------------------------------------------------"	
# 	echo -ne "LUT file:: "; tput setaf 2; echo " $ga_lut_file";tput sgr0;
# 	echo -ne "LUT file:: "; tput setaf 2; echo " $eu_lut_file";tput sgr0;
# 	echo -ne "CUT file:: "; tput setaf 2; echo " $eu_cut_file";tput sgr0;
# 	echo -ne "DATA file  "; tput setaf 6; echo " $data_path$file_prefix$domnb.root";tput sgr0;
# 	echo "LUT dir::   $lut_path"
# 	echo "CUT dir::   $cut_path"
# 	if [ "$TEST" -eq "1" ]; then 
# 
# 	tput setaf 1; echo "Testing: 1000000 events...."; tput sgr0;    
# 	fi
# 	echo "--------------------------------------------------------"
# 	
# 	unlink "$cut_path""cut_charged_particles.root"
# 	ln -s "$cut_path$eu_cut_file" "$cut_path""cut_charged_particles.root"
# 
# 	unlink "$lut_path""LUT_EUCLIDES.dat"
# 	ln -s "$lut_path$eu_lut_file" "$lut_path""LUT_EUCLIDES.dat"
# 
# 	unlink "$lut_path""LUT_GALILEO.dat"
# 	ln -s "$lut_path$ga_lut_file" "$lut_path""LUT_GALILEO.dat"	
# 	
#         if test -f "$data_path$file_prefix$domnb.root"
# 	    then	
# 		echo $rootcommand
# 		echo starting "$inputpath/$prefix$domnb.root"        
# 		rootcommand=automated_selector.C"($domnb,$domnb,$DIST,$TEST)" 
# 		echo "rootcommand $rootcommand"
# 		if [ $FLAG -eq 1 ]; then	    
# 		    root -l -b -q $rootcommand
# 		    else 
# 		    echo "<<<bash-script test>>"
# 		fi
# 	    fi
# 	    domnb=$(($domnb + 1))
# 	done
# 	FINISHED=1;
# fi
# 
# 
# 
# 
# 
