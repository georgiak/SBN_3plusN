#!/bin/bash

cd ..
rm fractiondata/pot3p1/analysed_3p1_both.dat
for ip in `seq -3.0 0.1 0.0`
do
	TOTAL=$(wc -l < fractiondata/pot3p1/"3p1_pot_"$ip"_both.dat" )
	echo "On POT"$ip
	for chi in `seq 1 2 100`
	do
		ANS=$(awk -v mychi="$chi" '{ if ($17 <= mychi) count++}  END {print count}' fractiondata/pot3p1/"3p1_pot_"$ip"_both.dat")
		PERCENT=$(calc $ANS/$TOTAL)
		echo 3p1 $ip $chi $ANS $TOTAL $PERCENT >> fractiondata/pot3p1/analysed_3p1_both.dat
	done

done
cd analysis_scripts
