#!/bin/bash
cd ..
for ip in `seq -3.0 0.1 0.0`
do
	echo "Starting 3p1 POT scan with POT of 10^"$ip
	./sbnfit --pot $ip --num 1 >> fractiondata/pot3p1/"3p1_pot_"$ip"_both.dat"
done
cd analysis_scripts
