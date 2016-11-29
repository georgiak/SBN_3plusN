#!/bin/bash
cd ..
for ip in `seq -1.4 0.2 0.6`
do
	echo "Starting 3p2 POT scan with POT of 10^"$ip
	./sbnfit2 --pot $ip --num 2 --wierd >> fractiondata/pot3p2/"3p2_pot_"$ip"_wierd.dat"
done
cd analysis_scripts
