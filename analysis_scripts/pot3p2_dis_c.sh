#!/bin/bash
cd ..
for ip in `seq 0.2 0.2 0.6`
do
	echo "Starting 3p2 POT scan with POT of 10^"$ip
	./sbnfit --pot $ip --num 2 --dis >> fractiondata/pot3p2/"3p2_pot_"$ip"_dis.dat"
done
cd analysis_scripts
