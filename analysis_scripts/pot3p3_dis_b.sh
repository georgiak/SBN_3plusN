#!/bin/bash
cd ..
for ip in `seq -1.4 0.2 0.0`
do
	echo "Starting 3p3 POT scan with POT of 10^"$ip
	./sbnfit --pot $ip --num 3 --dis >> fractiondata/pot3p3/"3p3_pot_"$ip"_dis.dat"
done
cd analysis_scripts
