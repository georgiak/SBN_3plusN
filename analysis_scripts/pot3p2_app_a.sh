#!/bin/bash
cd ..
for ip in `seq -3.0 0.2 -1.6`
do
	echo "Starting 3p2 POT scan with POT of 10^"$ip
	./sbnfit --pot $ip --num 2 --app >> fractiondata/pot3p2/"3p2_pot_"$ip"_app.dat"
done
cd analysis_scripts
