#!/bin/bash
cd ..
for ip in `seq -3.0 0.2 -1.6`
do
	echo "Starting 3p3 POT scan with POT of 10^"$ip
	./sbnfit --pot $ip --num 3 --app >> fractiondata/pot3p3/"3p3_pot_"$ip"_app.dat"
done
cd analysis_scripts
