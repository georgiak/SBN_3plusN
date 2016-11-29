#!/bin/bash
cd ..
rm fractiondata/pot3p1/3p1_pot_*
for ip in `seq -4.0 0.1 0.5`
do
	echo "Starting 3p1 POT scan with POT of 10^"$ip
	./sbnfit2 --pot $ip --num 1 --both >> fractiondata/pot3p1/"3p1_pot_"$ip"_both.dat"
	./sbnfit2 --pot $ip --num 1 --app >> fractiondata/pot3p1/"3p1_pot_"$ip"_app.dat"
	./sbnfit2 --pot $ip --num 1 --dis >> fractiondata/pot3p1/"3p1_pot_"$ip"_dis.dat"
	./sbnfit2 --pot $ip --num 1 --wierd >> fractiondata/pot3p1/"3p1_pot_"$ip"_wierd.dat"
done
cd analysis_scripts
