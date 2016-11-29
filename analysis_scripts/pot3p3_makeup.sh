#!/bin/bash
cd ..
for ip in `seq 0.2 0.2 0.6`
do
	./sbnfit --pot $ip --num 3 --app >> fractiondata/pot3p3/"3p3_pot_"$ip"_app.dat" &
	./sbnfit --pot $ip --num 3 --both >> fractiondata/pot3p3/"3p3_pot_"$ip"_both.dat" &
	./sbnfit --pot $ip --num 3 --dis >> fractiondata/pot3p3/"3p3_pot_"$ip"_dis.dat" &
done
cd analysis_scripts
