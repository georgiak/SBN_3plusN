#!/bin/bash
cd ..
rm fractiondata/NUBAR_MODE/pot3p1/3p1_pot_*

for ip in `seq -2.0 0.25 0.0`
do
for ipb in `seq -2.0 0.25 0.0`
do

	echo "Starting 3p1 POT scan with POT of 10^"$ip" and POTBAR 10^"$ipb
	./sbnfit --pot $ip:$ipb --num 1 --both --anti >> fractiondata/NUBAR_MODE/pot3p1/"3p1_pot_"$ip"_"$ipb"_both.bar.dat"
	./sbnfit --pot $ip:$ipb --num 1 --app --anti  >> fractiondata/NUBAR_MODE/pot3p1/"3p1_pot_"$ip"_"$ipb"_app.bar.dat"
	./sbnfit --pot $ip:$ipb --num 1 --dis --anti >> fractiondata/NUBAR_MODE/pot3p1/"3p1_pot_"$ip"_"$ipb"_dis.bar.dat"
done
done
cd analysis_scripts
