#!/bin/bash
cd ..

rm fractiondata/basic_analysed_3p2.dat

#run this once to produce grepped files
grep "#" fractiondata/3p2_app.dat > fractiondata/3p2_app.dat.grepped
grep "#" fractiondata/3p2_dis.dat > fractiondata/3p2_dis.dat.grepped
grep "#" fractiondata/3p2_both.dat > fractiondata/3p2_both.dat.grepped
grep "#" fractiondata/3p2_wierd.dat > fractiondata/3p2_wierd.dat.grepped

for chi in `seq 0 0.1 2.5 `
do
		TOTAL=$(wc -l < fractiondata/3p2_both.dat.grepped)
 		echo "On chi "$chi
		ANSD=$(awk -v mychi="$chi" -v count=0 '{ if ($17 >= 10^mychi) count++}  END {print count}' fractiondata/3p2_dis.dat.grepped)
		ANSB=$(awk -v mychi="$chi" -v count=0 '{ if ($17 >= 10^mychi) count++}  END {print count}' fractiondata/3p2_both.dat.grepped)
		ANSA=$(awk -v mychi="$chi" -v count=0 '{ if ($17 >= 10^mychi) count++}  END {print count}' fractiondata/3p2_app.dat.grepped)
		ANSW=$(awk -v mychi="$chi" -v count=0 '{ if ($17 >= 10^mychi) count++}  END {print count}' fractiondata/3p2_wierd.dat.grepped)
		echo 3p2 $chi $ANSD $ANSB $ANSA $TOTAL $ANSW  >> fractiondata/basic_analysed_3p2.dat

done


cd analysis_scripts
