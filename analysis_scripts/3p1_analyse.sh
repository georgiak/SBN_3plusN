#!/bin/bash
cd ..

rm fractiondata/basic_analysed_3p1.dat

grep "#" fractiondata/3p1_app.dat > fractiondata/3p1_app.dat.grepped
grep "#" fractiondata/3p1_dis.dat > fractiondata/3p1_dis.dat.grepped
grep "#" fractiondata/3p1_both.dat > fractiondata/3p1_both.dat.grepped
grep "#" fractiondata/3p1_wierd.dat > fractiondata/3p1_wierd.dat.grepped

for chi in `seq -1.0 0.1 2.5 `
do
		TOTAL=$(wc -l < fractiondata/3p1_both.dat.grepped)
 		echo "On Both POT"$ip
		ANSD=$(awk -v mychi=$chi -v count=0 '{ if ($17 >= 10^mychi) count++}  END {print count}' fractiondata/3p1_dis.dat.grepped)
		ANSB=$(awk -v mychi=$chi -v count=0 '{ if ($17 >= 10^mychi) count++}  END {print count}' fractiondata/3p1_both.dat.grepped)
		ANSA=$(awk -v mychi=$chi -v count=0 '{ if ($17 >= 10^mychi) count++}  END {print count}' fractiondata/3p1_app.dat.grepped)
		ANSW=$(awk -v mychi=$chi -v count=0 '{ if ($17 >= 10^mychi) count++}  END {print count}' fractiondata/3p1_wierd.dat.grepped)
		echo 3p1 $chi $ANSD $ANSB $ANSA $TOTAL $ANSW  >> fractiondata/basic_analysed_3p1.dat

done


cd analysis_scripts
