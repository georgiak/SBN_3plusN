#!/bin/bash
cd ..

rm fractiondata/basic_analysed_3p3.dat

#run this once to produce grepped files
grep "#" fractiondata/3p3_app.dat > fractiondata/3p3_app.dat.grepped
grep "#" fractiondata/3p3_dis.dat > fractiondata/3p3_dis.dat.grepped
grep "#" fractiondata/3p3_both.dat > fractiondata/3p3_both.dat.grepped

for chi in `seq 0 0.1 2.5 `
do
		TOTALD=$(wc -l < fractiondata/3p3_dis.dat.grepped)
		TOTALB=$(wc -l < fractiondata/3p3_both.dat.grepped)
		TOTALA=$(wc -l < fractiondata/3p3_app.dat.grepped)
 		echo "On chi "$chi
		ANSD=$(awk -v mychi="$chi" -v count=0 '{ if ($17 >= 10^mychi) count++}  END {print count}' fractiondata/3p3_dis.dat.grepped)
		ANSB=$(awk -v mychi="$chi" -v count=0 '{ if ($17 >= 10^mychi) count++}  END {print count}' fractiondata/3p3_both.dat.grepped)
		ANSA=$(awk -v mychi="$chi" -v count=0 '{ if ($17 >= 10^mychi) count++}  END {print count}' fractiondata/3p3_app.dat.grepped)
		echo 3p3 $chi $ANSD $ANSB $ANSA $TOTALD $TOTALB $TOTALA  >> fractiondata/basic_analysed_3p3.dat

done


cd analysis_scripts
