#!/bin/bash
cd ..

rm fractiondata/basic_analysed_3p3.dat

#run this once to produce grepped files
#grep "#" fractiondata/3p3_app.dat > fractiondata/3p3_app.dat.grepped
#grep "#" fractiondata/3p3_dis.dat > fractiondata/3p3_dis.dat.grepped
#grep "#" fractiondata/3p3_both.dat > fractiondata/3p3_both.dat.grepped
grep "#" fractiondata/3p3_wierd.dat > fractiondata/3p3_wierd.dat.grepped

for chi in `seq 0 0.1 2.5 `
do
		TOTALD=$(wc -l < fractiondata/3p3_dis.dat.grepped)
		TOTALB=$(wc -l < fractiondata/3p3_both.dat.grepped)
		TOTALA=$(wc -l < fractiondata/3p3_app.dat.grepped)
		TOTALW=$(wc -l < fractiondata/3p3_wierd.dat.grepped)
		echo "On chi "$chi
		ANSD=$(awk -v mychi="$chi" -v count=0 '{ if ($17 >= 10^mychi) count++}  END {print count}' fractiondata/3p3_dis.dat.grepped)
		ANSB=$(awk -v mychi="$chi" -v count=0 '{ if ($17 >= 10^mychi) count++}  END {print count}' fractiondata/3p3_both.dat.grepped)
		ANSA=$(awk -v mychi="$chi" -v count=0 '{ if ($17 >= 10^mychi) count++}  END {print count}' fractiondata/3p3_app.dat.grepped)
		ANSW=$(awk -v mychi="$chi" -v count=0 '{ if ($17 >= 10^mychi) count++}  END {print count}' fractiondata/3p3_wierd.dat.grepped)
		echo 3p3 $chi $ANSD $ANSB $ANSA $TOTALD $TOTALB $TOTALA $ANSW $TOTALW  >> fractiondata/basic_analysed_3p3.dat

done

#for chi in `seq 0 0.1 2.5 `
#do
#		TOTALD=$(wc -l < fractiondata/NUBAR_MODE/pot3p3/3p3_pot_0.00_-2.00_dis.bar.dat)
#		TOTALB=$(wc -l < fractiondata/NUBAR_MODE/pot3p3/3p3_pot_0.00_-2.00_both.bar.dat)
#		TOTALA=$(wc -l < fractiondata/NUBAR_MODE/pot3p3/3p3_pot_0.00_-1.00_app.bar.dat)
 #		echo "On chi "$chi
#		ANSD=$(awk -v mychi="$chi" -v count=0 '{ if ($17 >= 10^mychi) count++}  END {print count}' fractiondata/NUBAR_MODE/pot3p3/3p3_pot_0.00_-2.00_dis.bar.dat)
#		ANSB=$(awk -v mychi="$chi" -v count=0 '{ if ($17 >= 10^mychi) count++}  END {print count}' fractiondata/NUBAR_MODE/pot3p3/3p3_pot_0.00_-2.00_both.bar.dat)
#		ANSA=$(awk -v mychi="$chi" -v count=0 '{ if ($17 >= 10^mychi) count++}  END {print count}' fractiondata/NUBAR_MODE/pot3p3/3p3_pot_0.00_-1.00_app.bar.dat)
#		echo 3p3 $chi $ANSD $ANSB $ANSA $TOTALD $TOTALB $TOTALA  >> fractiondata/basic_analysed_3p3.dat
#
#done



cd analysis_scripts
