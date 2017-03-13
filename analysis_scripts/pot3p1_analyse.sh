#!/bin/bash
mkdir -p fractiondata/pot3p1
cd ..

rm fractiondata/pot3p1/analysed_3p1_both.dat
rm fractiondata/pot3p1/analysed_3p1_dis.dat
rm fractiondata/pot3p1/analysed_3p1_app.dat
rm fractiondata/pot3p1/analysed_3p1_wierd.dat

CHI3=14.1564
CHI5=31.8121

for ip in `seq -4.0 0.1 0.5`
do
	TOTAL=$(wc -l < fractiondata/pot3p1/"3p1_pot_"$ip"_both.dat" )
	echo "On Both POT"$ip
		ANS3=$(awk -v mychi="$CHI3" -v count=0 '{ if ($17 >= mychi) count++}  END {print count}' fractiondata/pot3p1/"3p1_pot_"$ip"_both.dat")
		ANS5=$(awk -v mychi="$CHI5" -v count=0 '{ if ($17 >= mychi) count++}  END {print count}' fractiondata/pot3p1/"3p1_pot_"$ip"_both.dat")
		echo 3p1 $ip $CHI3 $ANS3 $CHI5 $ANS5 $TOTAL  >> fractiondata/pot3p1/analysed_3p1_both.dat

done

for ip in `seq -4.0 0.1 0.5`
do
	TOTAL=$(wc -l < fractiondata/pot3p1/"3p1_pot_"$ip"_app.dat" )
	echo "On App POT"$ip
		ANS3=$(awk -v mychi="$CHI3" -v count=0 '{ if ($17 >= mychi) count++}  END {print count}' fractiondata/pot3p1/"3p1_pot_"$ip"_app.dat")
		ANS5=$(awk -v mychi="$CHI5" -v count=0 '{ if ($17 >= mychi) count++}  END {print count}' fractiondata/pot3p1/"3p1_pot_"$ip"_app.dat")
		echo 3p1 $ip $CHI3 $ANS3 $CHI5 $ANS5 $TOTAL  >> fractiondata/pot3p1/analysed_3p1_app.dat
	done


for ip in `seq -4.0 0.1 0.5`
do
	TOTAL=$(wc -l < fractiondata/pot3p1/"3p1_pot_"$ip"_dis.dat" )
	echo "On Dis POT"$ip
		ANS3=$(awk -v mychi="$CHI3" -v count=0 '{ if ($17 >= mychi) count++}  END {print count}' fractiondata/pot3p1/"3p1_pot_"$ip"_dis.dat")
		ANS5=$(awk -v mychi="$CHI5" -v count=0 '{ if ($17 >= mychi) count++}  END {print count}' fractiondata/pot3p1/"3p1_pot_"$ip"_dis.dat")
		echo 3p1 $ip $CHI3 $ANS3 $CHI5 $ANS5 $TOTAL  >> fractiondata/pot3p1/analysed_3p1_dis.dat
	done

for ip in `seq -4.0 0.1 0.5`
do
	TOTAL=$(wc -l < fractiondata/pot3p1/"3p1_pot_"$ip"_wierd.dat" )
	echo "On wierd POT"$ip
		ANS3=$(awk -v mychi="$CHI3" -v count=0 '{ if ($17 >= mychi) count++}  END {print count}' fractiondata/pot3p1/"3p1_pot_"$ip"_wierd.dat")
		ANS5=$(awk -v mychi="$CHI5" -v count=0 '{ if ($17 >= mychi) count++}  END {print count}' fractiondata/pot3p1/"3p1_pot_"$ip"_wierd.dat")
		echo 3p1 $ip $CHI3 $ANS3 $CHI5 $ANS5 $TOTAL  >> fractiondata/pot3p1/analysed_3p1_wierd.dat
	done





cd analysis_scripts
