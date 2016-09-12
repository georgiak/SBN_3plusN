#!/bin/bash
mkdir -p fractiondata/pot3p2
cd ..

rm fractiondata/pot3p2/analysed_3p2_both.dat
rm fractiondata/pot3p2/analysed_3p2_dis.dat
rm fractiondata/pot3p2/analysed_3p2_app.dat

CHI3=21.8466
CHI5=41.7798

for ip in `seq -3.0 0.2 0.0`
do
	TOTAL=$(wc -l < fractiondata/pot3p2/"3p2_pot_"$ip"_both.dat" )
	echo "On Both POT"$ip
		ANS3=$(awk -v mychi="$CHI3" -v count=0 '{ if ($17 >= mychi) count++}  END {print count}' fractiondata/pot3p2/"3p2_pot_"$ip"_both.dat")
		ANS5=$(awk -v mychi="$CHI5" -v count=0 '{ if ($17 >= mychi) count++}  END {print count}' fractiondata/pot3p2/"3p2_pot_"$ip"_both.dat")
		echo 3p2 $ip $CHI3 $ANS3 $CHI5 $ANS5 $TOTAL  >> fractiondata/pot3p2/analysed_3p2_both.dat

done

for ip in `seq -3.0 0.2 0.0`
do
	TOTAL=$(wc -l < fractiondata/pot3p2/"3p2_pot_"$ip"_app.dat" )
	echo "On App POT"$ip
		ANS3=$(awk -v mychi="$CHI3" -v count=0 '{ if ($17 >= mychi) count++}  END {print count}' fractiondata/pot3p2/"3p2_pot_"$ip"_app.dat")
		ANS5=$(awk -v mychi="$CHI5" -v count=0 '{ if ($17 >= mychi) count++}  END {print count}' fractiondata/pot3p2/"3p2_pot_"$ip"_app.dat")
		echo 3p2 $ip $CHI3 $ANS3 $CHI5 $ANS5 $TOTAL  >> fractiondata/pot3p2/analysed_3p2_app.dat
	done


for ip in `seq -3.0 0.2 0.0`
do
	TOTAL=$(wc -l < fractiondata/pot3p2/"3p2_pot_"$ip"_dis.dat" )
	echo "On Dis POT"$ip
		ANS3=$(awk -v mychi="$CHI3" -v count=0 '{ if ($17 >= mychi) count++}  END {print count}' fractiondata/pot3p2/"3p2_pot_"$ip"_dis.dat")
		ANS5=$(awk -v mychi="$CHI5" -v count=0 '{ if ($17 >= mychi) count++}  END {print count}' fractiondata/pot3p2/"3p2_pot_"$ip"_dis.dat")
		echo 3p2 $ip $CHI3 $ANS3 $CHI5 $ANS5 $TOTAL  >> fractiondata/pot3p2/analysed_3p2_dis.dat
	done




cd analysis_scripts
