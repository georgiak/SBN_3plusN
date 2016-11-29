#!/bin/bash
cd ..

mkdir -p fractiondata/NUBAR_MODE/pot3p2
rm fractiondata/NUBAR_MODE/pot3p2/analysed_3p2_both.bar.dat
#rm fractiondata/NUBAR_MODE/pot3p2/analysed_3p2_dis.bar.dat


CHI3=21.8466
CHI5=41.7798


for ip in `seq -2.00 0.25 0.50`
do
for ipb in `seq -2.00 0.25 0.50`
do
	TOTAL=$(wc -l < fractiondata/NUBAR_MODE/pot3p2/"3p2_pot_"$ip"_"$ipb"_both.bar.dat" )
	echo "On Both POT"$ip
		ANS3=$(awk -v mychi="$CHI3" -v count=0 '{ if ($17 >= mychi && NF <=18 ) count++}  END {print count}' fractiondata/NUBAR_MODE/pot3p2/"3p2_pot_"$ip"_"$ipb"_both.bar.dat")
		ANS5=$(awk -v mychi="$CHI5" -v count=0 '{ if ($17 >= mychi && NF <=18) count++}  END {print count}' fractiondata/NUBAR_MODE/pot3p2/"3p2_pot_"$ip"_"$ipb"_both.bar.dat")
		echo 3p2 $ip $ipb $CHI3 $ANS3 $CHI5 $ANS5 $TOTAL  >> fractiondata/NUBAR_MODE/pot3p2/analysed_3p2_both.bar.dat

done
done

cd analysis_scripts
