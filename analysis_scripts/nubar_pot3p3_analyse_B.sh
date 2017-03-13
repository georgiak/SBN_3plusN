#!/bin/bash
cd ..

mkdir -p fractiondata/NUBAR_MODE/pot3p3
#rm fractiondata/NUBAR_MODE/pot3p3/analysed_3p3_dis.bar.dat
rm fractiondata/NUBAR_MODE/pot3p3/analysed_3p3_app.bar.dat


CHI3=30.0973
CHI5=52.1917

for ip in `seq -1.00 0.25 0.50`
do
for ipb in `seq -1.00 0.25 0.50`
do
	TOTAL=$(wc -l < fractiondata/NUBAR_MODE/pot3p3/"3p3_pot_"$ip"_"$ipb"_app.bar.dat" )
	echo "On App POT"$ip
		ANS3=$(awk -v mychi="$CHI3" -v count=0 '{ if ($17 >= mychi && NF <=18) count++}  END {print count}' fractiondata/NUBAR_MODE/pot3p3/"3p3_pot_"$ip"_"$ipb"_app.bar.dat")
		ANS5=$(awk -v mychi="$CHI5" -v count=0 '{ if ($17 >= mychi && NF <=18) count++}  END {print count}' fractiondata/NUBAR_MODE/pot3p3/"3p3_pot_"$ip"_"$ipb"_app.bar.dat")
		echo 3p3 $ip $ipb $CHI3 $ANS3 $CHI5 $ANS5 $TOTAL  >> fractiondata/NUBAR_MODE/pot3p3/analysed_3p3_app.bar.dat

done
done

cd analysis_scripts
