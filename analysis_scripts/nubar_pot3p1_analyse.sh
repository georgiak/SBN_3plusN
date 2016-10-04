#!/bin/bash
cd ..

mkdir -p fractiondata/NUBAR_MODE/pot3p1
rm fractiondata/NUBAR_MODE/pot3p1/analysed_3p1_both.bar.dat
rm fractiondata/NUBAR_MODE/pot3p1/analysed_3p1_dis.bar.dat
rm fractiondata/NUBAR_MODE/pot3p1/analysed_3p1_app.bar.dat

CHI3=14.1564
CHI5=31.8121

for ip in `seq -2.00 0.25 0.50`
do
for ipb in `seq -2.00 0.25 0.50`
do
	TOTAL=$(wc -l < fractiondata/NUBAR_MODE/pot3p1/"3p1_pot_"$ip"_"$ipb"_both.bar.dat" )
	echo "On Both POT"$ip
		ANS3=$(awk -v mychi="$CHI3" -v count=0 '{ if ($17 >= mychi) count++}  END {print count}' fractiondata/NUBAR_MODE/pot3p1/"3p1_pot_"$ip"_"$ipb"_both.bar.dat")
		ANS5=$(awk -v mychi="$CHI5" -v count=0 '{ if ($17 >= mychi) count++}  END {print count}' fractiondata/NUBAR_MODE/pot3p1/"3p1_pot_"$ip"_"$ipb"_both.bar.dat")
		echo 3p1 $ip $ipb $CHI3 $ANS3 $CHI5 $ANS5 $TOTAL  >> fractiondata/NUBAR_MODE/pot3p1/analysed_3p1_both.bar.dat

done
done

for ip in `seq -2.00 0.25 0.50`
do
for ipb in `seq -2.00 0.25 0.50`
do
	TOTAL=$(wc -l < fractiondata/NUBAR_MODE/pot3p1/"3p1_pot_"$ip"_"$ipb"_app.bar.dat" )
	echo "On App POT"$ip
		ANS3=$(awk -v mychi="$CHI3" -v count=0 '{ if ($17 >= mychi) count++}  END {print count}' fractiondata/NUBAR_MODE/pot3p1/"3p1_pot_"$ip"_"$ipb"_app.bar.dat")
		ANS5=$(awk -v mychi="$CHI5" -v count=0 '{ if ($17 >= mychi) count++}  END {print count}' fractiondata/NUBAR_MODE/pot3p1/"3p1_pot_"$ip"_"$ipb"_app.bar.dat")
		echo 3p1 $ip $ipb $CHI3 $ANS3 $CHI5 $ANS5 $TOTAL  >> fractiondata/NUBAR_MODE/pot3p1/analysed_3p1_app.bar.dat

done
done

#for ip in `seq -2.00 0.25 0.50`
#do
#for ipb in `seq -2.00 0.25 0.50`
#do

#	TOTAL=$(wc -l < fractiondata/NUBAR_MODE/pot3p1/"3p1_pot_"$ip"_"$ipb"_dis.bar.dat" )
#	echo "On Dis POT"$ip
#		ANS3=$(awk -v mychi="$CHI3" -v count=0 '{ if ($17 >= mychi) count++}  END {print count}' fractiondata/NUBAR_MODE/pot3p1/"3p1_pot_"$ip"_"$ipb"_dis.bar.dat")
#		ANS5=$(awk -v mychi="$CHI5" -v count=0 '{ if ($17 >= mychi) count++}  END {print count}' fractiondata/NUBAR_MODE/pot3p1/"3p1_pot_"$ip"_"$ipb"_dis.bar.dat")
#		echo 3p1 $ip $ipb $CHI3 $ANS3 $CHI5 $ANS5 $TOTAL  >> fractiondata/NUBAR_MODE/pot3p1/analysed_3p1_dis.bar.dat
#done
#done



cd analysis_scripts
