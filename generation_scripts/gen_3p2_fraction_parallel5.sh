#!/bin/bash
cd ..

./sbnfit --fraction 1 --num 2 --dis >> fractiondata/3p2_dis.dat &
./sbnfit --fraction 2 --num 2 --dis >> fractiondata/3p2_dis.dat &
./sbnfit --fraction 3 --num 2 --dis >> fractiondata/3p2_dis.dat &
./sbnfit --fraction 4 --num 2 --dis >> fractiondata/3p2_dis.dat &
./sbnfit --fraction 5 --num 2 --dis >> fractiondata/3p2_dis.dat &


./sbnfit --fraction 1 --num 2 --both >> fractiondata/3p2_both.dat &
./sbnfit --fraction 2 --num 2 --both >> fractiondata/3p2_both.dat &
./sbnfit --fraction 3 --num 2 --both >> fractiondata/3p2_both.dat &
./sbnfit --fraction 4 --num 2 --both >> fractiondata/3p2_both.dat &
./sbnfit --fraction 5 --num 2 --both >> fractiondata/3p2_both.dat &

./sbnfit --fraction 1 --num 2 --app >> fractiondata/3p2_app.dat &
./sbnfit --fraction 2 --num 2 --app >> fractiondata/3p2_app.dat &
./sbnfit --fraction 3 --num 2 --app >> fractiondata/3p2_app.dat &
./sbnfit --fraction 4 --num 2 --app >> fractiondata/3p2_app.dat &
./sbnfit --fraction 5 --num 2 --app >> fractiondata/3p2_app.dat &

cd generation_scripts
