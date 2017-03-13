#!/bin/bash
cd ..

./sbnfit2 --fraction 1 --num 2 --dis >> fractiondata/3p2_dis.dat1 &
./sbnfit2 --fraction 2 --num 2 --dis >> fractiondata/3p2_dis.dat2 &
./sbnfit2 --fraction 3 --num 2 --dis >> fractiondata/3p2_dis.dat3 &
./sbnfit2 --fraction 4 --num 2 --dis >> fractiondata/3p2_dis.dat4 &
./sbnfit2 --fraction 5 --num 2 --dis >> fractiondata/3p2_dis.dat5 &
./sbnfit2 --fraction 6 --num 2 --dis >> fractiondata/3p2_dis.dat6 &
./sbnfit2 --fraction 7 --num 2 --dis >> fractiondata/3p2_dis.dat7 &
./sbnfit2 --fraction 8 --num 2 --dis >> fractiondata/3p2_dis.dat8 &

./sbnfit2 --fraction 1 --num 2 --both >> fractiondata/3p2_both.dat1 &
./sbnfit2 --fraction 2 --num 2 --both >> fractiondata/3p2_both.dat2 &
./sbnfit2 --fraction 3 --num 2 --both >> fractiondata/3p2_both.dat3 &
./sbnfit2 --fraction 4 --num 2 --both >> fractiondata/3p2_both.dat4 &
./sbnfit2 --fraction 5 --num 2 --both >> fractiondata/3p2_both.dat5 &
./sbnfit2 --fraction 6 --num 2 --both >> fractiondata/3p2_both.dat6 &
./sbnfit2 --fraction 7 --num 2 --both >> fractiondata/3p2_both.dat7 &
./sbnfit2 --fraction 8 --num 2 --both >> fractiondata/3p2_both.dat8 &

./sbnfit2 --fraction 1 --num 2 --app >> fractiondata/3p2_app.dat1 &
./sbnfit2 --fraction 2 --num 2 --app >> fractiondata/3p2_app.dat2 &
./sbnfit2 --fraction 3 --num 2 --app >> fractiondata/3p2_app.dat3 &
./sbnfit2 --fraction 4 --num 2 --app >> fractiondata/3p2_app.dat4 &
./sbnfit2 --fraction 5 --num 2 --app >> fractiondata/3p2_app.dat5 &
./sbnfit2 --fraction 6 --num 2 --app >> fractiondata/3p2_app.dat6 &
./sbnfit2 --fraction 7 --num 2 --app >> fractiondata/3p2_app.dat7 &
./sbnfit2 --fraction 8 --num 2 --app >> fractiondata/3p2_app.dat8 &

./sbnfit2 --fraction 1 --num 2 --wierd >> fractiondata/3p2_wierd.dat1 &
./sbnfit2 --fraction 2 --num 2 --wierd >> fractiondata/3p2_wierd.dat2 &
./sbnfit2 --fraction 3 --num 2 --wierd >> fractiondata/3p2_wierd.dat3 &
./sbnfit2 --fraction 4 --num 2 --wierd >> fractiondata/3p2_wierd.dat4 &
./sbnfit2 --fraction 5 --num 2 --wierd >> fractiondata/3p2_wierd.dat5 &
./sbnfit2 --fraction 6 --num 2 --wierd >> fractiondata/3p2_wierd.dat6 &
./sbnfit2 --fraction 7 --num 2 --wierd >> fractiondata/3p2_wierd.dat7 &
./sbnfit2 --fraction 8 --num 2 --wierd >> fractiondata/3p2_wierd.dat8 &


cd generation_scripts
