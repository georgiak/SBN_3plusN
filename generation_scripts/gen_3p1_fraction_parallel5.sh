#!/bin/bash
cd ..
#rm fractiondata/3p1*

./sbnfit2 --fraction 0 --num 1 --dis  >> fractiondata/3p1_dis.dat.new &

#./sbnfit2 --fraction 0 --num 1 --both  >> fractiondata/3p1_both.dat.new &

./sbnfit2 --fraction 0 --num 1 --app  >> fractiondata/3p1_app.dat.new &

./sbnfit2 --fraction 0 --num 1 --wierd >> fractiondata/3p1_wierd.dat.new &

cd generation_scripts
