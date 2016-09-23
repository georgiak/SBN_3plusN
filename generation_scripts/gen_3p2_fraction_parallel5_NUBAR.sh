#!/bin/bash
cd ..

./sbnfit --fraction 1 --num 2 --dis --anti >> fractiondata/NUBAR_MODE/3p2_dis.bar.dat.1 &
./sbnfit --fraction 2 --num 2 --dis --anti >> fractiondata/NUBAR_MODE/3p2_dis.bar.dat.2 &
./sbnfit --fraction 3 --num 2 --dis --anti >> fractiondata/NUBAR_MODE/3p2_dis.bar.dat.3 &
./sbnfit --fraction 4 --num 2 --dis --anti >> fractiondata/NUBAR_MODE/3p2_dis.bar.dat.4 &
./sbnfit --fraction 5 --num 2 --dis --anti >> fractiondata/NUBAR_MODE/3p2_dis.bar.dat.5 &


./sbnfit --fraction 1 --num 2 --both --anti >> fractiondata/NUBAR_MODE/3p2_both.bar.dat.1 &
./sbnfit --fraction 2 --num 2 --both --anti >> fractiondata/NUBAR_MODE/3p2_both.bar.dat.2 &
./sbnfit --fraction 3 --num 2 --both --anti >> fractiondata/NUBAR_MODE/3p2_both.bar.dat.3 &
./sbnfit --fraction 4 --num 2 --both --anti >> fractiondata/NUBAR_MODE/3p2_both.bar.dat.4 & 
./sbnfit --fraction 5 --num 2 --both --anti >> fractiondata/NUBAR_MODE/3p2_both.bar.dat.5 &

./sbnfit --fraction 1 --num 2 --app --anti >> fractiondata/NUBAR_MODE/3p2_app.bar.dat.1 &
./sbnfit --fraction 2 --num 2 --app --anti >> fractiondata/NUBAR_MODE/3p2_app.bar.dat.2 &
./sbnfit --fraction 3 --num 2 --app --anti >> fractiondata/NUBAR_MODE/3p2_app.bar.dat.3 &
./sbnfit --fraction 4 --num 2 --app --anti >> fractiondata/NUBAR_MODE/3p2_app.bar.dat.4 &
./sbnfit --fraction 5 --num 2 --app --anti >> fractiondata/NUBAR_MODE/3p2_app.bar.dat.5 &

cd generation_scripts
