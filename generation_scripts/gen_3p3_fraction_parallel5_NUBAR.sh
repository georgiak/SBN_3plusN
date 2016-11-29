#!/bin/bash
cd ..

#./sbnfit --fraction 1 --num 3 --dis --anti >> fractiondata/NUBAR_MODE/3p3_dis.bar.dat.1 &
#./sbnfit --fraction 2 --num 3 --dis --anti >> fractiondata/NUBAR_MODE/3p3_dis.bar.dat.2 &
#./sbnfit --fraction 3 --num 3 --dis --anti >> fractiondata/NUBAR_MODE/3p3_dis.bar.dat.3 &
#./sbnfit --fraction 4 --num 3 --dis --anti >> fractiondata/NUBAR_MODE/3p3_dis.bar.dat.4 &
#./sbnfit --fraction 5 --num 3 --dis --anti >> fractiondata/NUBAR_MODE/3p3_dis.bar.dat.5 &

./sbnfit --fraction 1 --num 3 --both --anti >> fractiondata/NUBAR_MODE/3p3_both.bar.dat.1 &
./sbnfit --fraction 2 --num 3 --both --anti >> fractiondata/NUBAR_MODE/3p3_both.bar.dat.2 &
./sbnfit --fraction 3 --num 3 --both --anti >> fractiondata/NUBAR_MODE/3p3_both.bar.dat.3 &
./sbnfit --fraction 4 --num 3 --both --anti >> fractiondata/NUBAR_MODE/3p3_both.bar.dat.4 & 
./sbnfit --fraction 5 --num 3 --both --anti >> fractiondata/NUBAR_MODE/3p3_both.bar.dat.5 & 
./sbnfit --fraction 6 --num 3 --both --anti >> fractiondata/NUBAR_MODE/3p3_both.bar.dat.6 &
./sbnfit --fraction 7 --num 3 --both --anti >> fractiondata/NUBAR_MODE/3p3_both.bar.dat.7 &
./sbnfit --fraction 8 --num 3 --both --anti >> fractiondata/NUBAR_MODE/3p3_both.bar.dat.8 &



cd generation_scripts
