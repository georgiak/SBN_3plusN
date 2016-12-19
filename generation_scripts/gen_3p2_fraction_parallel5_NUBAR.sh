#!/bin/bash
cd ..

#./sbnfit --fraction 1 --num 2 --dis --anti >> fractiondata/NUBAR_MODE/3p2_dis.bar.dat.1 &
#./sbnfit --fraction 2 --num 2 --dis --anti >> fractiondata/NUBAR_MODE/3p2_dis.bar.dat.2 &
#./sbnfit --fraction 3 --num 2 --dis --anti >> fractiondata/NUBAR_MODE/3p2_dis.bar.dat.3 &
#./sbnfit --fraction 4 --num 2 --dis --anti >> fractiondata/NUBAR_MODE/3p2_dis.bar.dat.4 &
#./sbnfit --fraction 5 --num 2 --dis --anti >> fractiondata/NUBAR_MODE/3p2_dis.bar.dat.5 &

./sbnfit2 --fraction 1 --num 2 --both --anti >> fractiondata/NUBAR_MODE/3p2_both.bar.dat.1 &
./sbnfit2 --fraction 2 --num 2 --both --anti >> fractiondata/NUBAR_MODE/3p2_both.bar.dat.2 &
./sbnfit2 --fraction 3 --num 2 --both --anti >> fractiondata/NUBAR_MODE/3p2_both.bar.dat.3 &
./sbnfit2 --fraction 4 --num 2 --both --anti >> fractiondata/NUBAR_MODE/3p2_both.bar.dat.4 & 
./sbnfit2 --fraction 5 --num 2 --both --anti >> fractiondata/NUBAR_MODE/3p2_both.bar.dat.5 & 
./sbnfit2 --fraction 6 --num 2 --both --anti >> fractiondata/NUBAR_MODE/3p2_both.bar.dat.6 &
./sbnfit2 --fraction 7 --num 2 --both --anti >> fractiondata/NUBAR_MODE/3p2_both.bar.dat.7 &
./sbnfit2 --fraction 8 --num 2 --both --anti >> fractiondata/NUBAR_MODE/3p2_both.bar.dat.8 &



cd generation_scripts
