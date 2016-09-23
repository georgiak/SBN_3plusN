#!/bin/bash
cd ..

./sbnfit --fraction 0 --num 1 --dis --anti >> fractiondata/NUBAR_MODE/3p1_dis.bar.dat &

./sbnfit --fraction 0 --num 1 --both --anti >> fractiondata/NUBAR_MODE/3p1_both.bar.dat &

./sbnfit --fraction 0 --num 1 --app --anti >> fractiondata/NUBAR_MODE/3p1_app.bar.dat &

cd generation_scripts
