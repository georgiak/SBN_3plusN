#!/bin/bash
cd ..
mkdir -p fractiondata/NUBAR_MODE/pot3p1
rm fractiondata/NUBAR_MODE/pot3p1/3p1_pot_*



parallel  --jobs 10 --joblog nu3p1both.log   './sbnfit2 --pot {1}:{2} --num 1 --both --anti >> fractiondata/NUBAR_MODE/pot3p1/3p1_pot_{1}_{2}_both.bar.dat ' ::: `seq -2.00 0.25 0.50` ::: `seq -2.00 0.25 0.50` 

parallel --jobs 10 --joblog nu3p1app.log   './sbnfit2 --pot {1}:{2} --num 1 --app --anti >> fractiondata/NUBAR_MODE/pot3p1/3p1_pot_{1}_{2}_app.bar.dat ' ::: `seq -2.00 0.25 0.50` ::: `seq -2.00 0.25 0.50` 


cd analysis_scripts
