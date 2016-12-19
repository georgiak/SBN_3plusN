#!/bin/bash
cd ..
mkdir -p fractiondata/NUBAR_MODE/pot3p3
rm fractiondata/NUBAR_MODE/pot3p3/3p3_pot_*


parallel  --jobs 16 --joblog nu3p3both.log   './sbnfit2 --pot {1}:{2} --num 3 --both --anti >> fractiondata/NUBAR_MODE/pot3p3/3p3_pot_{1}_{2}_both.bar.dat ' ::: `seq -2.00 0.25 0.50` ::: `seq -2.00 0.25 0.50` 

#parallel --jobs 16 --joblog nu3p3app.log   './sbnfit --pot {1}:{2} --num 3 --app --anti >> fractiondata/NUBAR_MODE/pot3p3/3p3_pot_{1}_{2}_app.bar.dat ' ::: `seq -1.00 0.25 0.50` ::: `seq -1.00 0.25 0.50` 


cd analysis_scripts
