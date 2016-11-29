#!/bin/bash
cd ..
mkdir -p fractiondata/NUBAR_MODE/pot3p2
rm fractiondata/NUBAR_MODE/pot3p2/3p2_pot_*



parallel  --jobs 16 --joblog nu3p2both.log   './sbnfit --pot {1}:{2} --num 2 --both --anti >> fractiondata/NUBAR_MODE/pot3p2/3p2_pot_{1}_{2}_both.bar.dat ' ::: `seq -2.00 0.25 0.50` ::: `seq -2.00 0.25 0.50` 

parallel --jobs 16 --joblog nu3p2app.log   './sbnfit --pot {1}:{2} --num 2 --app --anti >> fractiondata/NUBAR_MODE/pot3p2/3p2_pot_{1}_{2}_app.bar.dat ' ::: `seq -2.00 0.25 0.50` ::: `seq -2.00 0.25 0.50` 


cd analysis_scripts
