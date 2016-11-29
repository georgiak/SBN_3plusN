#!/bin/bash
cd ..



parallel  --jobs 16 --joblog nu3p3both.log2   './sbnfit --pot {1} --num 3 --both --anti >> fractiondata/pot3p3/3p3_pot_{1}_{2}_both.bar.dat ' ::: `seq -2.00 0.25 0.50`
parallel --jobs 16 --joblog nu3p3app.log2   './sbnfit --pot {1} --num 3 --app --anti >> fractiondata/pot3p3/3p3_pot_{1}_{2}_app.bar.dat ' ::: `seq -2.00 0.25 0.50`
parallel --jobs 16 --joblog nu3p3app.log2   './sbnfit --pot {1} --num 3 --app --anti >> fractiondata/pot3p3/3p3_pot_{1}_app.bar.dat ' ::: `seq -1.00 0.25 0.50`


cd analysis_scripts
