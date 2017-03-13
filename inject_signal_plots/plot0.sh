#!/bin/bash
cd ..
mkdir -p inject_signal_plots/plot0data


parallel  --jobs 4  './sbnfit2 --inject 0:0 --plotmode 1 --num 2 --app --phi {1} > inject_signal_plots/plot0data/plot0_phi_{1}.app.dat ' ::: `seq 0.0 0.1 6.2`

cd inject_signal_plots
