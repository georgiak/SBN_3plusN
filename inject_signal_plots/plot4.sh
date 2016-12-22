#!/bin/bash
cd ..
mkdir -p inject_signal_plots/plot4data


parallel  --jobs 16  'UE=$(bc <<< "scale=10;sqrt({1})"); ./sbnfit2 --inject 0:0  --plotmode 4 --num 2 --app --phi 4.7 --margin  --signal  0.3981707:1.0:0:$UE:$UE:0:1:1:0:4.7:0:0 > inject_signal_plots/plot4data/plot4_ue4um4_{1}.app.dat ' ::: `cat inject_signal_plots/logs`



parallel  --jobs 16  'UE=$(bc <<< "scale=10;sqrt({1})"); ./sbnfit2 --inject 0:0  --plotmode 4 --num 2 --app --phi 4.7 --anti  --margin  --signal  0.3981707:1.0:0:$UE:$UE:0:1:1:0:4.7:0:0 > inject_signal_plots/plot4data/plot4_ue4um4_{1}.app.bar.dat ' ::: `cat inject_signal_plots/logs`


cd inject_signal_plots
