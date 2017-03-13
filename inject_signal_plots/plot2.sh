#!/bin/bash
cd ..
mkdir -p inject_signal_plots/plot2data


./sbnfit2 --inject 0:0 --anti --num 2 --plotmode 2 --app  > inject_signal_plots/plot2data/plot2.app.bar.dat &
./sbnfit2 --inject 0.3:0 --num 2 --app --plotmode 2 > inject_signal_plots/plot2data/plot2.app.dat &


cd inject_signal_plots
