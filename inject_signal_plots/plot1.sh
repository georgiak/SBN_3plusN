#!/bin/bash
cd ..
mkdir -p inject_signal_plots/plot1data_16


#./sbnfit2 --inject 0:0 --anti --num 2 --plotmode 1 --app --phi 1.6 > inject_signal_plots/plot1data_16/plot1.app.bar.dat &
#./sbnfit2 --inject 0.3:0 --num 2 --app --plotmode 1 --phi 1.6 > inject_signal_plots/plot1data_16/plot1.app.double.dat &
#./sbnfit2 --inject 0:0 --num 2 --app --plotmode 1 --phi 1.6 > inject_signal_plots/plot1data_16/plot1.app.dat &

#./sbnfit2 --inject 0:0 --anti --num 2 --plotmode 1 --both --phi 1.6 > inject_signal_plots/plot1data_16/plot1.both.bar.dat &
#./sbnfit2 --inject 0.3:0 --num 2 --both --plotmode 1 --phi 1.6 > inject_signal_plots/plot1data_16/plot1.both.double.dat &
#./sbnfit2 --inject 0:0 --num 2 --both --plotmode 1 --phi 1.6 > inject_signal_plots/plot1data_16/plot1.both.dat &

./sbnfit2 --inject 0:0 --anti --num 2 --plotmode 1 --app --phi 1.6 --margin > inject_signal_plots/plot1data_16/plot1.app.bar.margin.dat &
./sbnfit2 --inject 0.3:0 --num 2 --app --plotmode 1 --phi 1.6 --margin > inject_signal_plots/plot1data_16/plot1.app.double.margin.dat &
./sbnfit2 --inject 0:0 --num 2 --app --plotmode 1 --phi 1.6 --margin > inject_signal_plots/plot1data_16/plot1.app.margin.dat &

./sbnfit2 --inject 0:0 --anti --num 2 --plotmode 1 --both --phi 1.6 --margin > inject_signal_plots/plot1data_16/plot1.both.bar.margin.dat &
./sbnfit2 --inject 0.3:0 --num 2 --both --plotmode 1 --phi 1.6 --margin > inject_signal_plots/plot1data_16/plot1.both.double.margin.dat &
./sbnfit2 --inject 0:0 --num 2 --both --plotmode 1 --phi 1.6 --margin > inject_signal_plots/plot1data_16/plot1.both.margin.dat &



cd inject_signal_plots
