#!/bin/bash
cd ..
mkdir -p inject_signal_plots/plot4data/m4_0.69_m5_0.9_phi_1.57/
mkdir -p inject_signal_plots/plot4data/m4_0.69_m5_0.9_phi_4.7/

#parallel  --jobs 16  'UE=$(bc <<< "scale=10;sqrt({1})"); ./sbnfit2 --inject 0:0  --plotmode 4 --num 2 --app --phi 4.7 --margin  --signal  0.3981707:1.0:0:$UE:$UE:0:1:1:0:4.7:0:0 > inject_signal_plots/plot4data/plot4_ue4um4_{1}.app.dat ' ::: `cat inject_signal_plots/logs`



#parallel  --jobs 16  'UE=$(bc <<< "scale=10;sqrt({1})"); ./sbnfit2 --inject 0:0  --plotmode 4 --num 2 --app --phi 4.7 --anti  --margin  --signal  0.3981707:1.0:0:$UE:$UE:0:1:1:0:4.7:0:0 > inject_signal_plots/plot4data/plot4_ue4um4_{1}.app.bar.dat ' ::: `cat inject_signal_plots/logs`

#one at 4.7

#with 
parallel  --jobs 16  'UE=$(bc <<< "scale=10;sqrt({1})"); ./sbnfit2 --inject 0.3:0  --plotmode 4 --num 2 --app --phi 1.57    --margin  --signal  0.691831:0.912011:0:$UE:$UE:0:1:1:0:1.57:0:0  > inject_signal_plots/plot4data/m4_0.69_m5_0.9_phi_1.57/plot4_true_1.57_ue4um4_{1}.app.dat ' ::: `cat inject_signal_plots/logs2`

parallel  --jobs 16  'UE=$(bc <<< "scale=10;sqrt({1})"); ./sbnfit2 --inject 0.3:0  --plotmode 4 --num 2 --app --phi 4.712    --margin  --signal  0.691831:0.912011:0:$UE:$UE:0:1:1:0:4.712:0:0  > inject_signal_plots/plot4data/m4_0.69_m5_0.9_phi_4.7/plot4_true_4.7_ue4um4_{1}.app.dat ' ::: `cat inject_signal_plots/logs2`

# and with anti
parallel  --jobs 16  'UE=$(bc <<< "scale=10;sqrt({1})"); ./sbnfit2 --inject 0:0  --plotmode 4 --num 2 --app --phi 1.57 --anti  --margin  --signal  0.691831:0.912011:0:$UE:$UE:0:1:1:0:1.57:0:0  > inject_signal_plots/plot4data/m4_0.69_m5_0.9_phi_1.57/plot4_true_1.57_ue4um4_{1}.app.bar.dat ' ::: `cat inject_signal_plots/logs2`

parallel  --jobs 16  'UE=$(bc <<< "scale=10;sqrt({1})"); ./sbnfit2 --inject 0:0  --plotmode 4 --num 2 --app --phi 4.712 --anti  --margin  --signal  0.691831:0.912011:0:$UE:$UE:0:1:1:0:4.712:0:0  > inject_signal_plots/plot4data/m4_0.69_m5_0.9_phi_4.7/plot4_true_4.7_ue4um4_{1}.app.bar.dat ' ::: `cat inject_signal_plots/logs2`

# one at 4.7
#parallel  --jobs 16  'UE=$(bc <<< "scale=10;sqrt({1})"); ./sbnfit2 --inject 0:0  --plotmode 4 --num 2 --app --phi 4.7 --anti  --margin  --signal  0.7:1.0:0:$UE:$UE:0:1:1:0:1.6:0:0 > inject_signal_plots/plot4data/m4_0.7_m5_1_phi_1.6/plot4_ue4um4_{1}.app.bar.dat ' ::: `cat inject_signal_plots/logs2`

cd inject_signal_plots
