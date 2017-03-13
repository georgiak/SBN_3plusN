#!/bin/bash
cd ..
mkdir -p inject_signal_plots/plot3data

# first one using BF values
parallel  --jobs 16  './sbnfit2 --inject 0:0 --num 2 --both --anti --phi {1} --plotmode 3 > inject_signal_plots/plot3data/plot3_phi_{1}.both.bf.bar.dat ' ::: `seq 0.0 0.1 6.2`
parallel  --jobs 16  './sbnfit2 --inject 0:0 --num 2 --both --anti --phi {1} --plotmode 3 > inject_signal_plots/plot3data/plot3_phi_{1}.both.bf.bar.dat ' ::: `seq 0.0 0.1 6.2`

# Next two using a large max value (not very large tho)
# 0.3981707:1.0:0:0.17:0.17:0:0.2:0.2:0:{1}:0:0 
parallel  --jobs 16  './sbnfit2 --inject 0:0 --plotmode 3 --num 2 --both --anti --phi {1}  --signal  0.3981707:1.0:0:0.17:0.17:0:0.2:0.2:0:{1}:0:0  > inject_signal_plots/plot3data/plot3_phi_{1}.both.max.bar.dat ' ::: `seq 0.0 0.1 6.2`
parallel  --jobs 16  './sbnfit2 --inject 0:0 --plotmode 3 --num 2 --both --anti --phi {1}  --signal  0.3981707:1.0:0:0.17:0.17:0:0.2:0.2:0:{1}:0:0 > inject_signal_plots/plot3data/plot3_phi_{1}.both.max.bar.dat ' ::: `seq 0.0 0.1 6.2`



# first one using BF values
parallel  --jobs 16  './sbnfit2 --inject 0:0 --num 2 --both --plotmode 3 --phi {1} --margin --anti > inject_signal_plots/plot3data/plot3_phi_{1}.both.bf.margin.bar.dat ' ::: `seq 0.0 0.1 6.2`
parallel  --jobs 16  './sbnfit2 --inject 0:0 --num 2 --both --plotmode 3 --phi {1} --margin --anti > inject_signal_plots/plot3data/plot3_phi_{1}.both.bf.margin.bar.dat ' ::: `seq 0.0 0.1 6.2`

# Next two using a large max value (not very large tho)
# 0.3981707:1.0:0:0.17:0.17:0:0.2:0.2:0:{1}:0:0 
parallel  --jobs 16  './sbnfit2 --inject 0:0 --num 2 --both --plotmode 3 --phi {1} --margin --anti --signal  0.3981707:1.0:0:0.17:0.17:0:0.2:0.2:0:{1}:0:0  > inject_signal_plots/plot3data/plot3_phi_{1}.both.max.margin.bar.dat ' ::: `seq 0.0 0.1 6.2`
parallel  --jobs 16  './sbnfit2 --inject 0:0 --num 2 --both --plotmode 3 --phi {1} --margin --anti --signal  0.3981707:1.0:0:0.17:0.17:0:0.2:0.2:0:{1}:0:0 > inject_signal_plots/plot3data/plot3_phi_{1}.both.max.margin.bar.dat ' ::: `seq 0.0 0.1 6.2`


cd inject_signal_plots
