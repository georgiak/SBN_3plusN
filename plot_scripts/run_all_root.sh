#!/bin/bash

echo "Begining Intrinsic Plotting"
echo '.x plot_sbnd_intrinsic.c++' | root -b -l
echo '.x plot_uboone_intrinsic.c++' | root -b -l
echo '.x plot_icarus_intrinsic.c++' | root -b -l

echo "Begining FullOsc electron App Plotting"

echo "Begining Muon Dis Plotting"

echo "Cleaning up"
rm *.d
rm *.so
echo done


