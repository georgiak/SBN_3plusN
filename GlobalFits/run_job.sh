#!/bin/bash
set -x
echo Start  `date`
echo Site:${GLIDEIN_ResourceName}
echo "the worker node is " `hostname` "OS: "  `uname -a`
whoami
id

cd $_CONDOR_SCRATCH_DIR

IFDH_OPTION=""

if [ -z $GROUP ]; then

# try to figure out what group the user is in
GROUP=`id -gn`

fi

export IFDH_DEBUG=1

case $GROUP in

lar1nd)
SCRATCH_DIR="/pnfs/lar1nd/scratch/users"
;;
esac

voms-proxy-info --all


source /grid/fermiapp/products/common/etc/setups.sh
#source /cvmfs/oasis.opensciencegrid.org/fermilab/products/common/etc/setups.sh
#source /cvmfs/oasis.opensciencegrid.org/fermilab/products/larsoft/setup
setup ifdhc

echo "Here is the your environment in this job: " > job_output_${CLUSTER}.${PROCESS}.log
env >> job_output_${CLUSTER}.${PROCESS}.log

echo "group = $GROUP"

if [ -z ${GRID_USER} ]; then
GRID_USER=`basename $X509_USER_PROXY | cut -d "_" -f 2`
fi

echo "GRID_USER = `echo $GRID_USER`"

umask 002

cd $_CONDOR_SCRATCH_DIR


####################################
###### setup your needed products here, e.g. geant4 etc...
####################################

source /grid/fermiapp/products/uboone/setup_uboone.sh; setup uboonecode v04_25_00 -qdebug:e7

####################################
#### This is where you copy all of your executable/necessary files to the worker node
#### ( If applicable )
####################################

ifdh cp -D /pnfs/lar1nd/scratch/users/dcianci/fits/globalFit .
ifdh cp -D /pnfs/lar1nd/scratch/users/dcianci/fits/jobOptions.txt .
ifdh cp -D /pnfs/lar1nd/scratch/users/dcianci/fits/data/miniboone_binboundaries.txt .
ifdh cp -D /pnfs/lar1nd/scratch/users/dcianci/fits/data/miniboone_binboundaries_disap.txt .
ifdh cp -D /pnfs/lar1nd/scratch/users/dcianci/fits/data/miniboone_binboundaries_lowe.txt .
ifdh cp -D /pnfs/lar1nd/scratch/users/dcianci/fits/data/miniboone_binboundaries_nubar.txt .
ifdh cp -D /pnfs/lar1nd/scratch/users/dcianci/fits/data/miniboone_binboundaries_nubar_lowe.txt .
ifdh cp -D /pnfs/lar1nd/scratch/users/dcianci/fits/data/miniboone_binboundaries_numu.txt .
ifdh cp -D /pnfs/lar1nd/scratch/users/dcianci/fits/data/miniboone_frac_shape_matrix_numu_disap.txt .
ifdh cp -D /pnfs/lar1nd/scratch/users/dcianci/fits/data/miniboone_frac_shape_matrix_numubar_disap.txt .
ifdh cp -D /pnfs/lar1nd/scratch/users/dcianci/fits/data/miniboone_full_fractcovmatrix.txt .
ifdh cp -D /pnfs/lar1nd/scratch/users/dcianci/fits/data/miniboone_full_fractcovmatrix_lowe.txt .
ifdh cp -D /pnfs/lar1nd/scratch/users/dcianci/fits/data/miniboone_full_fractcovmatrix_nubar.txt .
ifdh cp -D /pnfs/lar1nd/scratch/users/dcianci/fits/data/miniboone_full_fractcovmatrix_nubar_lowe.txt .
ifdh cp -D /pnfs/lar1nd/scratch/users/dcianci/fits/data/miniboone_full_fractcovmatrix_nunubar.txt .
ifdh cp -D /pnfs/lar1nd/scratch/users/dcianci/fits/data/miniboone_full_fractcovmatrix_nunubar_lowe.txt .
ifdh cp -D /pnfs/lar1nd/scratch/users/dcianci/fits/data/miniboone_nuebarbgr.txt .
ifdh cp -D /pnfs/lar1nd/scratch/users/dcianci/fits/data/miniboone_nuebarbgr_lowe.txt .
ifdh cp -D /pnfs/lar1nd/scratch/users/dcianci/fits/data/miniboone_nuebardata.txt .
ifdh cp -D /pnfs/lar1nd/scratch/users/dcianci/fits/data/miniboone_nuebardata_lowe.txt .
ifdh cp -D /pnfs/lar1nd/scratch/users/dcianci/fits/data/miniboone_nuebgr.txt .
ifdh cp -D /pnfs/lar1nd/scratch/users/dcianci/fits/data/miniboone_nuebgr_lowe.txt .
ifdh cp -D /pnfs/lar1nd/scratch/users/dcianci/fits/data/miniboone_nuedata.txt .
ifdh cp -D /pnfs/lar1nd/scratch/users/dcianci/fits/data/miniboone_nuedata_lowe.txt .
ifdh cp -D /pnfs/lar1nd/scratch/users/dcianci/fits/data/miniboone_numode_fullosc_ntuple.txt .
ifdh cp -D /pnfs/lar1nd/scratch/users/dcianci/fits/data/miniboone_numu.txt .
ifdh cp -D /pnfs/lar1nd/scratch/users/dcianci/fits/data/miniboone_numubar.txt .
ifdh cp -D /pnfs/lar1nd/scratch/users/dcianci/fits/data/miniboone_numubardata.txt .
ifdh cp -D /pnfs/lar1nd/scratch/users/dcianci/fits/data/miniboone_numubardata_disap.txt .
ifdh cp -D /pnfs/lar1nd/scratch/users/dcianci/fits/data/miniboone_numubarnuebarfullosc_ntuple.txt .
ifdh cp -D /pnfs/lar1nd/scratch/users/dcianci/fits/data/miniboone_numudata.txt .
ifdh cp -D /pnfs/lar1nd/scratch/users/dcianci/fits/data/miniboone_numudata_disap.txt .
ifdh cp -D /pnfs/lar1nd/scratch/users/dcianci/fits/data/miniboone_numunuefullosc_ntuple.txt .
ifdh cp -D /pnfs/lar1nd/scratch/users/dcianci/fits/data/miniboone_nunubarfullosc_ntuple.txt .
ifdh cp -D /pnfs/lar1nd/scratch/users/dcianci/fits/data/neutrino_frac_error_matrix.txt .
ifdh cp -D /pnfs/lar1nd/scratch/users/dcianci/fits/data/numi_fullosc.out .
ifdh cp -D /pnfs/lar1nd/scratch/users/dcianci/fits/data/numi_fullosc.txt .
ifdh cp -D /pnfs/lar1nd/scratch/users/dcianci/fits/data/numubardisap_ntuple.txt .
ifdh cp -D /pnfs/lar1nd/scratch/users/dcianci/fits/data/numudisap_ntuple.txt .

#######
####### launch executable
#######

chmod u+x globalFit
chmod u+x jobOptions.txt
chmod u+x miniboone_binboundaries.txt
chmod u+x miniboone_binboundaries_disap.txt
chmod u+x miniboone_binboundaries_lowe.txt
chmod u+x miniboone_binboundaries_nubar.txt
chmod u+x miniboone_binboundaries_nubar_lowe.txt
chmod u+x miniboone_binboundaries_numu.txt
chmod u+x miniboone_frac_shape_matrix_numu_disap.txt
chmod u+x miniboone_frac_shape_matrix_numubar_disap.txt
chmod u+x miniboone_full_fractcovmatrix.txt
chmod u+x miniboone_full_fractcovmatrix_lowe.txt
chmod u+x miniboone_full_fractcovmatrix_nubar.txt
chmod u+x miniboone_full_fractcovmatrix_nubar_lowe.txt
chmod u+x miniboone_full_fractcovmatrix_nunubar.txt
chmod u+x miniboone_full_fractcovmatrix_nunubar_lowe.txt
chmod u+x miniboone_nuebarbgr.txt
chmod u+x miniboone_nuebarbgr_lowe.txt
chmod u+x miniboone_nuebardata.txt
chmod u+x miniboone_nuebardata_lowe.txt
chmod u+x miniboone_nuebgr.txt
chmod u+x miniboone_nuebgr_lowe.txt
chmod u+x miniboone_nuedata.txt
chmod u+x miniboone_nuedata_lowe.txt
chmod u+x miniboone_numode_fullosc_ntuple.txt
chmod u+x miniboone_numu.txt
chmod u+x miniboone_numubar.txt
chmod u+x miniboone_numubardata.txt
chmod u+x miniboone_numubardata_disap.txt
chmod u+x miniboone_numubarnuebarfullosc_ntuple.txt
chmod u+x miniboone_numudata.txt
chmod u+x miniboone_numudata_disap.txt
chmod u+x miniboone_numunuefullosc_ntuple.txt
chmod u+x miniboone_nunubarfullosc_ntuple.txt
chmod u+x neutrino_frac_error_matrix.txt
chmod u+x numi_fullosc.out
chmod u+x numi_fullosc.txt
chmod u+x numubardisap_ntuple.txt
chmod u+x numudisap_ntuple.txt

./globalFit
chmod u+x globFit.root

#######
####### Copy results back
#######

ifdh mkdir /lar1nd/data/users/dcianci/output_all

ifdh cp globFit.root /lar1nd/data/users/dcianci/output_all/globFit_${PROCESS}.root

exit 0
