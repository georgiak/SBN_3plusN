cp /lar1nd/app/users/dcianci/SBN_3plusN/GlobalFits/inputs/jobOptionsDiag.txt jobOptions.txt
cp jobOptions.txt jobOptions1.txt

sed -i -e 's/CCFRProcess=0/CCFRProcess=1/g' jobOptions1.txt
mv jobOptions1.txt /pnfs/lar1nd/scratch/users/dcianci/fits/jobOptions.txt
cp run_job_diag.sh run_job_diag1.sh
sed -i -e 's/VMODE=""/VMODE="CCFR"/g' run_job_diag1.sh
jobsub_submit -G lar1nd --role=Analysis -N 100 --resource-provides=usage_model=DEDICATED,OPPORTUNISTIC --OS=SL6 file://`pwd`/run_job_diag.sh
cp jobOptions.txt jobOptions1.txt
cp run_job_diag.sh run_job_diag1.sh

sed -i -e 's/CDHSProcess=0/CDHSProcess=1/g' jobOptions1.txt
mv jobOptions1.txt /pnfs/lar1nd/scratch/users/dcianci/fits/jobOptions.txt
cp run_job_diag.sh run_job_diag1.sh
sed -i -e 's/VMode=""/VMode="CDHSProcess/g' run_job_diag1.sh
jobsub_submit -G lar1nd --role=Analysis -N 100 --resource-provides=usage_model=DEDICATED,OPPORTUNISTIC --OS=SL6 file://`pwd`/run_job_diag.sh
cp jobOptions.txt jobOptions1.txt
cp run_job_diag.sh run_job_diag1.sh

sed -i -e 's/CHOOZProcess=0/CHOOZProcess=1/g' jobOptions1.txt
mv jobOptions1.txt /pnfs/lar1nd/scratch/users/dcianci/fits/jobOptions.txt
cp run_job_diag.sh run_job_diag1.sh
sed -i -e 's/VMode=""/VMode="CHOOZProcess/g' run_job_diag1.sh
jobsub_submit -G lar1nd --role=Analysis -N 100 --resource-provides=usage_model=DEDICATED,OPPORTUNISTIC --OS=SL6 file://`pwd`/run_job_diag.sh
cp jobOptions.txt jobOptions1.txt
cp run_job_diag.sh run_job_diag1.sh

sed -i -e 's/KARMENProcess=0/KARMENProcess=1/g' jobOptions1.txt
mv jobOptions1.txt /pnfs/lar1nd/scratch/users/dcianci/fits/jobOptions.txt
cp run_job_diag.sh run_job_diag1.sh
sed -i -e 's/VMode=""/VMode="KARMENProcess/g' run_job_diag1.sh
jobsub_submit -G lar1nd --role=Analysis -N 100 --resource-provides=usage_model=DEDICATED,OPPORTUNISTIC --OS=SL6 file://`pwd`/run_job_diag.sh
cp jobOptions.txt jobOptions1.txt
cp run_job_diag.sh run_job_diag1.sh

sed -i -e 's/LSNDProcess=0/LSNDProcess=1/g' jobOptions1.txt
mv jobOptions1.txt /pnfs/lar1nd/scratch/users/dcianci/fits/jobOptions.txt
cp run_job_diag.sh run_job_diag1.sh
sed -i -e 's/VMode=""/VMode="LSNDProcess/g' run_job_diag1.sh
jobsub_submit -G lar1nd --role=Analysis -N 100 --resource-provides=usage_model=DEDICATED,OPPORTUNISTIC --OS=SL6 file://`pwd`/run_job_diag.sh
cp jobOptions.txt jobOptions1.txt
cp run_job_diag.sh run_job_diag1.sh

sed -i -e 's/NOMADProcess=0/NOMADProcess=1/g' jobOptions1.txt
mv jobOptions1.txt /pnfs/lar1nd/scratch/users/dcianci/fits/jobOptions.txt
cp run_job_diag.sh run_job_diag1.sh
sed -i -e 's/VMode=""/VMode="NOMADProcess/g' run_job_diag1.sh
jobsub_submit -G lar1nd --role=Analysis -N 100 --resource-provides=usage_model=DEDICATED,OPPORTUNISTIC --OS=SL6 file://`pwd`/run_job_diag.sh
cp jobOptions.txt jobOptions1.txt
cp run_job_diag.sh run_job_diag1.sh

sed -i -e 's/MBProcess=0/MBProcess=1/g' jobOptions1.txt
mv jobOptions1.txt /pnfs/lar1nd/scratch/users/dcianci/fits/jobOptions.txt
cp run_job_diag.sh run_job_diag1.sh
sed -i -e 's/VMode=""/VMode="MBProcess/g' run_job_diag1.sh
jobsub_submit -G lar1nd --role=Analysis -N 100 --resource-provides=usage_model=DEDICATED,OPPORTUNISTIC --OS=SL6 file://`pwd`/run_job_diag.sh
cp jobOptions.txt jobOptions1.txt
cp run_job_diag.sh run_job_diag1.sh

sed -i -e 's/MBProcessNubar=0/MBProcessNubar=1/g' jobOptions1.txt
mv jobOptions1.txt /pnfs/lar1nd/scratch/users/dcianci/fits/jobOptions.txt
cp run_job_diag.sh run_job_diag1.sh
sed -i -e 's/VMode=""/VMode="MBProcessNubar/g' run_job_diag1.sh
jobsub_submit -G lar1nd --role=Analysis -N 100 --resource-provides=usage_model=DEDICATED,OPPORTUNISTIC --OS=SL6 file://`pwd`/run_job_diag.sh
cp jobOptions.txt jobOptions1.txt
cp run_job_diag.sh run_job_diag1.sh

sed -i -e 's/NUMIProcess=0/NUMIProcess=1/g' jobOptions1.txt
mv jobOptions1.txt /pnfs/lar1nd/scratch/users/dcianci/fits/jobOptions.txt
cp run_job_diag.sh run_job_diag1.sh
sed -i -e 's/VMode=""/VMode="NUMIProcess/g' run_job_diag1.sh
jobsub_submit -G lar1nd --role=Analysis -N 100 --resource-provides=usage_model=DEDICATED,OPPORTUNISTIC --OS=SL6 file://`pwd`/run_job_diag.sh
cp jobOptions.txt jobOptions1.txt
cp run_job_diag.sh run_job_diag1.sh

sed -i -e 's/MINOSProcess=0/MINOSProcess=1/g' jobOptions1.txt
mv jobOptions1.txt /pnfs/lar1nd/scratch/users/dcianci/fits/jobOptions.txt
cp run_job_diag.sh run_job_diag1.sh
sed -i -e 's/VMode=""/VMode="MINOSProcess/g' run_job_diag1.sh
jobsub_submit -G lar1nd --role=Analysis -N 100 --resource-provides=usage_model=DEDICATED,OPPORTUNISTIC --OS=SL6 file://`pwd`/run_job_diag.sh
cp jobOptions.txt jobOptions1.txt
cp run_job_diag.sh run_job_diag1.sh

sed -i -e 's/GALLIUMProcess=0/GALLIUMProcess=1/g' jobOptions1.txt
mv jobOptions1.txt /pnfs/lar1nd/scratch/users/dcianci/fits/jobOptions.txt
cp run_job_diag.sh run_job_diag1.sh
sed -i -e 's/VMode=""/VMode="GALLIUMProcess/g' run_job_diag1.sh
jobsub_submit -G lar1nd --role=Analysis -N 100 --resource-provides=usage_model=DEDICATED,OPPORTUNISTIC --OS=SL6 file://`pwd`/run_job_diag.sh
cp jobOptions.txt jobOptions1.txt
cp run_job_diag.sh run_job_diag1.sh

sed -i -e 's/XSECProcess=0/XSECProcess=1/g' jobOptions1.txt
mv jobOptions1.txt /pnfs/lar1nd/scratch/users/dcianci/fits/jobOptions.txt
cp run_job_diag.sh run_job_diag1.sh
sed -i -e 's/VMode=""/VMode="XSECProcess/g' run_job_diag1.sh
jobsub_submit -G lar1nd --role=Analysis -N 100 --resource-provides=usage_model=DEDICATED,OPPORTUNISTIC --OS=SL6 file://`pwd`/run_job_diag.sh
cp jobOptions.txt jobOptions1.txt
cp run_job_diag.sh run_job_diag1.sh

sed -i -e 's/MBDISProcess=0/MBDISProcess=1/g' jobOptions1.txt
mv jobOptions1.txt /pnfs/lar1nd/scratch/users/dcianci/fits/jobOptions.txt
cp run_job_diag.sh run_job_diag1.sh
sed -i -e 's/VMode=""/VMode="MBDISProcess/g' run_job_diag1.sh
jobsub_submit -G lar1nd --role=Analysis -N 100 --resource-provides=usage_model=DEDICATED,OPPORTUNISTIC --OS=SL6 file://`pwd`/run_job_diag.sh
cp jobOptions.txt jobOptions1.txt
cp run_job_diag.sh run_job_diag1.sh

sed -i -e 's/MBDISProcessNubar=0/MBDISProcessNubar=1/g' jobOptions1.txt
mv jobOptions1.txt /pnfs/lar1nd/scratch/users/dcianci/fits/jobOptions.txt
cp run_job_diag.sh run_job_diag1.sh
sed -i -e 's/VMode=""/VMode="MBDISProcessNubar/g' run_job_diag1.sh
jobsub_submit -G lar1nd --role=Analysis -N 100 --resource-provides=usage_model=DEDICATED,OPPORTUNISTIC --OS=SL6 file://`pwd`/run_job_diag.sh
cp jobOptions.txt jobOptions1.txt
cp run_job_diag.sh run_job_diag1.sh
