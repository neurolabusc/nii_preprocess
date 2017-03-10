#!/bin/bash
#go to nii_preprocess folder
if [[ ! -e /home/crlab/nii_preprocess ]]; then
            mkdir /home/crlab/nii_preprocess
fi
cd /home/crlab/nii_preprocess
#commit with message
git add .
d=`date`
h=`hostname`
git commit -a -m "$h-$d"
#git commit -a
#add files in nii_preprocess folder
#git add .
git push origin +master
