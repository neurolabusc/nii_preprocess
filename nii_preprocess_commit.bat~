#!/bin/bash
#go to nii_preprocess folder
if [[ ! -e /usr/local/nii_preprocess ]]; then
            mkdir /usr/local/nii_preprocess
fi
cd /usr/local/nii_preprocess
#commit with message
d=`date`
h=`hostname`
git commit -a -m "$h-$d"
#git commit -a
#add files in nii_preprocess folder
#git add .
git push origin master
