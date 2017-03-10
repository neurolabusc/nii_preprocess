#!/bin/bash
#This script processes DTI images using FSL tools
#
#NB : THIS USES EDDY_CORRECT WHICH IS LESS CAPABLE THAN EDDY
#     YOU SHOULD ONLY USE THIS SCRIPT FOR HALF-SPHERE SAMPLING
#     OTHERWISE, USE dti_1_eddy_cuda.sh
#
# dti_1_eddy_correct.sh <dti> <dtir>
#  dti  : diffusion image (4D NIfTI img.nii.gz with img.bvec, img.bval)
#Processing steps
# 1: undistort DTI data"
#     -use both TOPUP and EDDY if both <dti> and <dtir> are provided
#     -if only <dti> provided, EDDY_CORRECT data
# 2: diffusion modeling
#      -DTIFIT used to create fractional anisotropy (FA)
#      -FSLMATHS used to threshold FA (mask for subsequent analyses)
# ./dti_1_eddy_correct.sh "/Volumes/catdd_projects/Ext_Collab/Chris_Rorden/topup/test/DTIA_LM1001"

if [ $# -lt 1 ]; then
	echo "Please call this script with input arguments"
	echo "Example: single DTI scan"
	echo ' ./dti_1_eddy_correct.sh "DTIA_LM1001"'
	exit 1
fi
dti=`$FSLDIR/bin/remove_ext $1`
pth=$(dirname $dti)
if [ "$pth" = "." ]; then
	pth=$(pwd)
fi
cd $pth
if [ $# -gt 1 ]; then
	echo "Error: For multiple images use dti_1_eddy_cuda.sh"
	exit 1
fi

nvol=$(fslnvols $dti)
if [ $nvol -lt 7 ]; then
 echo "Error: DTI has $nvol volumes: must be 4D image with at least 7 volumes $dti"
 exit 1
fi
echo "Filenames dti= $dti"
echo "minimum BValue is $minBval (volume $minBvalIdx)"
#process dti data
dti_bvec=${dti}.bvec
dti_bval=${dti}.bval
minBval=$(awk '{m=$1;for(i=1;i<=NF;i++)if($i>=0 && $i<m)m=$i;print m}' $dti_bval)
minBvalIdx=$(awk '{m=$1; idx=1;for(i=1;i<=NF;i++){if($i>=0 && $i<m){m=$i;idx=i}}print idx}' $dti_bval)

# FSL indices volumes from 0
minBvalIdx0=`expr $minBvalIdx - 1`

echo "minimum BValue is $minBval volume $minBvalIdx"
if [ $minBval -gt 10 ] #undistortion assumes first volume has b-value of zero
then
	echo "Error: no B=0 volumes"
	exit 1
fi

dti_b=${dti}b #brain-extracted dti
dti_u=${dti}u #undistorted dti
dti_fa=${dti}_FA #dti fractional anisotropy map
dti_faThr=${dti}_FA_thr #thresholded dti fractional anisotropy map
dti_faEro=${dti}_FA_ero #eroded dti fractional anisotropy map

echo "1 EDDY_CORRECT: full-sphere bvecs should use dti_1_eddy_cuda.sh"
#eddy_correct $dti $dti_u 0
eddy_correct $dti $dti_u $minBvalIdx0
bet $dti_u $dti_b  -f 0.2 -n -R -m
dti_b=${dti}b_mask #masked brain-extracted dti

#########
echo "2 DTIFIT + THRESHOLD FA MAPS : Compute anisotropy"
dtifit --data=$dti_u --out=$dti --mask=$dti_b --bvecs=$dti_bvec --bvals=$dti_bval
fslmaths $dti_fa -ero $dti_faEro
fslmaths $dti_fa -ero -thr 0.15 -bin $dti_faThr

#########
#	echo "3 BEDPOSTX : Model fibers"
#	bed_dir=${pth}/bedpost
#	#xfms/eye.mat #file created when bedpost finishes
#	mkdir $bed_dir
#	imcp $dti_u $bed_dir/data
#	imcp $dti_faThr $bed_dir/nodif_brain_mask
#	cp $dti_bvecm $bed_dir/bvecs
#	cp $dti_bvalm $bed_dir/bvals
#	START=$(date +%s)
#	bedpostx $bed_dir
#	bed_done=${bed_dir}.bedpostX/monitor
#	$bed_done
#	END=$(date +%s)
#	DIFF=$(( $END - $START ))
#	echo "bedpost took $DIFF seconds"

exit 0
