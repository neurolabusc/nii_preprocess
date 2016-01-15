#!/bin/bash
#This script processes DTI images using FSL tools
# dti.sh <dti> <dtir> <lesion> <anat>
#  dti  : diffusion image (4D NIfTI img.nii.gz with img.bvec, img.bval)
#  dtir : (optional) as dti but reversed phase-encoding polarity
#         if <dti> but no <dtir> then eddy_correct used instead of topup+eddy
#  lesion : (optional) binary map indicating region of injury (NIfTI format)
#  anat   : (optional) anatomical image used to draw lesion
#           if <lesion> provided but no <anat>, one assumes lesion in <dti> space
#Processing steps
# 1/2/3: run dtix.sh
# 	1: undistort DTI data"
#     -use both TOPUP and EDDY if both <dti> and <dtir> are provided
#     -if only <dti> provided, EDDY_CORRECT data
# 	2: diffusion modeling
#      -DTIFIT used to create fractional anisotropy (FA)
#      -FSLMATHS used to threshold FA (mask for subsequent analyses)
# 	3: model fibers
#      -BEDPOSTX
# 4: warp regions of interest to native DTI space
#      -If both <lesion> and <anat> are provided, FLIRT used to warp <lesion> to <dti>
#      -Region of interest warped from standard to dti space using FLIRT
# 5: estimate connections between regions
#      -FSLMATHS creates one image per region
#      -PROBTRACKX evaluates connections

#For details on FSL scripting see
#http://fsl.fmrib.ox.ac.uk/fslcourse/lectures/scripting/all.htm 
# 
# ./dti.sh "DTIA_LM1001" "" "LS_LM1001" "T1_LM1001"
# ./dti.sh "/Users/rorden/lxb/test/DTIA_LM1001" "" "/Users/rorden/lxb/test/LS_LM1001" "/Users/rorden/lxb/test/T1_LM1001"
#./dti.sh "/Volumes/catdd_projects/Ext_Collab/Chris_Rorden/fx/test/DTIA_LM1001" "/Volumes/catdd_projects/Ext_Collab/Chris_Rorden/fx/test/DTIP_LM1001" "/Volumes/catdd_projects/Ext_Collab/Chris_Rorden/fx/test/LS_LM1001" "/Volumes/catdd_projects/Ext_Collab/Chris_Rorden/fx/test/T1_LM1001"
BASEDIR=$(dirname $0)
if [ "$BASEDIR" = "." ]; then 
	BASEDIR=$(pwd)
fi

if [ $# -lt 1 ]; then
	echo "Please call this script with 4 input arguments"
	echo "Example: single DTI scan from healthy participant"
	echo ' ./dti.sh "DTIA_LM1001"'	
	echo "Example: opposite phase DTI scans"
	echo ' ./dti.sh "DTIA_LM1001" "DTIP_LM1001" '
	echo "Example: scan from patient with lesion drawn on DTI"
	echo ' ./dti.sh "DTIA_LM1001" "DTIP_LM1001" "LS_LM1001"'
	echo "Example: scan from patient with lesion drawn on T1"
	echo ' ./dti.sh "DTIA_LM1001" "DTIP_LM1001" "LS_LM1001" "T1_LM1001"'
	echo "Example: scan from patient with lesion drawn on T1 (single DTI)"
	echo ' ./dti.sh "DTIA_LM1001" "" "LS_LM1001" "T1_LM1001"'
	exit 1
fi
dti=`$FSLDIR/bin/remove_ext $1`
pth=$(dirname $dti)
if [ "$pth" = "." ]; then 
	pth=$(pwd)
fi
cd $pth
if [ $# -gt 1 ]; then
	dtir=`$FSLDIR/bin/remove_ext $2`
else
	dtir=""
fi 
if [ $# -gt 2 ]; then
	lesion=`$FSLDIR/bin/remove_ext $3`
else
	lesion=""
fi 
if [ $# -gt  3 ]; then
	anat=`$FSLDIR/bin/remove_ext $4`
else
	anat=""
fi
echo "Filenames dti=\"$dti\" dtir=\"$dtir\" lesion=\"$lesion\" anat=\"$anat\""

#preprocess

${BASEDIR}/dti_core.sh $dti $dtir

if [ $? -ne 0 ]; then 
 echo "Error: dti_core script reported a failure"
 exit 1
fi

echo "4 FLIRT: Warp images to DTI space"
#warp lesion to DTI if it was drawn on a different modality
dti_faEro=${dti}_FA_ero #eroded dti fractional anisotropy map
if [ ${#lesion} -ne 0 ] &&  [ ${#anat} -ne 0 ]; then 
	echo " flirt anat -> dti, apply transform to lesion"
	anat_m=${anat}.mt
	anat_w=${anat}w
	flirt -in $anat -ref $dti_faEro -out $anat_w -omat $anat_m -bins 256 -cost corratio -searchrx -45 45 -searchry -45 45 -searchrz -45 45 -dof 12 -interp trilinear
	lesion_w=${lesion}w
	flirt -in $lesion -ref $dti_faEro -out $lesion_w -applyxfm -init $anat_m -interp nearestneighbour
	anat_w=${dti_fa} #<--- remove this line for T1 normalization
else
	lesion_w=${lesion}
	anat_w=${dti_fa}
fi

#warp template to native DTI space 
template_anat=${BASEDIR}/mni152_2009_1mm.nii #<---change this for T1 normalization --> template_anat=${BASEDIR}/sct1_1mm.nii
template_roi=${BASEDIR}/jhu1mm.nii
template_roiM=${dti}_std2native.mt
template_roiW=${dti}_roi #region of interest warped to DTI
template_roiWThr=${template_roiW}_thr #FA-thresholded region of interest warped to DTI
dti_faThr=${dti}_FA_thr #thresholded dti fractional anisotropy map
if [ ${#lesion} -ne 0 ]; then
	mask=${lesion}_msk
	fslmaths $lesion_w -binv $mask
	echo " flirt roi_anat -> dti, apply transform to roi (lesion masked)"
	flirt -in $template_anat -ref $dti_faEro -refweight $mask -out $template_roiW -omat $template_roiM -bins 256 -cost corratio -searchrx -45 45 -searchry -45 45 -searchrz -45 45 -dof 12 -interp trilinear
 	flirt -in $template_roi -ref $dti_faEro -out $template_roiW -applyxfm -init $template_roiM -interp nearestneighbour 
else
	echo " flirt roi_anat -> dti, apply transform to roi"
	flirt -in $template_anat -ref $dti_faEro -out $template_roiW -omat $template_roiM -bins 256 -cost corratio -searchrx -45 45 -searchry -45 45 -searchrz -45 45 -dof 12 -interp trilinear
 	flirt -in $template_roi -ref $dti_faEro -out $template_roiW -applyxfm -init $template_roiM -interp nearestneighbour 
fi
#I would have liked a dilation, but these do not work well with indexed images: dilF is biased, dilM creates fake indices at borders, and dilD creates some odd results (ties?)
# fslmaths $template_roiW -dilM -mas $dti_faThr $template_roiWThr
fslmaths $template_roiW -mas $dti_faThr $template_roiWThr

echo "5 PROBTRACKX: Create seed data"
#prepare a seed list for probtrackx
mask_dir=${pth}/masks
seedlist=${pth}/seeds.txt
rm -f $seedlist	
mkdir -p $mask_dir	
for ((i=1; i<=189; i+=1)); do 
	maski=${mask_dir}/roi$i.nii.gz
	fslmaths $template_roiWThr -thr $i -uthr $i $maski
	max=$(fslstats $maski -M)
	if awk -v x=$max -v y=0 'BEGIN { exit (x > y) ? 0 : 1 }'; then   
		printf "%s\n" "$maski" >> $seedlist
	else
	   rm $maski
	fi
done

probtrackx2 --network -x $seedlist  -l --onewaycondition -c 0.2 -S 2000 --steplength=0.5 -P 5000 --fibthresh=0.01 --distthresh=0.0 --sampvox=0.0 --forcedir --opd -s /Users/rorden/lxb/sh/bedpost.bedpostX/merged -m /Users/rorden/lxb/sh/bedpost.bedpostX/nodif_brain_mask  --dir=/Users/rorden/lxb/sh/prob
exit 0
