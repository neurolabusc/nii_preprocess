#!/bin/bash
#This script processes DTI images using FSL tools
# dti_2warp_template.sh <fa> <t1> <et1> 
#  dti : DTI (must have FA map named {dti}_FA_ero
#  t1 : scalp stripped T1
#  et1 : (optional) T1 with lesion filled in 
#Processing steps
# 1: scalp strip et1 using t1 (creating bet1)
# 2: coreg t1 -> fa, apply to bet1 (creating wbet1)
# 3: coreg templatet1 -> wbet1, apply to template (creating wtemplate)
#For details on FSL scripting see
#http://fsl.fmrib.ox.ac.uk/fslcourse/lectures/scripting/all.htm 

BASEDIR=$(dirname $0)
if [ "$BASEDIR" = "." ]; then 
	BASEDIR=$(pwd)
fi

pushd `dirname $0` > /dev/null
scriptpath=`pwd -P`
popd > /dev/null

if [ $# -lt 2 ]; then
	echo "Please call this script with 4-5 input arguments"
	echo "Example: DTI scan healthy stroke participant"
	echo ' ./warp_template.sh "1DTI_P084_1_FA_ero" "bT1_P084.nii"'	
	echo "Example: DTI scan from stroke participant"
	echo ' ./warp_template.sh "1DTI_P084_1_FA_ero" "bT1_P084.nii" "eT1_P084.nii"'	
	exit 1
fi
template=${scriptpath}/jhu
if [ $(${FSLDIR}/bin/imtest "$template") = 0 ] ; then 
	echo "ERROR: image '$template' doesn't exist"; 
	exit 1
fi
templatet1=${scriptpath}/jhu_t1
if [ $(${FSLDIR}/bin/imtest "$templatet1") = 0 ] ; then 
	echo "ERROR: image '$templatet1' doesn't exist"; 
	exit 1
fi
#template_nodir=`basename $template`

#template=`$FSLDIR/bin/remove_ext $1`
#templatet1=`$FSLDIR/bin/remove_ext $2`
dti=`$FSLDIR/bin/remove_ext $1`
t1=`$FSLDIR/bin/remove_ext $2`
pth=$(dirname $dti)
if [ "$pth" = "." ]; then 
	pth=$(pwd)
fi
cd $pth
dti=$(basename "$dti" )
fa=${dti}_FA_ero
if [ $# -gt  2 ]; then
	#warp t1 to match fa, apply to et1
	et1=`$FSLDIR/bin/remove_ext $3`
	bet1=t1$fa
	fslmaths "$t1" -bin -mul "$et1" "$bet1"
	xmat=${t1}.mt
	
	echo "flirt -in \"$t1\" -ref \"$fa\" -omat \"$xmat\" -usesqform   -searchrx -20 20 -searchry -20 20 -searchrz -20 20 "
	flirt -in "$t1" -ref "$fa" -omat "$xmat" -usesqform  -searchrx -20 20 -searchry -20 20 -searchrz -20 20
	rt1=r$bet1
	echo "flirt -in \"$bet1\" -ref \"$fa\" -out \"$rt1\" -applyxfm -init \"$xmat\" -interp trilinear"	
	flirt -in "$bet1" -ref "$fa" -out "$rt1" -applyxfm -init "$xmat" -interp trilinear
	
else
	#warp t1 to match fa
	rt1=rt1$fa
	echo "flirt -in \"$t1\" -ref \"$fa\" -out \"$rt1\""
	flirt -in "$t1" -ref "$fa" -out "$rt1"	
fi
#warp templatet1 to match results of step 1, apply to template
xmat=${fa}.mt
echo "flirt -in \"$templatet1\" -ref \"$rt1\" -omat \"$xmat\""
flirt -in "$templatet1" -ref "$rt1" -omat "$xmat"
rtemplate=${dti}_roi
echo "flirt -in \"$template\" -ref \"$fa\" -out \"$rtemplate\" -applyxfm -init \"$xmat\" -interp trilinear"	
flirt -in "$template" -ref "$fa" -out "$rtemplate" -applyxfm -init "$xmat" -interp nearestneighbour
exit 0
