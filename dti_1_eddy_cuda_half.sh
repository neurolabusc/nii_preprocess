#!/bin/bash
#This script processes DTI images using FSL tools
# dti_core.sh <dti> <dtir>
#  dti  : diffusion image (4D NIfTI img.nii.gz with img.bvec, img.bval)
#  dtir : (optional) as dti but reversed phase-encoding polarity
#         if <dti> but no <dtir> then eddy_correct used instead of topup+eddy
#Processing steps
# 1: undistort DTI data"
#     -use both TOPUP and EDDY if both <dti> and <dtir> are provided
#     -if only <dti> provided, EDDY_CORRECT data
# 2: diffusion modeling
#      -DTIFIT used to create fractional anisotropy (FA)
#      -FSLMATHS used to threshold FA (mask for subsequent analyses)
# 3: model fibers
#      -BEDPOSTX
# ./dti_1_eddy.sh "/Volumes/catdd_projects/Ext_Collab/Chris_Rorden/topup/test/DTIA_LM1001" "/Volumes/catdd_projects/Ext_Collab/Chris_Rorden/topup/test/DTIP_LM1001"

if [ $# -lt 1 ]; then
	echo "Please call this script with input arguments"
	echo "Example: single DTI scan from healthy participant"
	echo ' ./dti_1_eddy.sh "DTIA_LM1001"'
	echo "Example: opposite phase DTI scans"
	echo ' ./dti_1_eddy.sh "DTIA_LM1001" "DTIP_LM1001" '
	exit 1
fi

#much better way to detect CUDA version - CR and DPR 20200327
if ! [ -x "$(which nvcc)" ] ; then
    echo "Unable to find nvcc: check CUDA is installed correctly"
    exit 1
fi
#cudav=$(nvcc --version | grep "release" | awk '{ print $6 }' | cut -c2- ) 

eddyExeName=eddy_cuda
if ! [ -x "$(command -v $eddyExeName)" ]; then
  echo 'You may need to run configure_eddy.sh
  echo ' https://github.com/neurolabusc/gpu_test/tree/master/eddy_setup'
  echo 'Error: ' $eddyExeName ' is not installed.' >&2
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
nvol=$(fslnvols $dti)
if [ $nvol -lt 7 ]; then
 echo "Error: DTI has $nvol volumes: must be 4D image with at least 7 volumes $dti"
 exit 1
fi
echo "Filenames dti= $dti dtir=$dtir"

#process dti data
dti_bvec=${dti}.bvec
dti_bval=${dti}.bval
if [ ${#dtir} -eq 0 ]; then  #only given a single DTI
	dti_bvecm=${dti}.bvec #dti values if no dtir, else merged dti/dtir
	dti_bvalm=${dti}.bval #dti values if no dtir, else merged dti/dtir
else
	dti_bvecm=${dti}both.bvec
	dti_bvalm=${dti}both.bval
fi
echo $dti_bval
minBval=$(awk '{m=$1;for(i=1;i<=NF;i++)if($i>=0 && $i<m)m=$i;print m}' $dti_bval)
minBvalIdx=$(awk '{m=$1; idx=1;for(i=1;i<=NF;i++){if($i>=0 && $i<m){m=$i;idx=i}}print idx}' $dti_bval)
#FSL indices volumes from 0
minBvalIdx0=`expr $minBvalIdx - 1`
echo "minimum BValue is $minBval (volume $minBvalIdx)"
if [ $minBval -gt 10 ] #undistortion assumes first volume has b-value of zero
then
	echo "Error: no B=0 volumes"
	exit 1
fi

# read -a bval0 <$dti_bval
# #next line used to be "$bval0 -ne 0", changed since CMRR correctly reports small B values
# if [ $bval0 -gt 10 ] #undistortion assumes first volume has b-value of zero
# then
# 	echo "Error: first bvalue should be 0 not $bval0"
#	exit 1
# fi
dti_b=${dti}b #brain-extracted dti
dti_u=${dti}u #undistorted dti
dti_fa=${dti}_FA #dti fractional anisotropy map
dti_faThr=${dti}_FA_thr #thresholded dti fractional anisotropy map
dti_faEro=${dti}_FA_ero #eroded dti fractional anisotropy map

#########
echo "*** EDDY: half-sphere using --slm=linear https://github.com/MRtrix3/mrtrix3/issues/1138"
if [ ${#dtir} -eq 0 ]; then  #only given a single DTI
        #https://www.jiscmail.ac.uk/cgi-bin/webadmin?A2=FSL;dd9c03b7.1504
	dti_b0=${dti}b0
	#echo fslroi $dti $dti_b0 $minBvalIdx0 1
	fslroi $dti $dti_b0 $minBvalIdx0 1
	#echo bet $dti_b0 $dti_b  -f 0.2 -R -n -m
	bet $dti_b0 $dti_b  -f 0.2 -R -n -m
	dti_b=${dti}b_mask #masked brain-extracted dti
	#create acq_param: dummies as we will not run TOPUP
	dti_txt=${dti}_acq_param.txt
	printf "0 1 0 0.03388\n0 -1 0 0.03388\n" > $dti_txt
	#create index files: all = 1 as we will not run TOPUP
	dti_txt2=${dti}_index.txt
	nvol=$(fslnvols $dti)
	indx=""
	for ((i=1; i<=nvol; i+=1)); do indx="$indx 1"; done
	echo $indx > $dti_txt2
	echo $eddyExeName  --slm=linear --imain=$dti --mask=$dti_b --acqp=$dti_txt --index=$dti_txt2 --bvecs=$dti_bvec --bvals=$dti_bval --repol --out=$dti_u
	time $eddyExeName  --slm=linear --imain=$dti --mask=$dti_b --acqp=$dti_txt --index=$dti_txt2 --bvecs=$dti_bvec --bvals=$dti_bval --repol --out=$dti_u
	#cr 2017: use rotated vectors
	dti_bvec=${dti}u.eddy_rotated_bvecs
else #dual DTI: run topup
	echo "1 TOPUP+EDDY: undistort DTI data"
	#http://fsl.fmrib.ox.ac.uk/fsl/fslwiki/TOPUP/ExampleTopupFollowedByApplytopup
	dti_b0=${dti}b0
	dtir_b0=${dtir}b0
	both_b0=${dti}b0m #merged
	#fslroi $dti $dti_b0 0 1
	#fslroi $dti $dti_b0 $minBvalIdx0 $minBvalIdx
	fslroi $dti $dti_b0 $minBvalIdx0 1
	#fslroi $dtir $dtir_b0 0 1
	dtir_bval=${dtir}.bval
	minBval=$(awk '{m=$1;for(i=1;i<=NF;i++)if($i>=0 && $i<m)m=$i;print m}' $dtir_bval)
	if [ $minBval -gt 10 ] #undistortion assumes first volume has b-value of zero
	then
		echo "Error: reversed polarity scan has no B=0 volumes"
		exit 1
	fi
	minBvalIdx=$(awk '{m=$1; idx=1;for(i=1;i<=NF;i++){if($i>=0 && $i<m){m=$i;idx=i}}print idx}' $dtir_bval)
	minBvalIdx0=`expr $minBvalIdx - 1`
	#fslroi $dtir $dtir_b0 $minBvalIdx0 $minBvalIdx
	fslroi $dtir $dtir_b0 $minBvalIdx0 1
	fslmerge -t $both_b0 $dti_b0 $dtir_b0
	echo 'topup assuming flip in 2nd dimension and a readout of 0.3388x'
	dti_txt=${dti}_acq_param.txt
	printf "0 1 0 0.03388\n0 -1 0 0.03388\n" > $dti_txt
	dti_t=${dti}tp #topup dti
	dti_tb0=${dti}tp #topup b0 map
	#topup requires about 6 minutes for LIME data
	echo "topup --imain=$both_b0 --datain=$dti_txt --config=b02b0.cnf --out=$dti_t  --iout=$dti_tb0"
	time topup --imain=$both_b0 --datain=$dti_txt --config=b02b0.cnf --out=$dti_t  --iout=$dti_tb0
	#applytopup --imain=$dti,$dtir --inindex=1,2 --datain=$dti_txt --topup=$dti_t --out=$dti_u
	#http://fsl.fmrib.ox.ac.uk/fsl/fslwiki/EDDY/UsersGuide
	fslmaths $dti_tb0 -Tmean $dti_tb0
	bet $dti_tb0 $dti_b  -f 0.2 -R -n -m
	dti_b=${dti}b_mask #masked brain-extracted dti
	dti_txt2=${dti}_index.txt
	nvol=$(fslnvols $dti)
	nvolr=$(fslnvols $dtir)
	indx=""
	for ((i=1; i<=nvol; i+=1)); do indx="$indx 1"; done
	for ((i=1; i<=nvolr; i+=1)); do indx="$indx 2"; done
	echo $indx > $dti_txt2
	paste ${dti}.bvec ${dtir}.bvec > $dti_bvecm
	paste ${dti}.bval ${dtir}.bval > $dti_bvalm
	dti_merge=${dti}both #merged
	fslmerge -t $dti_merge $dti $dtir
	#For GPU eddy_cuda7.0 instead of eddy_openmp
	#eddy takes about 11 minutes for LIME data
	#CR 3/2017: with FSL 5.0.10 and later, use "--repol": https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/eddy/UsersGuide
	#  https://fsl.fmrib.ox.ac.uk/fsldownloads/patches/eddy-patch-fsl-5.0.9/centos6/
	echo $eddyExeName --slm=linear --imain=$dti_merge --mask=$dti_b --acqp=$dti_txt --index=$dti_txt2 --bvecs=$dti_bvecm --bvals=$dti_bvalm --topup=$dti_t --repol --out=$dti_u
	time $eddyExeName --slm=linear --imain=$dti_merge --mask=$dti_b --acqp=$dti_txt --index=$dti_txt2 --bvecs=$dti_bvecm --bvals=$dti_bvalm --topup=$dti_t --repol --out=$dti_u
	#use merged dataset for dtifit
	#dti_bvec=$dti_bvecm
	#cr 2017: use rotated vectors
	dti_bvec=${dti}u.eddy_rotated_bvecs
	dti_bval=$dti_bvalm
fi #if single DTI else dual DTI

#########
echo "2 DTIFIT + THRESHOLD FA MAPS : Compute anisotropy"
echo dtifit --data=$dti_u --out=$dti --mask=$dti_b --bvecs=$dti_bvec --bvals=$dti_bval
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
