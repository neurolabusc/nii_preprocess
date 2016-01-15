#!/bin/bash
#This script processes DTI images using FSL tools
# dti_3_tract.sh <dti>
#  dti  : diffusion image (4D NIfTI img.nii.gz with img.bvec, img.bval)
#Processing steps
# 	1: BEDPOST
# 	2: PROBTRACKX
#   5: QUANTIFY FIBERS

#For details on FSL scripting see
#http://fsl.fmrib.ox.ac.uk/fslcourse/lectures/scripting/all.htm 
# 
# ./dti_3_tract.sh "DTIA_LM1001" 

BASEDIR=$(dirname $0)
if [ "$BASEDIR" = "." ]; then 
	BASEDIR=$(pwd)
fi

if [ $# -lt 1 ]; then
	echo "Please call this script with input arguments"
	echo "Example: single DTI scan from healthy participant"
	echo ' ./dti_travis.sh "DTIA_LM1001"'	
	exit 1
fi
dti=`$FSLDIR/bin/remove_ext $1`
pth=$(dirname $dti)
if [ "$pth" = "." ]; then 
	pth=$(pwd)
fi
cd $pth
echo "Filenames dti=\"$dti\" "

SLEEP_INTERVAL=2
dti_b=${dti}b #brain-extracted dti
dti_u=${dti}u #undistorted dti
dti_fa=${dti}_FA #dti fractional anisotropy map
dti_faThr=${dti}_FA_thr #thresholded dti fractional anisotropy map
dti_bvecm=${dti}both.bvec
dti_bvalm=${dti}both.bval
if ! [ -a $dti_bvecm ]; then
	echo "Could not find \"$dti_bvecm\" (assume single session)."
	dti_bvecm=${dti}.bvec
	dti_bvalm=${dti}.bval
	if ! [ -a $dti_bvecm ]; then
		echo "Could not find \"$dti_bvecm\" (aborting)!"
		exit 1
	fi
fi
template_roiW=${dti}_roi #region of interest warped to DTI
template_roiWThr=${template_roiW}_thr #FA-thresholded region of interest warped to DTI
if ! [ -a $template_roiWThr ]; then
	fslmaths $template_roiW -mas $dti_faThr $template_roiWThr
	echo "Creating thresholded image \"$template_roiWThr\""
fi

seedlist=${pth}/seeds.txt
#paste can put a tab between merged files - remove them
#cat $dti_bvecm | tr -d "\t\r" > $dti_bvecm
#cat $dti_bvalm | tr -d "\t\r" > $dti_bvalm

###########################################################################
## BEDPOST works as is assuming that fsl_sub is configured correctly
## and SGE_ROOT is set
###########################################################################
if false; then #ssssssss

	if ! [ -a $template_roiWThr ]; then
	echo
	echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
	echo "1) BEDPOSTX: Model fibers"
	bed_dir=${pth}/bedpost
	mkdir $bed_dir
	imcp $dti_u $bed_dir/data
	[ $FSLOUTPUTTYPE = "NIFTI_GZ" ] && [ -f $bed_dir/data.nii ] && rm $bed_dir/data.nii
	imcp $dti_faThr $bed_dir/nodif_brain_mask
	[ $FSLOUTPUTTYPE = "NIFTI_GZ" ] && [ -f $bed_dir/nodif_brain_mask.nii ] && rm $bed_dir/nodif_brain_mask.nii
	cp $dti_bvecm $bed_dir/bvecs
	cp $dti_bvalm $bed_dir/bvals
	START=$(date +%s)
	bedpostx $bed_dir
	# wait while bedpost runs
	while [ ! -f ${bed_dir}.bedpostX/xfms/eye.mat ]; do
		sleep $SLEEP_INTERVAL
	done
	sleep $SLEEP_INTERVAL
	END=$(date +%s)
	DIFF=$(( $END - $START ))
	echo ".........."
	echo "bedpost took $DIFF seconds"
	fi 

fi #ssssssss


echo
echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
echo "2) PROBTRACKX: Create seed data"
#prepare a seed list for probtrackx
mask_dir=${pth}/masks
[ -f $seedlist ] && rm $seedlist
[ -d $mask_dir ] && rm -rf $mask_dir 2>/dev/null
mkdir -p $mask_dir	
for i in $(seq -w 1 189); do
	maski=${mask_dir}/$i
	fslmaths $template_roiWThr -thr $i -uthr $i -bin $maski -odt char
	mu=$(fslstats $maski -M)
	if awk -v x=$mu -v y=0 'BEGIN { exit (x > y) ? 0 : 1 }'; then   
	    printf "%s\n" "$maski" >> $seedlist
	else
	   rm ${maski}.nii*
	fi
done
# setup probtrackx
bed_dir=${pth}/bedpost
#TODO: Determine optimum number of jobs
n_jobs=50 
num_samples=5000
brain_mask=${bed_dir}.bedpostX/nodif_brain_mask
merged=${bed_dir}.bedpostX/merged
prob_opts="-s $merged -m $brain_mask --opd --pd -l -c 0.2 --distthresh=0"
prob_dir=${pth}/probtrackx
rm -rf $prob_dir 2>/dev/null
mkdir $prob_dir
mkdir $prob_dir/logs
tmp_dir=$prob_dir/tmp
START=$(date +%s)
# for each seed region...
while read seed_mask; do
    echo " - tracking from seed_mask: $seed_mask"
    # create a temporary staging area and set parameters
    rm -rf $tmp_dir 2>/dev/null
    mkdir $tmp_dir
    vox_count=$(fslstats $seed_mask -V | cut -d ' ' -f1)
    #TODO: handle case where num_samples % n_jobs != 1
    samples_per_run=$(( $num_samples / $n_jobs ))
    #TODO: est time should be based on voxels * samples_per_run
    est_time=1
    est_time2=2
    # write each proposed job as a line in the command file
    for j_idx in $(seq -w 1 $n_jobs); do
        echo "${FSLDIR}/bin/probtrackx2 -x $seed_mask --dir=$tmp_dir --forcedir -o $j_idx -P $samples_per_run $prob_opts" >> $tmp_dir/commands.txt
    done
    echo "${FSLDIR}/bin/probtrackx2 -x $seed_mask --dir=$prob_dir --forcedir -o $seed_mask -P $num_samples $prob_opts" >> $prob_dir/commands.txt
    probtrackxid=`${FSLDIR}/bin/fsl_sub -l $prob_dir/logs -N probtrackx -T $est_time -t $tmp_dir/commands.txt`
    sentinelid=`${FSLDIR}/bin/fsl_sub -l $prob_dir/logs -j $probtrackxid -T $est_time2 -N sentinel touch $tmp_dir/complete`
    # wait until all fdt images have been generated
    complete=0
    while [ $complete = 0 ]; do
        [ -f $tmp_dir/complete ] && complete=1
        sleep $SLEEP_INTERVAL
    done
    # merge the temporary fdt's into a single fdt
    fdt_out=$prob_dir/`basename $seed_mask`
    fslmaths $template_roiWThr -mul 0 $fdt_out
    for j_idx in $(seq -w 1 $n_jobs); do
        fslmaths $fdt_out -add $tmp_dir/$j_idx $fdt_out
    done
    # remove the temporary staging space
    rm -rf $tmp_dir 2>/dev/null
done < $seedlist
END=$(date +%s)
DIFF=$(( $END - $START ))
echo ".........."
echo "PROBTRACKX took $DIFF seconds"


exit 0
#DO NOT RUN QUANT - FAILS ON OSX: "declare -A" not supported
echo
echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
echo "3) QUANTIFICATION: Count Fibers"
fiber_count_mat_file=$prob_dir/fiber_count.mat
density_mat_file=$prob_dir/density.mat
unset fiber_count_mat
unset density_mat
declare -A fiber_count_mat
declare -A density_mat
START=$(date +%s)
for i in $( seq -w 1 189); do
    echo " - populating row($i)"
    i_fdt=$prob_dir/$i
    i_mask=$mask_dir/$i
    if [ `imtest $i_fdt` = 1 ]; then
        i_exists=1
        i_nvox=$(fslstats $i_mask -V | cut -d ' ' -f1)
    else
        i_exists=0
    fi
    for j in $( seq -w $i 189 ); do
        j_fdt=$prob_dir/$j
        j_mask=$mask_dir/$j
        if [ $i -ne $j ] && [ $i_exists = 1 ] && [ `imtest $j_fdt` = 1 ]; then
            j_nvox=$(fslstats $j_mask -V | cut -d ' ' -f1)
            ij_mean=$(fslstats $i_fdt -k $j_mask -M 2>/dev/null)
            ji_mean=$(fslstats $j_fdt -k $i_mask -M 2>/dev/null)
            ij_sum=$(printf "%0.f" $(echo "$ij_mean * $j_nvox" | bc -l))
            ji_sum=$(printf "%0.f" $(echo "$ji_mean * $i_nvox" | bc -l))
            fiber_count=$(( $ij_sum + $ji_sum ))
            normalizing_factor=$( echo "( $i_nvox + $j_nvox ) * ( $num_samples + 1 )" | bc -l )
            density=$( awk "BEGIN {printf \"%e\", $fiber_count/$normalizing_factor}" )
        else
            fiber_count=0
            density=$(printf "%e" 0)
        fi
        fiber_count_mat[$i,$j]=$fiber_count
        fiber_count_mat[$j,$i]=$fiber_count
        density_mat[$i,$j]=$density
        density_mat[$j,$i]=$density
    done
done
# write out matrices
echo " - writing out matrices"
for i in $( seq -w 1 189 ); do
    for j in $( seq -w 1 189 ); do
        printf "%d " ${fiber_count_mat[$i,$j]} >> $fiber_count_mat_file
        printf "%e " ${density_mat[$i,$j]} >> $density_mat_file
    done
    printf "\n" >> $fiber_count_mat_file
    printf "\n" >> $density_mat_file
done
# remove trailing whitespace at the end of each matrix row
sed -i 's/ $//g' $fiber_count_mat_file
sed -i 's/ $//g' $density_mat_file

END=$(date +%s)
DIFF=$(( $END - $START ))
echo ".........."
echo "QUANTIFICATION took $DIFF seconds"
echo
exit 0
