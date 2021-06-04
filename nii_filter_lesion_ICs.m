function [prefix] = nii_filter_lesion_ICs (lesion_fname, fmri_fname, TR)
fsl_path = ([getenv('FSLDIR') '/bin']);
fsl_outtype = getenv('FSLOUTPUTTYPE');
if isempty(fsl_outtype)
  fsl_outtype = 'NIFTI';
  setenv ('FSLOUTPUTTYPE', 'NIFTI');
end
fprintf ('Identifying and removing independent components associated with lesion\n');
% 0. try to get the TR from the fMRI NIFTI header
fmri_hdr = spm_vol (fmri_fname);
if ~exist ('TR') || (TR == 0)
    temp = fieldnames (fmri_hdr(1).private);
    if ~isempty (fmri_hdr(1).private.timing) && ~isempty (find (strcmp (temp, 'timing')))
        TR = fmri_hdr(1).private.timing.tspace;
    else
        if strcmp (rest_fname(1:4), 'dswa')
            TR = 1.85;
        else
            TR = 1.65;
        end
        warning ('%s: TR not specified in the header, assuming %f\n', fmri_fname, TR);
    end
end
% 1. run MELODIC ICA
systemSub ('sh -c ". ${FSLDIR}/etc/fslconf/fsl.sh"');
command = sprintf ('%s/melodic -i %s -o melodic_output --tr=%.2f', fsl_path, fmri_fname, TR);
systemSub (command);
if strcmp (fsl_outtype, 'NIFTI_GZ')
    system ('gunzip -f melodic_output/melodic_IC.nii.gz');
end
% 2. find the ICs that overlap with the lesion mask
% overlap is measured with Jaccard index: length (intersection) / length (union)
% significant overlap: Jaccard index > 5%
lesion_hdr = spm_vol (lesion_fname);
lesion_img = spm_read_vols (lesion_hdr);
all_ic_hdr = spm_vol ('melodic_output/melodic_IC.nii');
all_ic_img = spm_read_vols (all_ic_hdr);
for i = 1:size (all_ic_img, 4)
    ic_img = squeeze (all_ic_img (:, :, :, i));
    if exist('prctile','file')
        thresh = prctile (ic_img (:), 97.5); % threshold at p<=0.05, 2-tailed
    else
        thresh = nii_prctile (ic_img (:), 97.5); % threshold at p<=0.05, 2-tailed
    end
    [~, r_ic_img] = nii_reslice_target (all_ic_hdr(i), ic_img, lesion_hdr, 1);
    t_ic_img = (r_ic_img >= thresh | r_ic_img <= -thresh);
    jaccard(i) = sum (t_ic_img(:) & lesion_img(:)) / sum (t_ic_img (:) | lesion_img(:));
end
lesion_ic = find (jaccard > 0.05);
if isempty (lesion_ic) % no ICs overlap with lesion
    prefix = '';
    return;
end
lesion_ic_str = sprintf ('%d,', lesion_ic);
lesion_ic_str = lesion_ic_str (1:length(lesion_ic_str)-1);
% 3. filter the lesion ICs from the fMRI data
prefix = 'm'; % for "melodic"
[p, n, x] = fileparts (fmri_fname);
if ~isempty (p)
    p = [p filesep];
end
filtered_fname = [p prefix n x];
command = sprintf ('%s/fsl_regfilt  -i %s -o %s -d melodic_output/melodic_mix -f %s\n', fsl_path, fmri_fname, filtered_fname, lesion_ic_str);
systemSub (command);
if strcmp (fsl_outtype, 'NIFTI_GZ')
    system (['gunzip -f ' filtered_fname '.gz']);
end
% 4. clear MELODIC output
systemSub ('rm -r melodic_output'); %% commented out for now
%end nii_filter_lesion_ICs()

function status  = systemSub (cmd)
% Save library paths
MatlabPath = getenv('LD_LIBRARY_PATH');
% Make Matlab use system libraries
setenv('LD_LIBRARY_PATH',getenv('PATH'));
fprintf('%s\n',cmd);
status = system(cmd);
% Reassign old library paths
setenv('LD_LIBRARY_PATH',MatlabPath);
%systemSub