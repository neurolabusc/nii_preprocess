function nii_rest_batch (Restnames, T1name, TRsec,SliceOrder);
%preprocesses resting state data using SPM
%  Restnames  : 4D Resting State image(s)
%  T1name     : T1-weighted anatomical scan from participant
%  TRsec      : Time between volumes in seconds
%  SliceOrder : see nii_sliceTime 0=auto-detect,1=ascend,2=descend
%Examples
% nii_rest_batch('rs.nii','t1.nii');
% nii_rest_batch('REST_LM1019.nii','T1_LM1019.nii',1.85, 2);

doSliceTime = true; %determine whether you want slice time correction 
doDeleteTemporary = true; %remove images generated during the middle stages of pipeline
resliceMM = 3; %reslicing resolution, e.g. if 3 then data will be 3x3x3 mm
if exist('spm')~=2; fprintf('%s requires SPM\n',mfilename); return; end;
if nargin <1 %no input: select ASL file[s]
 Restnames = spm_select(inf,'image','Select 4D resting state volumes[s]');
end
if nargin <2 %no input: select T1 file[s]
 T1name = spm_select(1,'image','Select T1 anatomical scan');
end
if (~exist('TRsec')) || (TRsec == 0)
        Filename = deblank (Restnames(1,:));
        TRsec =  getTRSub(Filename);
        if (TRsec ==0) 
            fprintf('%s error unable to determine TR\n',mfilename);
            return;
        end;
end
if ~exist('SliceOrder')
   SliceOrder = 2; %2 = NIFTI_SLICE_SEQ_DEC - see http://nifti.nimh.nih.gov/pub/dist/src/niftilib/nifti1.h
   fprintf('%s warning: assuming slice order is DESCENDING\n', mfilename);
   %SliceOrder = 0; %0 = auto detect
end
fprintf('Assuming time between volumes (TR) is %.3f seconds\n', TRsec);
spm('Defaults','fMRI');
spm_jobman('initcfg'); % useful in SPM8 only
clear matlabbatch
tic; %start timer
% 0 - ONCE PER INDIVIDUAL normalize T1 scan using unified normalization-segmentation...
[pthT1,namT1,extT1,vol] = spm_fileparts( deblank (T1name));
c1name = fullfile(pthT1,['c1', namT1, extT1]); %segmented gray matter
c2name = fullfile(pthT1,['c2', namT1, extT1]); %segmented white matter
c3name = fullfile(pthT1,['c3', namT1, extT1]); %segmented CSF
mT1name = priorSegnormSub(T1name); %use precomputed segmentation if available
if length(mT1name) < 1 %no previous segmentation
    mT1name = segnormSub(T1name);
    segnormwriteSub(T1name,{mT1name,c1name, c2name, c3name}, resliceMM); %warp bias corrected image
end
%REMAINING STEPS ONCE PER RESTING STATE SESSION
nses = length(Restnames(:,1));
for ses = 1 : nses
    fMRIname = deblank (Restnames(ses,:));
    [pth,nam,ext,vol] = spm_fileparts( deblank (fMRIname));
    prefix = ''; %prefixes applied to filename with operations, e.g. 'r'esliced, 's'moothed
    nVol = numVolumesSub(fMRIname);
    if nVol < 2 %this script uses 4d data
        fprintf('ERROR: %s requires 4D datasets volumes\n',mfilename);
        return;
    end;
    %1 adjust for slice time differences:
    if doSliceTime, %do we adjust for timing of 2D slices within 3D volume, uses spm_slice_timing
        prefix = [nii_sliceTime(fullfile(pth,[prefix, nam, ext]),SliceOrder, TRsec), prefix];
    end;
    %2 motion correct data:
    prefixRealign = prefix;
    meanfmri = realignSub(fullfile(pth,[prefix, nam, ext])); %does not change prefix - only new matrix
    %3 coregister fMRI to T1, then warp fMRI to standard space
    coregEstSub(T1name, meanfmri, fullfile(pth,[prefix, nam, ext]))
    segnormwriteSub(T1name, meanfmri, resliceMM);  % spm_write_sn
    prefix = [segnormwriteSub(T1name,  getsesvolsSub(fullfile(pth,[prefix, nam, ext])), resliceMM), prefix];
    %4 generate a brain mask - must load from disk, may save to memory
    maskName = makeMaskSub(strvcat(fullfile(pthT1,['wc1', namT1, extT1]),fullfile(pthT1,['wc2', namT1, extT1]),fullfile(pthT1,['wc3', namT1, extT1])),0.05);
    %5 detrend for motion and white matter fluctuations AND SMOOTH
    prefix = [nii_detrend(fullfile(pth,[prefix, nam, ext]),maskName,fullfile(pthT1,['wc2', namT1, extT1]),fullfile(pth,['rp_' prefixRealign, nam, '.txt']),3,0.8,6), prefix];
    % note: detrend PRIOR to bandpass http://www.ncbi.nlm.nih.gov/pubmed/23747457
    %6 Smooth data - no longer required (integrated with detrend)
    % prefix = [smoothSub(fullfile(pth,[prefix, nam, ext]), 6), prefix];   
    %7 do ALFF
    alfName = alffSub(fullfile(pth,[prefix, nam, ext]), TRsec, maskName);
    %8 remove frequencies
    prefix = [nii_temporalFilter(fullfile(pth,[prefix, nam, ext]), TRsec, true), prefix];
    % note: bandpass AFTER detrend http://www.ncbi.nlm.nih.gov/pubmed/23747457 
    %remove temp files
    if (doDeleteTemporary) && (length(prefix) > 1)
        for i = 2:length(prefix)
            fprintf('Deleting  %s\n',fullfile(pth,[prefix(i:length(prefix)), nam, ext]) );
            delete(fullfile(pth,[prefix(i:length(prefix)), nam, ext]) );
        end;  
    end; %if delete intermediate
end; %for each session
fprintf('Done processing %d sessions in %0.3fsec\n',nses, toc);
%%END rest_batch()

%%%%%%% SUBFUNCTIONS FOLLOW
function [bias_t1] = priorSegnormSub(t1);
%returns name of bias file if normalization was already computed
[pth,nam,ext,vol] = spm_fileparts( deblank (t1));
bias_t1 = fullfile(pth,['m', nam, ext]);
segmat =fullfile(pth,[ nam, '_seg_sn.mat']);
if  (exist(segmat) == 2)
    fprintf('Skipping normalization (will use existing values from %s)\n',segmat);
    return;
end
bias_t1 = [];
%end priorSegnormSub()

function [bias_t1] = segnormSub(t1);
%estimate unified segmentation-normalization for image t1
[pth,nam,ext,vol] = spm_fileparts( deblank (t1));
bias_t1 = fullfile(pth,['m', nam, ext]);
fprintf('Unified segmentation normalization of %s\n',t1);
matlabbatch{1}.spm.spatial.preproc.data = {t1};
matlabbatch{1}.spm.spatial.preproc.output.GM = [0 0 1];
matlabbatch{1}.spm.spatial.preproc.output.WM = [0 0 1];
matlabbatch{1}.spm.spatial.preproc.output.CSF = [0 0 1];
matlabbatch{1}.spm.spatial.preproc.output.biascor = 1;
matlabbatch{1}.spm.spatial.preproc.output.cleanup = 2;
matlabbatch{1}.spm.spatial.preproc.opts.tpm = {fullfile(spm('Dir'),'tpm','grey.nii');fullfile(spm('Dir'),'tpm','white.nii');fullfile(spm('Dir'),'tpm','csf.nii')};
matlabbatch{1}.spm.spatial.preproc.opts.ngaus = [2;2;2;4];
matlabbatch{1}.spm.spatial.preproc.opts.regtype = 'mni';
matlabbatch{1}.spm.spatial.preproc.opts.warpreg = 1;
matlabbatch{1}.spm.spatial.preproc.opts.warpco = 25;
matlabbatch{1}.spm.spatial.preproc.opts.biasreg = 0.0001;
matlabbatch{1}.spm.spatial.preproc.opts.biasfwhm = 60;
matlabbatch{1}.spm.spatial.preproc.opts.samp = 3;
matlabbatch{1}.spm.spatial.preproc.opts.msk = {''};
spm_jobman('run',matlabbatch);
%end segnormSub()

function [prefix] = segnormwriteSub(t1,images2warp, mm);
%reslice images2warp based on prior segmentation-normalization of image t1
prefix = 'w'; %'w'arped data
[pth,nam,ext,vol] = spm_fileparts( deblank (t1));
fprintf('Applying unified segmentation normalization parameters from %s to %d image[s], resliced to %fmm\n',t1,length(images2warp(:,1)),mm);
matlabbatch{1}.spm.spatial.normalise.write.subj.matname = {fullfile(pth,[ nam, '_seg_sn.mat'])};
if isa(images2warp,'char')
    matlabbatch{1}.spm.spatial.normalise.write.subj.resample = {images2warp};
else
    matlabbatch{1}.spm.spatial.normalise.write.subj.resample = images2warp;
end;
matlabbatch{1}.spm.spatial.normalise.write.roptions.preserve = 0;
matlabbatch{1}.spm.spatial.normalise.write.roptions.bb = [-78 -112 -50; 78 76 85];
matlabbatch{1}.spm.spatial.normalise.write.roptions.vox = [mm mm mm];
matlabbatch{1}.spm.spatial.normalise.write.roptions.interp = 1;
matlabbatch{1}.spm.spatial.normalise.write.roptions.wrap = [0 0 0];
matlabbatch{1}.spm.spatial.normalise.write.roptions.prefix = prefix;
spm_jobman('run',matlabbatch);
%%END segnormwriteSub()

function [nVol] = numVolumesSub (imgName);
%returns number of volumes (timepoints, directions) in a NIfTI image
[pth,nam,ext,vol] = spm_fileparts( deblank (imgName));
imgName = fullfile(pth,[nam, ext]); %remove volume is 'img.nii,1' -> 'img.nii'
hdr = spm_vol([imgName]); 
nVol = length(hdr);
%end numVolumesSub()

function [tr] =  getTRSub(fMRIname);
%Returns Repeat Time in seconds for volume fMRIname
hdr = spm_vol(fMRIname);
if isfield(hdr(1,1).private.timing,'tspace')
  tr = hdr(1,1).private.timing.tspace;
else
  fprintf('%s error: unable to determine TR for image %s\n',mfilename,fmriname); 
end
%end getTRSub()

function meanfmri = realignSub(fMRIname);
[pth,nam,ext,vol] = spm_fileparts( deblank (fMRIname));
meanfmri = fullfile(pth,['mean', nam, ext]);
fprintf('Realigning data (motion correction), creating mean image named %s\n',meanfmri);
matlabbatch{1}.spm.spatial.realign.estwrite.data = {getsesvolsSub(fMRIname)}';
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.quality = 0.9;
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.sep = 4;
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.fwhm = 5;
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.rtm = 1;
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.interp = 2;
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.wrap = [0 0 0];
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.weight = '';
matlabbatch{1}.spm.spatial.realign.estwrite.roptions.which = [0 1];
matlabbatch{1}.spm.spatial.realign.estwrite.roptions.interp = 4;
matlabbatch{1}.spm.spatial.realign.estwrite.roptions.wrap = [0 0 0];
matlabbatch{1}.spm.spatial.realign.estwrite.roptions.mask = 1;
matlabbatch{1}.spm.spatial.realign.estwrite.roptions.prefix = 'r';
spm_jobman('run',matlabbatch);
%end  mocoSub()

function [imgNameVol1] = getVol1Sub(imgName);
% input: nifti image, output: name of FIRST volume
%  example 4D file with 3 volumes input= 'img.nii', output= 'img.nii,1'
[pth,nam,ext,vol] = spm_fileparts( deblank (imgName));
imgNameVol1 = fullfile(pth,[nam, ext, ',1']);
% end getVol1Sub()

function coregEstSub(t1, meanfmri, fMRIname);
fprintf('Coregistering %s to match %s\n',meanfmri,t1);
matlabbatch{1}.spm.spatial.coreg.estimate.ref = {t1};
matlabbatch{1}.spm.spatial.coreg.estimate.source = {meanfmri};
matlabbatch{1}.spm.spatial.coreg.estimate.other = getsesvolsSub(fMRIname);
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];
spm_jobman('run',matlabbatch);
%end subfunction coregEstSub()

function maskName = makeMaskSub(imgNames, thresh);
%creates a binary mask based on a series of images (e.g. white and gray matter tissue probability maps)
%   imgNames   : list of image(s) to sum together
%   thresh     : (optional) mask includes voxels that exceed this threshold, default = 0
%Example 
% nii_makeMask('wc1t1.nii',0.5);  
% nii_makeMask(strvcat('wc1t1.nii','wc2t1.nii'),0.2);
% nii_makeMask(strvcat('wc1t1.nii','wc2t1.nii','wc3t1.nii'),0.1);
if ~exist('imgNames') %no image
 imgNames = spm_select(inf,'image','Select images to binarize');
end;
if ~exist('thresh') %no image
 thresh = 0;
end;
%load first image
hdr = spm_vol (deblank (imgNames(1,:)));
inimg = spm_read_vols (hdr);
%add additional images
nImgs = length(imgNames(:,1));
if nImgs > 1
    for i = 1 : nImgs
        addhdr = spm_vol (deblank (imgNames(i,:)));
        addimg = spm_read_vols (addhdr);
        inimg = inimg + addimg;
    end
end
%save output
[pth,nam,ext] = spm_fileparts(deblank (imgNames(1,:)));
maskName = fullfile(pth,['b' nam ext]);
outimg = zeros(size(inimg));
tmp = find((inimg)> thresh);
outimg(tmp) = 1;
hdr.fname = maskName;
spm_write_vol(hdr,outimg);
%end nii_makeMask()

function [sesvols] = getsesvolsSub(sesvol1);
% input: single volume from 4D volume, output: volume list
%  example 4D file with 3 volumes input= 'img.nii', output= {'img.nii,1';'img.nii,2';'img.nii,3'}
	[pth,nam,ext,vol] = spm_fileparts( deblank (sesvol1));
	sesname = fullfile(pth,[nam, ext]);
	hdr = spm_vol(sesname);
	nvol = length(hdr);
	if (nvol < 2), fprintf('Error 4D fMRI data required %s\n', sesname); return; end;
    sesvols = cellstr([sesname,',1']);
    for vol = 2 : nvol
        sesvols = [sesvols; [sesname,',',int2str(vol)]  ];
    end;
%end getsesvolsSub()

function outName = alffSub(name4D, kTRsec, maskName,namePrefix);
%Identify low frequency flutuations 
% name4D     : filename of 4D NIFTI image
% TRsec      : Repeat time (seconds)
% nameMask   : (optional) name for masking image
% namePrefix : (optional) text appended to output image name, default 'alf_'
%Adapts f_alff from REST toolbox for NIfTI data
% Ref: Zou QH, Zhu CZ, Yang Y, Zuo XN, Long XY, Cao QJ, Wang YF, Zang YF (2008) An improved approach to detection of amplitude of low-frequency fluctuation (ALFF) for resting-state fMRI: fractional ALFF. Journal of neuroscience methods 172:137-141.
% see f_alff.m for authors, http://restfmri.net 
HighCutoff =0.08;			% Band Info for fALFF computing
LowCutoff =0.01;				% Band Info for fALFF computing
if ~exist('name4D')  
 name4D = spm_select(1,'image','Select 4D resting state session');
end;
[pth,nam,ext,vol] = spm_fileparts( name4D);
name4D = fullfile(pth,[ nam, ext]); %'img.nii,1' -> 'img.nii'
% load the 4D image
nii_hdr = spm_vol (name4D);
data4D = spm_read_vols (nii_hdr);
[nDim1 nDim2 nDim3 nDimTimePoints]=size(data4D);
% Get the frequency index
sampleFreq 	 = 1/kTRsec;
sampleLength = nDimTimePoints;
paddedLength = sampleLength; %no need to pad with Matlab
if (LowCutoff >= sampleFreq/2) % All high included
    idx_LowCutoff = paddedLength/2 + 1;
else % high cut off, such as freq > 0.01 Hz
    idx_LowCutoff = ceil(LowCutoff * paddedLength * kTRsec + 1);
    % Change from round to ceil: idx_LowCutoff = round(LowCutoff *paddedLength *ASamplePeriod + 1);
end
if (HighCutoff>=sampleFreq/2)||(HighCutoff==0) % All low pass
    idx_HighCutoff = paddedLength/2 + 1;
else % Low pass, such as freq < 0.08 Hz
    idx_HighCutoff = fix(HighCutoff *paddedLength *kTRsec + 1);
    % Change from round to fix: idx_HighCutoff	=round(HighCutoff *paddedLength *ASamplePeriod + 1);
end
data4D=reshape(data4D,[],nDimTimePoints)';
data4D = 2*abs(fft(data4D))/sampleLength;
fALFF_2D = sum(data4D(idx_LowCutoff:idx_HighCutoff,:)) ./ sum(data4D(2:(paddedLength/2 + 1),:));
fALFF_2D(~isfinite(fALFF_2D))=0;
fALFFBrain = reshape(fALFF_2D,nDim1, nDim2, nDim3);
%mask data
if ~exist('maskName') || isempty (maskName)
    fprintf('ALFF image is not masked\n');
else
    mask_hdr = spm_vol (maskName);
    mask_img = spm_read_vols (mask_hdr);
    %mask_img = (mask_img > 0);
    fALFFBrain(mask_img(:)==0)=0; %mask
end
%save data
outHdr = spm_vol (fullfile(pth,[ nam, ext, ',1']));
if ~exist('namePrefix')
    namePrefix = 'alf_';
end
outName = fullfile(pth,[namePrefix nam ext]);
outHdr.fname = outName;
outHdr.pinfo = [1;0;0];
outHdr.dt    =[16,0]; %32-bit real datatype
spm_write_vol(outHdr,fALFFBrain);
%end alffSub()