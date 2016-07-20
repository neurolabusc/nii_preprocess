function nii_rest (imgs, TRsec, SliceOrder, doSliceTime)
%preprocesses resting state data using SPM
%  imgs : (optional) structure of image names (imgs.T1, imgs.Rest imgs.SE, imgs.SErev)
%  TRsec      : Time between volumes in seconds
%  SliceOrder : see nii_sliceTime 0=auto-detect,1=ascend,2=descend
%Examples
% nii_rest; %use GUI
% imgs.T1 = 'T1.nii'; imgs.Rest = 'Rest.nii'; imgs.SE = 'RestSE.nii'; imgs.SErev = 'RestSErev.nii';
% nii_rest(imgs);
%
% History:
%
% TH 03/16: now pass in imgs structure rather than separate images for T1, rest, spin echo, and spin echo rev. 

if nargin < 5
    doSliceTime = true;
end
isSPM12orNewerSub;
%doSliceTime = true; %determine whether you want slice time correction 
doDeleteTemporary = false; %remove images generated during the middle stages of pipeline
resliceMM = 3; %reslicing resolution, e.g. if 3 then data will be 3x3x3 mm
if exist('spm','file')~=2; fprintf('%s requires SPM\n',mfilename); return; end;
if nargin < 1 %no input: select imgs[s]
 imgs.Rest = spm_select(inf,'image','Select 4D resting state volumes[s]');
 imgs.T1 = spm_select(1,'image','Select T1 anatomical scan');
 imgs.SE = spm_select(1,'image','(Optional) select spin-echo scan for undistortion');
 imgs.SE = spm_select(1,'image','(Optional) select reversed-phase spin-echo scan for undistortion');
end
Restnames = imgs.Rest;
T1name = imgs.T1;
if (~exist('TRsec','var')) || (TRsec == 0)
    Filename = deblank (Restnames(1,:));
    TRsec = getTRSub(Filename);
    if (TRsec ==0) 
            TRsec = str2double(cell2mat(inputdlg('Repitition time (TR, sec)?:', 'Timing', 1,{'2'})));
    end;
end
if ~exist('SliceOrder','var') || (SliceOrder == 0) %user did not specify slice order
    SliceOrder = getSliceOrderSub(Filename);
    if (SliceOrder == 0) %unable to auto-detect
        SliceOrder = str2double(cell2mat(inputdlg('Slice order (1=ascend,2=descend,3asc.int,4=desc.int,5=asc.int2,6=desc.int2)?', 'Timing', 1,{'1'})));
    end
end
fprintf('Assuming time between volumes (TR) is %.3f seconds\n', TRsec);
spm('Defaults','fMRI');
spm_jobman('initcfg'); % useful in SPM8 only
clear matlabbatch
tic; %start timer
% 0 - ONCE PER INDIVIDUAL normalize T1 scan using unified normalization-segmentation...
[defname, wc1, wc2, wc3] = segnormSub(T1name, resliceMM);   
%segnormwriteSub(defname,{mT1name,c1name, c2name, c3name}, resliceMM); %warp bias corrected image
%REMAINING STEPS ONCE PER RESTING STATE SESSION
nses = length(Restnames(:,1));
for ses = 1 : nses
    fMRIname = deblank (Restnames(ses,:));
    [pth,nam,ext] = spm_fileparts( deblank (fMRIname));
    prefix = ''; %prefixes applied to filename with operations, e.g. 'r'esliced, 's'moothed
    nVol = numVolumesSub(fMRIname);
    if nVol < 2 %this script uses 4d data
        fprintf('ERROR: %s requires 4D datasets volumes\n',mfilename);
        return;
    end;
    %1 adjust for slice time differences:
    if doSliceTime, %do we adjust for timing of 2D slices within 3D volume, uses spm_slice_timing
        prefix = [nii_sliceTime(fullfile(pth,[prefix, nam, ext]),SliceOrder, TRsec), prefix]; %#ok<AGROW>
    end;
    %2 motion correct data:
    prefixRealign = prefix;
    if isfield(imgs,'SE') && isfield(imgs,'SErev') && ~isempty(imgs.SE) && ~isempty(imgs.SErev)
        meanfmri = realignSub(fullfile(pth,[prefix, nam, ext]), true); %does not change prefix - only new matrix
        prefix = ['r' prefix]; %#ok<AGROW> %realigned data has prefix 'r'
        % nii_bold_undistort(prefix, 'fmri.nii', 'se.nii', 'seRev.nii');
        rfMRIname = prefixSub(prefix, fMRIname);
        nii_bold_undistort(rfMRIname, meanfmri, imgs.SE, imgs.SErev);
        prefix = ['u' prefix]; %#ok<AGROW> %realigned data has prefix 'r'
        % char /media/FAT500/update/M2120/meanaxfRest_AP_M2120_15.nii
        meanfmri = prefixSub('u', meanfmri);
        %fMRIname = prefixSub(prefix, fMRIname);
    else
        meanfmri = realignSub(fullfile(pth,[prefix, nam, ext])); %does not change prefix - only new matrix
    end
    %3 coregister fMRI to T1, then warp fMRI to standard space
    coregEstSub(T1name, meanfmri, prefix, fMRIname)
    segnormwriteSub(defname, meanfmri, resliceMM);  % spm_write_sn
    prefix = [segnormwriteSub(defname,  getsesvolsSub(fullfile(pth,[prefix, nam, ext])), resliceMM), prefix]; %#ok<AGROW>
    %4 generate a brain mask - must load from disk, may save to memory
    maskName = makeMaskSub(strvcat(wc1,wc2,wc3),0.05); %#ok<*REMFF1>
    %5 detrend for motion and white matter fluctuations AND SMOOTH
    prefix = [nii_detrend(fullfile(pth,[prefix, nam, ext]),maskName,wc2,fullfile(pth,['rp_' prefixRealign, nam, '.txt']),3,0.8,6), prefix]; %#ok<AGROW>
    % note: detrend PRIOR to bandpass http://www.ncbi.nlm.nih.gov/pubmed/23747457
    %6 Smooth data - no longer required (integrated with detrend)
    % prefix = [smoothSub(fullfile(pth,[prefix, nam, ext]), 6), prefix];   
    %7 do ALFF
    alffSub(fullfile(pth,[prefix, nam, ext]), TRsec, maskName);
    %8 remove frequencies
    prefix = [nii_temporalFilter(fullfile(pth,[prefix, nam, ext]), TRsec, true), prefix]; %#ok<AGROW>
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
%%END nii_rest()

%%%%%%% SUBFUNCTIONS FOLLOW

function isSPM12orNewerSub
%check that SPM is installed and is at least release 6225
if exist('spm','file') ~= 2, error('Please install SPM12 or later'); end;
[v,r] = spm('Ver','',1); r = str2double(r); %#ok<ASGLU>
if r < 6225, error('Please update your copy of SPM'); end;
%end isSPM12orNewer()

function [defname, wc1, wc2, wc3]  = segnormSub(t1, mm)
%apply new segment - return name of warping matrix
[pth,nam,ext] = spm_fileparts( deblank (t1));
defname = fullfile(pth, ['y_e', nam, ext]);
c1name = fullfile(pth, ['c1e', nam, ext]);
c2name = fullfile(pth, ['c2e', nam, ext]);
c3name = fullfile(pth, ['c3e', nam, ext]);
wc1 = fullfile(pth, ['wc1e', nam, ext]);
wc2 = fullfile(pth, ['wc2e', nam, ext]);
wc3 = fullfile(pth, ['wc3e', nam, ext]);
if ~exist(defname, 'file')
    defname = fullfile(pth, ['y_', nam, ext]);
    c1name = fullfile(pth, ['c1', nam, ext]);
    c2name = fullfile(pth, ['c2', nam, ext]);
    c3name = fullfile(pth, ['c3', nam, ext]);
    wc1 = fullfile(pth, ['wc1', nam, ext]);
    wc2 = fullfile(pth, ['wc2', nam, ext]);
    wc3 = fullfile(pth, ['wc3', nam, ext]);
end
if exist(defname, 'file') && exist(c1name,'file') && exist(c2name,'file') && exist(c3name,'file')
    fprintf('Skipping normalization (images exist) %s\n', defname);
    segnormwriteSub(defname, {c1name, c2name, c3name}, mm);
    return;
end
template = fullfile(spm('Dir'),'tpm','TPM.nii');
if ~exist(template,'file')
    error('Unable to find template named %s',template);
end
fprintf('NewSegment of %s\n', t1);
matlabbatch{1}.spm.spatial.preproc.channel(1).vols = {t1};
matlabbatch{1}.spm.spatial.preproc.channel(1).biasreg = 0.001;
matlabbatch{1}.spm.spatial.preproc.channel(1).biasfwhm = 60;
matlabbatch{1}.spm.spatial.preproc.channel(1).write = [0 1];
matlabbatch{1}.spm.spatial.preproc.tissue(1).tpm = {[template ',1']};
matlabbatch{1}.spm.spatial.preproc.tissue(1).ngaus = 1;
matlabbatch{1}.spm.spatial.preproc.tissue(1).native = [1 1];
matlabbatch{1}.spm.spatial.preproc.tissue(1).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(2).tpm = {[template ',2']};
matlabbatch{1}.spm.spatial.preproc.tissue(2).ngaus = 1;
matlabbatch{1}.spm.spatial.preproc.tissue(2).native = [1 1];
matlabbatch{1}.spm.spatial.preproc.tissue(2).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(3).tpm = {[template ',3']};
matlabbatch{1}.spm.spatial.preproc.tissue(3).ngaus = 2;
matlabbatch{1}.spm.spatial.preproc.tissue(3).native = [1 0];
matlabbatch{1}.spm.spatial.preproc.tissue(3).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(4).tpm = {[template ',4']};
matlabbatch{1}.spm.spatial.preproc.tissue(4).ngaus = 3;
matlabbatch{1}.spm.spatial.preproc.tissue(4).native = [1 0];
matlabbatch{1}.spm.spatial.preproc.tissue(4).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(5).tpm = {[template ',5']};
matlabbatch{1}.spm.spatial.preproc.tissue(5).ngaus = 4;
matlabbatch{1}.spm.spatial.preproc.tissue(5).native = [1 0];
matlabbatch{1}.spm.spatial.preproc.tissue(5).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(6).tpm = {[template ',6']};
matlabbatch{1}.spm.spatial.preproc.tissue(6).ngaus = 2;
matlabbatch{1}.spm.spatial.preproc.tissue(6).native = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(6).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.warp.mrf = 1;
matlabbatch{1}.spm.spatial.preproc.warp.cleanup = 1;
matlabbatch{1}.spm.spatial.preproc.warp.reg = [0 0.001 0.5 0.05 0.2];
matlabbatch{1}.spm.spatial.preproc.warp.affreg = 'mni';
matlabbatch{1}.spm.spatial.preproc.warp.fwhm = 0;
matlabbatch{1}.spm.spatial.preproc.warp.samp = 3;
matlabbatch{1}.spm.spatial.preproc.warp.write = [1 1];
spm_jobman('run',matlabbatch);
segnormwriteSub(defname, {c1name, c2name, c3name}, mm);
%end newSegSub()

function  prefix = segnormwriteSub(defname, targetname, mm, bb)
%reslice img using pre-existing new-segmentation deformation field
prefix = 'w';
if isempty(targetname) || isempty(defname), return; end;
if ~exist('bb','var'), bb = [-78 -112 -50; 78 76 85]; end; 
if ~exist('mm','var'), mm = 3; end;
if ~exist(defname,'file')
    error('Unable to find new-segment deformation image %s',defname);
end
%targetname = {strcat(char(targetname))};
if ischar(targetname)
    targetname = cellstr(targetname); %orig
    %targetname = {strcat(char(targetname))};
end
if size(targetname,2) > size(targetname,1) %more columns than rows
    %then transpose so that now more rows than columns
    targetname = targetname';
end
if ischar(defname)
    defname = cellstr(defname);
end
[pth,nam,ext] = spm_fileparts(targetname{1}); %orig
%[pth,nam,ext] = spm_fileparts(targetname{1}(1,:));
wtar = fullfile(pth, [prefix, nam, ext]); 
matlabbatch{1}.spm.spatial.normalise.write.subj.def = defname;
matlabbatch{1}.spm.spatial.normalise.write.subj.resample = targetname;
matlabbatch{1}.spm.spatial.normalise.write.woptions.bb = bb;
matlabbatch{1}.spm.spatial.normalise.write.woptions.vox = [mm mm mm];
matlabbatch{1}.spm.spatial.normalise.write.woptions.interp = 1; %4; //trilinear avoids ringing
if exist(wtar,'file')
    fprintf('Skipping write normalize (file exists) %s\n', wtar);
else
    spm_jobman('run',matlabbatch);
    if ~exist('binarize','var') || ~binarize, return; end;
    hdr = spm_vol(tar);
    img = spm_read_vols(hdr);
    mn = min(img(:));
    mx = max(img(:));
    thresh = ((mx-mn)*0.5) + mn;
    spm_write_vol(hdr,+(img > thresh));
end
%end newSegWriteSub()

function [nVol] = numVolumesSub (imgName)
%returns number of volumes (timepoints, directions) in a NIfTI image
[pth,nam,ext] = spm_fileparts( deblank (imgName));
imgName = fullfile(pth,[nam, ext]); %remove volume is 'img.nii,1' -> 'img.nii'
hdr = spm_vol(imgName); 
nVol = length(hdr);
%end numVolumesSub()

function slice_order = getSliceOrderSub(fMRIname)
[pth,nam,ext] = spm_fileparts( deblank (fMRIname));
fMRIname = fullfile(pth,[ nam, ext]); %remove volume, e.g. "fMRI.nii,1"-> "fMRI.nii"
fid = fopen(fMRIname);
fseek(fid,122,'bof');
slice_order = fread(fid,1,'uint8');
fclose(fid);
%end getSliceOrderSub()

function tr =  getTRSub(fMRIname)
%Returns Repeat Time in seconds for volume fMRIname
hdr = spm_vol(fMRIname);
if isfield(hdr(1,1).private.timing,'tspace')
  tr = hdr(1,1).private.timing.tspace;
else
  fprintf('%s error: unable to determine TR for image %s\n',mfilename,fmriname); 
end
%end getTRSub()

function meanImgOut = realignSub(inName, isReslice)
%if isReslice = true then resliced data is created with 'r' prefix
[pth,nam,ext] = spm_fileparts( deblank (inName));
meanImgOut = fullfile(pth,['mean', nam, ext]);
if ~exist(meanImgOut,'file')
    fprintf('Realigning data (motion correction), creating mean image named %s\n',meanImgOut);
    matlabbatch{1}.spm.spatial.realign.estwrite.data = {getsesvolsSub(inName)}';
    matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.quality = 0.9;
    matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.sep = 4;
    matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.fwhm = 5;
    matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.rtm = 1;
    matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.interp = 2;
    matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.wrap = [0 0 0];
    matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.weight = '';
    if exist('isReslice','var') && isReslice
       matlabbatch{1}.spm.spatial.realign.estwrite.roptions.which = [2 1]; 
    else
        matlabbatch{1}.spm.spatial.realign.estwrite.roptions.which = [0 1];
    end
    matlabbatch{1}.spm.spatial.realign.estwrite.roptions.which = [2 1];
    matlabbatch{1}.spm.spatial.realign.estwrite.roptions.interp = 4;
    matlabbatch{1}.spm.spatial.realign.estwrite.roptions.wrap = [0 0 0];
    matlabbatch{1}.spm.spatial.realign.estwrite.roptions.mask = 1;
    matlabbatch{1}.spm.spatial.realign.estwrite.roptions.prefix = 'r';
    spm_jobman('run',matlabbatch);
else
    disp(['File already exists: ' meanImgOut]);
end
%end  realignSub()

%function [imgNameVol1] = getVol1Sub(imgName)
% input: nifti image, output: name of FIRST volume
%  example 4D file with 3 volumes input= 'img.nii', output= 'img.nii,1'
%[pth,nam,ext] = spm_fileparts( deblank (imgName));
%imgNameVol1 = fullfile(pth,[nam, ext, ',1']);
% end getVol1Sub()

function coregEstSub(ref, imgToCoreg, prefix, otherImgs)
if ~exist('prefix','var'), prefix =''; end
if ~exist('otherImgs','var'), otherImgs = []; end
%coregister fmri data to match T1 image
fprintf('Coregistering %s to match %s\n',imgToCoreg,ref);
%fMRIses = getsesvolsSubHier(prefix, fmriname);
matlabbatch{1}.spm.spatial.coreg.estimate.ref = {ref};
matlabbatch{1}.spm.spatial.coreg.estimate.source = {imgToCoreg};
if isempty(otherImgs) 
    matlabbatch{1}.spm.spatial.coreg.estimate.other = {''};
else
    otherNam = prefixSub(prefix, otherImgs);
    if ~exist(otherNam,'file'), error('Coreg unable to find %s', otherNam); end;
    otherImgsSes = getsesvolsSub(otherNam);
    matlabbatch{1}.spm.spatial.coreg.estimate.other = otherImgsSes;
end;
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];
spm_jobman('run',matlabbatch);
%end coregEstSub()

function maskName = makeMaskSub(imgNames, thresh)
%creates a binary mask based on a series of images (e.g. white and gray matter tissue probability maps)
%   imgNames   : list of image(s) to sum together
%   thresh     : (optional) mask includes voxels that exceed this threshold, default = 0
%Example 
% nii_makeMask('wc1t1.nii',0.5);  
% nii_makeMask(strvcat('wc1t1.nii','wc2t1.nii'),0.2);
% nii_makeMask(strvcat('wc1t1.nii','wc2t1.nii','wc3t1.nii'),0.1);
if ~exist('imgNames','var') %no image specified
 imgNames = spm_select(inf,'image','Select images to binarize');
end;
if ~exist('thresh','var') %no image
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
outimg(inimg > thresh) = 1;
%tmp = find((inimg)> thresh);
%outimg(tmp) = 1;
hdr.fname = maskName;
spm_write_vol(hdr,outimg);
%end nii_makeMask()

function [sesvols] = getsesvolsSub(sesvol1)
% input: single volume from 4D volume, output: volume list
%  example 4D file with 3 volumes input= 'img.nii', output= {'img.nii,1';'img.nii,2';'img.nii,3'}
	[pth,nam,ext] = spm_fileparts( deblank (sesvol1));
	sesname = fullfile(pth,[nam, ext]);
	hdr = spm_vol(sesname);
	nvol = length(hdr);
	if (nvol < 2), fprintf('Error 4D fMRI data required %s\n', sesname); return; end;
    sesvols = cellstr([sesname,',1']);
    for vol = 2 : nvol
        sesvols = [sesvols; [sesname,',',int2str(vol)]  ]; %#ok<AGROW>
    end;
%end getsesvolsSub()

function outName = alffSub(name4D, kTRsec, maskName,namePrefix)
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
if ~exist('name4D','var')  
 name4D = spm_select(1,'image','Select 4D resting state session');
end;
[pth,nam,ext] = spm_fileparts( name4D);
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
if ~exist('maskName','var') || isempty (maskName)
    fprintf('ALFF image is not masked\n');
else
    mask_hdr = spm_vol (maskName);
    mask_img = spm_read_vols (mask_hdr);
    %mask_img = (mask_img > 0);
    fALFFBrain(mask_img(:)==0)=0; %mask
end
%save data
outHdr = spm_vol (fullfile(pth,[ nam, ext, ',1']));
if ~exist('namePrefix','var')
    namePrefix = 'alf_';
end
outName = fullfile(pth,[namePrefix nam ext]);
outHdr.fname = outName;
outHdr.pinfo = [1;0;0];
outHdr.dt    =[16,0]; %32-bit real datatype
spm_write_vol(outHdr,fALFFBrain);
%end alffSub()

function nam = prefixSub (pre, nam)
[p, n, x] = spm_fileparts(nam);
nam = fullfile(p, [pre, n, x]);
%end prefixSub()