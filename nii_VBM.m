function prefix = nii_VBM (imgs, TRsec, SliceOrder)
%preprocesses resting state data using SPM
%  imgs : structure of image names (imgs.T1, imgs.Rest)
%  TRsec      : Time between volumes in seconds
%  SliceOrder : see nii_sliceTime -1=no, 0=auto-detect,1=ascend,2=descend
%  doSliceTime : compute slice timing (true/false)
%Examples
% %linear Rest->T1, nonlinear T1->MNI: use when good linear fit between Rest and T1
%  imgs.T1 = 'T1_M2129_LIME.nii';
%  imgs.Rest = 'Rest_M2129_LIME copy.nii';
%  nii_rest(imgs);
% %nonlinear Rest->MNI: use when Rest and T1
%  imgs.T1 = [];
%  imgs.Rest = 'Rest_M2129_LIME copy.nii';
%  nii_rest(imgs);
% %linear Rest->T2, T2->T1, nonlinear T1->MNI, T2->wT1, nonlinear Rest->wT2
%  imgs.T2 = 'T2_M2129_LIME.nii'
%  imgs.T1 = 'T1_M2129_LIME.nii'
%  imgs.Rest = 'Rest.nii'
%  TRsec = 1.65;
%  nii_rest (imgs, TRsec);

%Roger Comment this in for Jill Data due to failed AutoDetect Slice Order of 5
%SliceOrder = 5;

if isempty(which('spm')) || ~strcmp(spm('Ver'),'SPM12'), error('SPM12 required'); end;
if ~exist('nii_batch12','file'), error('Make sure nii_batch12 is in path'); end;
if ~exist('imgs','var') %no input: select imgs[s]
    imgs.Rest = spm_select(inf,'image','Select 4D resting state volumes[s]');
    if isempty(imgs.Rest), return; end;
    imgs.T1 = spm_select(1,'image','(Optional) Select T1 anatomical scan');
    imgs.SE = spm_select(1,'image','(Optional) select spin-echo scan for undistortion');
    imgs.SE = spm_select(1,'image','(Optional) select reversed-phase spin-echo scan for undistortion');
    imgs.T2 = spm_select(1,'image','(Optional) Select T2 pathological scan');
end
if ~exist('TRsec','var') %no input: select imgs[s]
    TRsec = 0;
    fprintf('Will attempt to detect TR from image file.\n');
end
if ~exist('SliceOrder','var') %no input: select imgs[s]
    SliceOrder = 0;
    fprintf('Will attempt to detect Slice Order from image file.\n');
end
if ~isfield(imgs, 'T1'), error('T1 scan required'); end;
if ~isfield(imgs, 'Rest'), error('Rest scan required'); end;
if isfield(imgs, 'T2')
    p.t2name = imgs.T2;
else
    p.t2name = [];
end
clear matlabbatch
tic; %start timer
p.setOrigin = true;
p.t1name = imgs.T1; %'t1.nii'
p.TRsec = TRsec; %repeat time off 10 seconds
p.slice_order = SliceOrder;
p.FWHM = 6;
for ses = 1 : length(imgs.Rest(:,1));
    p.fmriname = deblank(imgs.Rest(ses,:));
    [pth,nam,ext] = spm_fileparts(p.fmriname);
    [prefix, TRsec, so, maskName] = nii_batch12(p);
    % warning('demo code: skipping preprocessing'); 
    % prefix = 'sw';
    % TRsec = 1.65;
    % maskName = 'wbmeanRest.nii';
    %ALFF - simple detrending only
    alffSub(fullfile(pth,[prefix, nam, ext]), TRsec, maskName, '', true);
    %we need white matter or CSF tissue map for advanced detrending
    if true
        wc2name = wc2(p.t1name); % use WM for detrending
    else
        disp ('using CSF for detrending...');
        wc2name = wc3(p.t1name); % use CSF for detrending %%% GY
    end
    if isempty(wc2name), continue; end; 
    % we need motion parameters for detrending
    mocoTxt = fullfile(pth,['rp_', nam, '.txt']);
    if ~exist(mocoTxt, 'file'), fprintf('Unable to find motion parameters %s\n', mocoTxt); continue; end;
    %Nonlinear detrend
    % prefix = [nii_detrend(fullfile(pth,[prefix, nam, ext]),maskName,wc2name,mocoTxt,3,0.8,6), prefix]; %#ok<AGROW>
    %we now smooth in the pre-processing step
    prefix = [nii_detrend(fullfile(pth,[prefix, nam, ext]),maskName,wc2name,mocoTxt,3,0.8,0), prefix]; %#ok<AGROW>
    %ALFF - after complex detrending
    alffSub(fullfile(pth,[prefix, nam, ext]), TRsec, maskName, '', true);
    %remove frequencies that are too high or low to be resting state
    prefix = [nii_temporalFilter(fullfile(pth,[prefix, nam, ext]), TRsec, true), prefix]; %#ok<AGROW>        
    % note: bandpass AFTER detrend http://www.ncbi.nlm.nih.gov/pubmed/23747457 
    % GY, Feb 3, 2016: running ICA and removing lesion-driven ICs
    if isfield (imgs, 'Lesion')
        if ~isempty (imgs.Lesion)
            [lpth, lnam, lext] = fileparts(imgs.Lesion);
            lesName = fullfile(lpth,['ws', lnam ,lext]);
            if ~exist(lesName,'file'), 
                lesName = fullfile(lpth,['wsr', lnam ,lext]); %resliced to match t1 
            end;
            if ~exist(lesName,'file'), 
                lesName = fullfile(lpth,['wsr', lnam ,lext]);   
            end;
            prefix = [nii_filter_lesion_ICs(lesName, fullfile(pth,[prefix, nam, ext]), TRsec), prefix]; %#ok<AGROW> 
        end
    end        
end; %for each session
fprintf('Done processing sessions in %0.3fsec\n', toc);
%end nii_rest()

function wc2name = wc2(t1name)
wc2name = [];
if isempty(t1name), fprintf('Skipping resting state analyses: no T1'); return; end;
wc2name = prefixSub ('wc2e', t1name);
if exist(wc2name, 'file'), return; end;
wc2name = prefixSub ('wc2', t1name);
if exist(wc2name, 'file'), return; end;
error('Unable to find white matter tissue map %s\n', wc2name);
%end wc2()

function wc2name = wc3(t1name)
wc2name = [];
if isempty(t1name), fprintf('Skipping resting state analyses: no T1'); return; end;
wc2name = prefixSub ('wc3e', t1name);
if exist(wc2name, 'file'), return; end;
wc2name = prefixSub ('wc3', t1name);
if exist(wc2name, 'file'), return; end;
error('Unable to find CSF map %s\n', wc2name);
%end wc3()

function outName = alffSub(name4D, kTRsec, maskName,namePrefix, isDetrend, divideByMask)
%Identify low frequency flutuations 
% name4D     : filename of 4D NIFTI image
% TRsec      : Repeat time (seconds)
% nameMask   : (optional) name for masking image
% namePrefix : (optional) text appended to output image name, default 'alf_'
% isDetrend  : (optional) apply linear detrend
% divideByMask : (optional) pALFF typically normalizes intensity by dividing by all voxels in brain mask, set to 'false' to save raw data 
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
if isDetrend
    data4D = detrendSub(data4D);
    if false %following is a test routine to examine de-trending performance
        outHdr = nii_hdr(1);
        outName = fullfile(pth,['detrend_', nam, ext]);
        outHdr.fname = outName;
        outHdr.pinfo = [1;0;0];
        outHdr.dt    =[16,0]; %32-bit real datatype
        for vol = 1 : size(data4D,4)
            outHdr.n(1)=vol;
            spm_write_vol(outHdr,data4D(:, :, :, vol));
        end;
    end
end
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
outHdr = nii_hdr(1);
if ~exist('namePrefix','var') || isempty(namePrefix)
    namePrefix = 'alf_';
end
outHdr.pinfo = [1;0;0];
outHdr.dt    =[16,0]; %32-bit real datatype
if exist('divideByMask','var') && (~divideByMask)
    outName = fullfile(pth,[namePrefix nam ext]);
    outHdr.fname = outName;
    spm_write_vol(outHdr,fALFFBrain);
    return;
end
%save fractional
fALFFBrain = divideMeanWithinMaskSub(fALFFBrain, mask_img);
outName = fullfile(pth,['p' namePrefix nam ext]);
outHdr.fname = outName;
spm_write_vol(outHdr,fALFFBrain);
%end alffSub()

function img = detrendSub(img)
%remove linear trends from data
dim = size(img);
if (dim(4) < 2), error('Not a 4D image %s', fnm); end;
img = reshape(img, prod(dim(1:3)), dim(4))';
ave = mean(img);
img = detrend(img);
img = img';
img = bsxfun(@plus,img,ave');
img = reshape(img, dim(1), dim(2), dim(3), dim(4));
%end detrendSub()

function img = divideMeanWithinMaskSub(img,msk)
msk(isnan(msk)) = 0;
ok = img(msk ~= 0);
ave = mean(ok(:));
img = img/ave;
img(msk == 0) = 0;
%end divideMeanWithinMaskSub()

function nam = prefixSub (pre, nam)
[p, n, x] = spm_fileparts(nam);
nam = fullfile(p, [pre, n, x]);
%end prefixSub()

