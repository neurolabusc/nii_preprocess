function [prefix, TRsec, slice_order, meanname] = nii_batch12 (p)
%preprocess and analyze fMRI data using standard settings
% p
%   structure for preprocessing
%  p.fmriname : name of 4D fMRI volumes
%  p.fmriname, t1name, TRsec, slice_order, phase, magn
%  p.t1name   : name of anatomical scan (-1 to skip)
%  p.TRsec    : TR for fMRI data, 0=auto
%  p.slice_order : EPI slice order 0=auto,1=[1234],2=[4321],3=[1324],4=[4231], 5=[2413],6=[3142]
%  p.phase  (optional) name of fieldmap phase image
%  p.magn : (optional) name of fieldmap magnitude image
%   (optional) structure for statistical analysis
%  p.names : condition names
%  p.onset : condition onset times
%  p.duration : duration of events (either a single value or array matching p.onset)
%  p.mocoRegress : true or false: should motion parameters be modeled?
%  p.t2name   : name of anatomical scan (typically not provided)
%  p.FWHM : full-idth half maximum for smoothing (defaults to 6mm)
%  p.resliceMM : resampled resolution (defaults to 2mm)
%Versionsf
%   7/7/2016 added "normIntensitySub" to deal with Philips data
%Examples
% TO DO

[fmriname, t1name, TRsec, slice_order, phase, magn, prefix, isNormalized, t2name, FWHM, resliceMM] = validatePreprocSub(p);
deleteImagesSub('swa', fmriname, true); %delete previous images & mat files
deleteImagesSub('wbmean', fmriname, true); %delete previous images & mat files
if isempty(spm_figure('FindWin','Graphics')), spm fmri; end; %launch SPM if it is not running
spm_jobman('initcfg'); % useful in SPM8 only
%0.) set origin fMRI
if isfield(p,'setOrigin') && (p.setOrigin == true)
    setOriginSub(strvcat(fmriname, phase, magn), 3);   %#ok<DSTRVCT> %align images to MNI space 
end
%0.) set origin T1
if ~isempty(t1name) && (~isNormalized) && isfield(p,'setOrigin') && (p.setOrigin == true)
       setOriginSub(t1name, 1); % align images to MNI space
end
%1.) motion correct, used fieldmap if specifiedclass
if isfield(p,'SE') && isfield(p,'SErev') && ~isempty(p.SE) && ~isempty(p.SErev)
    error('Spin echo undistortion not yet integrated!');
else
    [meanname, prefix] = mocoFMSub(prefix, fmriname, phase, magn);
end;
%2.) brain extact mean (for better coregistration)
%if ~isNormalized %bet can fail with patient scans
if ~isempty(t1name) %if we have a T1, coregister BET fMRI to BET T1, if we do not have T1 we will keep the scalp (it is in TPM.nii)
   meanname = betSubFSL(meanname); %brain extract mean image for better coreg
    %meanname = betSub(meanname); %brain extract mean image for better coreg
end
%end
%3.) slice-time correction
prefix = slicetimeSub(prefix, fmriname, TRsec, slice_order); %slice-time correct
%4.) estimate normalization, coregister and reslice fMRI
if isNormalized %if the image was previously normalized use existing paramters
    if isempty(t2name)
        prefix = normWriteSub(t1name, meanname, prefix, fmriname, resliceMM);
    else
        prefix = normT2Sub(t1name, t2name, meanname, prefix, fmriname, resliceMM);
    end
else
    prefix = normNewSegSub(t1name, meanname, prefix, fmriname, resliceMM);
    
end
normWriteWhiteMatterSub(t1name, resliceMM); %we want a wc2* white matter prob map for detrending
normWriteCsfSub(t1name, resliceMM); %GY: we want a wc3* CSF prob map for detrending

meanname = prefixSub('w', meanname);
%5.) blur data
prefix = smoothSub(FWHM, prefix, fmriname); %smooth images
%-- get rid of images we don't need
deleteImagesSub(prefix, fmriname); %delete intermediate images

%6.) (GY, Oct 2017): detrend 
% maybe remove this thing from nii_rest because we do it in batch12 anyway?
if isfield (p, 'lesion_name') &&  ~isempty (p.lesion_name) % i.e. we preprocess task fmri, not rest
    [pth,nam,ext] = spm_fileparts(fmriname);
    mocoTxt = fullfile(pth,['rp_', nam, '.txt']);    
    prefix = [nii_detrend(fullfile(pth,[prefix, nam, ext]),meanname,'',mocoTxt,1), prefix]; %#ok<AGROW>    
    % in future: add WM regression??? -GY
end
%7.) (GY, Oct 5, 2017): running ICA and removing lesion-driven ICs
if isfield (p, 'lesion_name')
    if ~isempty (p.lesion_name)
        [lpth, lnam, lext] = fileparts(p.lesion_name);
        lesName = fullfile(lpth,['wsr', lnam ,lext]);
        if ~exist(lesName,'file')
            lesName = fullfile(lpth,['ws', lnam ,lext]);
        end
        if ~exist(lesName,'file'), error('Catastrophic failure: can not find %s', lesName); end;
        [pth, nam, ext] = fileparts(fmriname);
        prefix = [nii_filter_lesion_ICs(lesName, fullfile(pth,[prefix, nam, ext]), p.TRsec), prefix]; %#ok<AGROW> 
    end
end

%8.) compute statistics
if ~isfield(p,'onsets'), return; end; %only if user provides details
stat_1st_levelSub (prefix, fmriname, TRsec, p);
%end nii_batch()




%---------- LOCAL FUNCTIONS FOLLOW

function prefix = normT2Sub(t1, t2, meanname, prefix, fmriname, resliceMM)
t2 = maskSub(t2,  t1);
newSegWriteSub(t1, t2, '', resliceMM);
t2 = prefixSub('w', t2); %we will match fMRI scan to normalized T2
coregEstSub(t2, meanname, prefix, fmriname); %make sure fMRI is aligned with T2
prefix = oldNormSub(t2, meanname, fmriname, prefix);
%end normT2Sub()

function mskfnm = maskSub(fnm, t1)
%zero voxels in fnm that are zero in fnmMask
pad = 'e'; %enatiomorphic
efnm = prefixSub(['c1', pad], t1);
if ~exist(efnm, 'file')
   pad = ''; %no enantiomorphic
   fnm = prefixSub(['c1', pad], t1);
   if ~exist(fnm, 'file')
    error('Unable to find gray matter estimate "%s" or "%s"\n', fnm, efnm);
   end
end
hdr = spm_vol(prefixSub(['c1', pad], t1)); %GM
imgM = spm_read_vols(hdr);
hdr = spm_vol(prefixSub(['c2', pad], t1)); %WM
imgM = imgM + spm_read_vols(hdr);
hdr = spm_vol(prefixSub(['c3', pad], t1)); %CSF
imgM = imgM + spm_read_vols(hdr);
img(imgM < 0.05) = 0; %zero voxels that are less than 5% GM+WM+CSF
hdr = spm_vol(deblank(fnm));
img = spm_read_vols(hdr);
img(imgM == 0) = 0;
mskfnm = prefixSub('b', fnm);
hdr.fname = mskfnm;  
spm_write_vol(hdr,img);
%end maskSub()

function prefix = oldNormSub(template, meanfmri, fmri, prefix)
warpses = getsesvolsSubFlat(prefix, fmri);
mbatch{1}.spm.tools.oldnorm.estwrite.subj.source = {meanfmri};
mbatch{1}.spm.tools.oldnorm.estwrite.subj.wtsrc = '';
%mbatch{1}.spm.tools.oldnorm.estwrite.subj.resample = {meanfmri};
mbatch{1}.spm.tools.oldnorm.estwrite.subj.resample = [{meanfmri}; warpses];
mbatch{1}.spm.tools.oldnorm.estwrite.eoptions.template = {template};
%mbatch{1}.spm.tools.oldnorm.estwrite.eoptions.weight = ''; %TODO: CHECK THIS HELPS
%mbatch{1}.spm.tools.oldnorm.estwrite.eoptions.weight = {fullfile(spm('Dir'),'toolbox','FieldMap','brainmask.nii')};
btemplate = extract2Sub(template);
mbatch{1}.spm.tools.oldnorm.estwrite.eoptions.weight = {btemplate}; 

mbatch{1}.spm.tools.oldnorm.estwrite.eoptions.smosrc = 4;
mbatch{1}.spm.tools.oldnorm.estwrite.eoptions.smoref = 4; %blur the template, too!
mbatch{1}.spm.tools.oldnorm.estwrite.eoptions.regtype = 'mni';
mbatch{1}.spm.tools.oldnorm.estwrite.eoptions.cutoff = 25;
mbatch{1}.spm.tools.oldnorm.estwrite.eoptions.nits = 16;
mbatch{1}.spm.tools.oldnorm.estwrite.eoptions.reg = 1;
mbatch{1}.spm.tools.oldnorm.estwrite.roptions.preserve = 0;
mbatch{1}.spm.tools.oldnorm.estwrite.roptions.bb = [NaN NaN NaN; NaN NaN NaN];
mbatch{1}.spm.tools.oldnorm.estwrite.roptions.vox = [NaN NaN NaN];
mbatch{1}.spm.tools.oldnorm.estwrite.roptions.interp = 1;
mbatch{1}.spm.tools.oldnorm.estwrite.roptions.wrap = [0 0 0];
mbatch{1}.spm.tools.oldnorm.estwrite.roptions.prefix = 'w';
spm_jobman('run',mbatch);
prefix = ['w' prefix];
%end oldNormSub()

function fnm = extract2Sub(fnm)
%subroutine to extract brain from surrounding scalp

fprintf('Brain extraction of %s\n', fnm);
[pth,nam,ext] = spm_fileparts(fnm);
%load headers
hdr = spm_vol(fnm);%bias corrected T1
img = spm_read_vols(hdr);
img(img ~= 0) = 1;
fnm = fullfile(pth,['b',  nam, ext]);
hdr.fname = fnm;
hdr.pinfo = [0; 0; 352];
hdr.dt(2) = 2; %16-bit precision more than sufficient uint8=2; int16=4; int32=8; float32=16; float64=64
spm_write_vol(hdr,img);
%end extractSub()

function [meanname, prefix] = mocoFMSub(prefix, fmriname, phase, magn) %motion correct with field map
if isempty(phase) || isempty(magn)
    meanname = mocoSub(prefix, fmriname); %motion correct
    return;
end;
FieldMapSub(fmriname, phase, magn);
FieldMapMocoSub(fmriname,phase);
prefix = strcat('u',prefix); %'u'nwarped
meanname = prefixSub(['mean', prefix], fmriname);
%end mocoFMSub()

function FieldMapMocoSub(fmriname,phase)
%spm('Defaults','fMRI');
%spm_jobman('initcfg');
%clear matlabbatch
fprintf('Motion correction and fieldmap unwarping\n');
nsessions = length(fmriname(:,1));
for s = 1 : nsessions
    [pth,nam,ext,~] = spm_fileparts( deblank (fmriname(s,:)));
    sesname = addpth(fullfile(pth,[nam, ext]));
    fMRIses = getsesvolsSub(sesname);
    matlabbatch{1}.spm.spatial.realignunwarp.data(s).scans = fMRIses;
    %get fieldmap voxel displacement name
    [pth,nam,ext,~] = spm_fileparts(phase);
    if length(pth) < 1; pth = pwd;end;
    if nsessions < 2
        longname = fullfile(pth,['vdm5_sc', nam,  ext]);
    else
        longname = fullfile(pth,['vdm5_sc', nam, '_z',  int2str(s), ext]);
    end;
    if exist(longname,'file')~=2; fprintf('Warning: unable to find fieldmap voxel displacement map %s - cd to approrpiate working directory?\n',longname); end;
    longname = [longname, ',1']; %#ok<AGROW>
    matlabbatch{1}.spm.spatial.realignunwarp.data(s).pmscan = {longname};
end;
matlabbatch{1}.spm.spatial.realignunwarp.eoptions.quality = 0.9;
matlabbatch{1}.spm.spatial.realignunwarp.eoptions.sep = 4;
matlabbatch{1}.spm.spatial.realignunwarp.eoptions.fwhm = 5;
matlabbatch{1}.spm.spatial.realignunwarp.eoptions.rtm = 0;
matlabbatch{1}.spm.spatial.realignunwarp.eoptions.einterp = 2;
matlabbatch{1}.spm.spatial.realignunwarp.eoptions.ewrap = [0 0 0];
matlabbatch{1}.spm.spatial.realignunwarp.eoptions.weight = '';
matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.basfcn = [12 12];
matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.regorder = 1;
matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.lambda = 100000;
matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.jm = 0;
matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.fot = [4 5];
matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.sot = [];
matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.uwfwhm = 4;
matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.rem = 1;
matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.noi = 5;
matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.expround = 'Average';
matlabbatch{1}.spm.spatial.realignunwarp.uwroptions.uwwhich = [2 1];
matlabbatch{1}.spm.spatial.realignunwarp.uwroptions.rinterp = 4;
matlabbatch{1}.spm.spatial.realignunwarp.uwroptions.wrap = [0 0 0];
matlabbatch{1}.spm.spatial.realignunwarp.uwroptions.mask = 1;
matlabbatch{1}.spm.spatial.realignunwarp.uwroptions.prefix = 'u';
spm_jobman('run',matlabbatch);
%end FieldMapMocoSub()

function FieldMapSub(fmriname, phase, magnitude)
%Input FMRI scans, and fieldmaps (both phase and magnitude images)
%  computes undistortion maps
%Examples
% nii_fieldmap(strvcat('fMRIrun1.nii','fMRIrun2.nii'),'phase.nii','mag.nii')
% nii_fieldmap('fMRI.nii','phase.nii','mag.nii');
%spm('Defaults','fMRI');
%spm_jobman('initcfg');
%clear matlabbatch
%typical settings to change:
te1 = 5.19; %short echo time for fieldmap, ms
te2 = 7.65; %long echo time for fieldmap, ms
frmi_readout = 0.49 * 34; %EPI readout time for fMRI scan
%Typical readout times at MCBI
%  no iPAT = 35.2ms (echo spacing of 0.55ms * 64 lines)
%  x2 iPAT = 18.24ms (echo spacing of 0.57ms * 32 lines)
blip_dir = -1; %readout direction: -1 or +1
fprintf('Fieldmap assumptions (from %s)\n',which(mfilename));
fprintf (' Short/Long TEs %0.2f/%0.2fms, EPI readout time %0.2fms, Blip direction %d\n',  te1, te2, frmi_readout, blip_dir);
%run with the selected settings...
magnitude = addpth (magnitude);
phase = addpth(phase);
matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj.session.epi = {fmriname};
matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj.magnitude = {magnitude};
matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj.phase = {phase};
matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj.anat = {magnitude};
nsessions = length(fmriname(:,1));
for s = 1 : nsessions
    sesname = addpth( deblank (fmriname(s,:)));
    matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj.session(s).epi = {sesname};
end;
matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj.defaults.defaultsval.et = [te1 te2];
matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj.defaults.defaultsval.maskbrain = 1;
matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj.defaults.defaultsval.blipdir = blip_dir;
matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj.defaults.defaultsval.tert = frmi_readout;
matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj.defaults.defaultsval.epifm = 0;
matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj.defaults.defaultsval.ajm = 0;
matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj.defaults.defaultsval.uflags.method = 'Mark3D';
matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj.defaults.defaultsval.uflags.fwhm = 10;
matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj.defaults.defaultsval.uflags.pad = 0;
matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj.defaults.defaultsval.uflags.ws = 1;
template = fullfile(spm('Dir'),'canonical','avg152T1.nii');
if ~exist(template,'file')
    error('Unable to find fieldmap template image %s',template);
end
matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj.defaults.defaultsval.mflags.template = {template};
matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj.defaults.defaultsval.mflags.fwhm = 5;
matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj.defaults.defaultsval.mflags.nerode = 2;
matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj.defaults.defaultsval.mflags.ndilate = 4;
matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj.defaults.defaultsval.mflags.thresh = 0.5;
matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj.defaults.defaultsval.mflags.reg = 0.02;
matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj.matchvdm = 1;
matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj.sessname = 'z';
matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj.writeunwarped = 1;
matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj.matchanat = 0;
spm_jobman('run',matlabbatch);
%end FieldMapSub()

function isSPM12orNewerSub
%check that SPM is installed and is at least release 6225
if exist('spm','file') ~= 2, error('Please install SPM12 or later'); end;
[v,r] = spm('Ver','',1); r = str2double(r); %#ok<ASGLU>
if r < 6225, error('Please update your copy of SPM'); end;
%end isSPM12orNewer()

function prefix = normNewSegSub(t1, meanname, prefix, fmriname, resliceMM)
%warp data to standard space with new segment
% t1 : filename of T1-weighted image
% meanname : name of mean fMRI data
% prefix   : prefix appended to name of fMRI data
% fmriname : base name(s) of 4D fMRI images
%Example
% normNewSeg('T1.nii','meanfmri.nii','fmri.nii','');
% normNewSeg('T1.nii','meanfmri.nii','fmri.nii','a'); %normalize 'afmri.nii'
if ~exist('meanname','var'), meanname = ''; end;
if ~exist('prefix','var'), prefix = ''; end;
if ~exist('fmriname','var'), fmriname = ''; end;
if ~exist('resliceMM','var'), resliceMM = 2; end;
if isempty(t1)
   prefix = normSub( meanname, prefix, fmriname, resliceMM);
   return;
end
newSegSub(t1); %normalize images
extractSub(0.01, t1, prefixSub('c1', t1), prefixSub('c2', t1));
prefix = normWriteSub(t1, meanname, prefix, fmriname, resliceMM);
%end normNewSeg()

function prefix = normSub( meanname, prefix, fmriname, resliceMM)
mbatch{1}.spm.spatial.normalise.estwrite.subj.vol = {meanname};
warpses = getsesvolsSubFlat(prefix, fmriname);
warpses = [warpses; {meanname}];
fprintf('Normalizing %s to match TPM\n', meanname);
mbatch{1}.spm.spatial.normalise.estwrite.subj.resample = warpses;
mbatch{1}.spm.spatial.normalise.estwrite.eoptions.biasreg = 0.0001;
mbatch{1}.spm.spatial.normalise.estwrite.eoptions.biasfwhm = 60;
template = fullfile(spm('Dir'),'tpm','TPM.nii');
if ~exist(template,'file')
    error('Unable to find template named %s',template);
end
mbatch{1}.spm.spatial.normalise.estwrite.eoptions.tpm = {template};
mbatch{1}.spm.spatial.normalise.estwrite.eoptions.affreg = 'mni';
mbatch{1}.spm.spatial.normalise.estwrite.eoptions.reg = [0 0.001 0.5 0.05 0.2];
mbatch{1}.spm.spatial.normalise.estwrite.eoptions.fwhm = 0;
mbatch{1}.spm.spatial.normalise.estwrite.eoptions.samp = 3;
mbatch{1}.spm.spatial.normalise.estwrite.woptions.bb = [-78 -112 -70; 78 76 85];
mbatch{1}.spm.spatial.normalise.estwrite.woptions.vox = [resliceMM resliceMM resliceMM];
mbatch{1}.spm.spatial.normalise.estwrite.woptions.interp = 4;
mbatch{1}.spm.spatial.normalise.estwrite.woptions.prefix = 'w'; %use spm_update if you get an error regarding woptions prefix
spm_jobman('run',mbatch);
prefix = ['w' prefix];
%end normSub()

function prefix = normWriteSub(t1, meanname, prefix, fmriname, resliceMM)
coregEstSub(prefixSub('b', t1), meanname, prefix, fmriname); %make sure fMRI is aligned with T1
newSegWriteSub(t1, t1, '', 0.9); %reslice anatomical with 0.9mm isotropic
newSegWriteSub(t1, prefixSub('b', t1), '', 0.9); %reslice anatomical with 0.9mm isotropic
newSegWriteSub(t1, meanname, '', resliceMM);
prefix = newSegWriteSub(t1, fmriname, prefix, resliceMM);
% end normWriteSub()

function nam = prefixSub (pre, nam)
[p, n, x] = spm_fileparts(nam);
nam = fullfile(p, [pre, n, x]);
%end prefixSub()

function t1Bet = extractSub(thresh, t1, c1, c2, c3)
%subroutine to extract brain from surrounding scalp
% t1: anatomical scan to be extracted
% c1: gray matter map
% c2: white matter map
% c3: [optional] spinal fluid map
fprintf('Brain extraction of %s\n', t1);
[pth,nam,ext] = spm_fileparts(t1);
%load headers
mi = spm_vol(t1);%bias corrected T1
gi = spm_vol(c1);%Gray Matter map
wi = spm_vol(c2);%White Matter map
%load images
m = spm_read_vols(mi);
g = spm_read_vols(gi);
w = spm_read_vols(wi);
if nargin > 4 && ~isempty(c3)
   ci = spm_vol(c3);%CSF map
   c = spm_read_vols(ci);
   w = c+w;
end;
w = g+w;
if thresh <= 0
    m=m.*w;
else
    mask= zeros(size(m));
    for px=1:length(w(:)),
      if w(px) >= thresh
        mask(px) = 255;
      end;
    end;
    spm_smooth(mask,mask,1); %feather the edges
    mask = mask / 255;
    m=m.*mask;
end;
mi.fname = fullfile(pth,['b',  nam, ext]);
mi.dt(1) = 4; %16-bit precision more than sufficient uint8=2; int16=4; int32=8; float32=16; float64=64
spm_write_vol(mi,m);
t1Bet = mi.fname;
%end extractSub()

function normWriteWhiteMatterSub(t1name, resliceMM) %we want a wc2* white matter prob map for detrending
if isempty(t1name), return; end;
eChar = 'e';
[pth,nam,ext] = spm_fileparts(t1name);
c2name = fullfile(pth, ['c2', eChar, nam, ext]);
if ~exist(c2name, 'file')
	eChar = '';
	c2name = fullfile(pth, ['c2', eChar, nam, ext]);
    if ~exist(c2name, 'file')
        fprintf('Unable to find white matter tissue map %s\n', c2name);
        return;
    end
end 
newSegWriteSub(t1name, c2name, '', resliceMM);
%end normWriteWhiteMatterSub()

function normWriteCsfSub(t1name, resliceMM) %GY: we also want a wc3* CSF prob map for detrending
if isempty(t1name), return; end;
eChar = 'e';
[pth,nam,ext] = spm_fileparts(t1name);
c3name = fullfile(pth, ['c3', eChar, nam, ext]);
if ~exist(c3name, 'file')
	eChar = '';
	c3name = fullfile(pth, ['c3', eChar, nam, ext]);
    if ~exist(c3name, 'file')
        fprintf('Unable to find CSF tissue map %s\n', c3name);
        return;
    end
end 
newSegWriteSub(t1name, c3name, '', resliceMM);
%end normWriteWhiteMatterSub()


function  prefix = newSegWriteSub(t1name, warpname, prefix, resliceMM, bb)
%reslice img using pre-existing new-segmentation deformation field
if isempty(warpname) || isempty(t1name), return; end;
[pth,nam,ext] = spm_fileparts(t1name);
defname = fullfile(pth,['y_' nam ext]);
if ~exist(defname,'file')
    defname = fullfile(pth,['y_e' nam ext]);
    if ~exist(defname,'file')
        defname2 = fullfile(pth,['y_' nam ext]);
        error('Unable to find new-segment deformation image %s or %s',defname2, defname);
    end
end
warpses = getsesvolsSubFlat(prefix, warpname);
matlabbatch{1}.spm.spatial.normalise.write.subj.def = {defname};
matlabbatch{1}.spm.spatial.normalise.write.subj.resample = warpses;
if ~exist('bb','var')
	matlabbatch{1}.spm.spatial.normalise.write.woptions.bb = [-78 -112 -70; 78 76 85];
else
	matlabbatch{1}.spm.spatial.normalise.write.woptions.bb = bb;
end
matlabbatch{1}.spm.spatial.normalise.write.woptions.vox = [resliceMM resliceMM resliceMM];
matlabbatch{1}.spm.spatial.normalise.write.woptions.interp = 4;
spm_jobman('run',matlabbatch);
prefix = ['w' prefix];
%end newsegwritesub()

function newSegSub(t1, t2)
%apply new segment - return name of warping matrix
template = fullfile(spm('Dir'),'tpm','TPM.nii');
if ~exist(template,'file')
    error('Unable to find template named %s',template);
end
fprintf('NewSegment of %s\n', t1);
matlabbatch{1}.spm.spatial.preproc.channel(1).vols = {t1};
matlabbatch{1}.spm.spatial.preproc.channel(1).biasreg = 0.001;
matlabbatch{1}.spm.spatial.preproc.channel(1).biasfwhm = 60;
matlabbatch{1}.spm.spatial.preproc.channel(1).write = [1 1];
if nargin > 1 && ~isempty(t2)
    matlabbatch{1}.spm.spatial.preproc.channel(2).vols = {t2};
    matlabbatch{1}.spm.spatial.preproc.channel(2).biasreg = 0.0001;
    matlabbatch{1}.spm.spatial.preproc.channel(2).biasfwhm = 60;
    matlabbatch{1}.spm.spatial.preproc.channel(2).write = [0 0];
end;
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
%end newSegSub()

function deleteImagesSub(prefix, fmriname, isDeleteMat)
%delete intermediate images, e.g. if prefix is 'swa' then delete 'wa' and 'a' images
if length(prefix) < 2, return; end;
for s = 1 : length(fmriname(:,1))
    for i = 2 : length(prefix)
        [pth,nam,ext] = spm_fileparts( deblank (fmriname(s,:)));
        nam = fullfile(pth,[prefix(i:length(prefix)), nam, ext]);
        if exist(nam, 'file')
            fprintf('Deleting  %s\n',nam );
            delete(nam);
        end;
        if exist('isDeleteMat','var') && isDeleteMat
            nam = fullfile(pth,['rp_', prefix(i:length(prefix)), nam, '.txt']);
            if exist(nam, 'file')
                fprintf('Deleting  %s\n',nam );
                delete(nam);
            end;
            nam = fullfile(pth,[prefix(i:length(prefix)), nam, '.mat']);
            if exist(nam, 'file')
                fprintf('Deleting  %s\n',nam );
                delete(nam);
            end;
        end
    end;
end;
%end deleteImagesSub()

function prefix = smoothSub(FWHMmm, prefix, fmriname)
%blur images
if FWHMmm <= 0, return; end;
fprintf('Smoothing with a %gmm FWHM Gaussian kernel\n',FWHMmm);
matlabbatch{1}.spm.spatial.smooth.data = getsesvolsSubFlat(prefix, fmriname);
matlabbatch{1}.spm.spatial.smooth.fwhm = [FWHMmm FWHMmm FWHMmm];
matlabbatch{1}.spm.spatial.smooth.dtype = 0;
matlabbatch{1}.spm.spatial.smooth.im = 0;
matlabbatch{1}.spm.spatial.smooth.prefix = 's';
spm_jobman('run',matlabbatch);
prefix = ['s' prefix];
%smoothSub()

function coregEstSub(t1, meanfmri, prefix, fmriname)
%coregister fmri data to match T1 image
fprintf('Coregistering %s to match %s\n',meanfmri,t1);

%fMRIses = getsesvolsSubHier(prefix, fmriname);
fMRIses = getsesvolsSubFlat(prefix, fmriname);
matlabbatch{1}.spm.spatial.coreg.estimate.ref = {t1};
matlabbatch{1}.spm.spatial.coreg.estimate.source = {meanfmri};
matlabbatch{1}.spm.spatial.coreg.estimate.other = fMRIses;
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];
spm_jobman('run',matlabbatch);
%save('coreg.mat', 'matlabbatch')
%end coregEstSub()

function [fmriname, t1name, TRsec, slice_order, phase, magn, prefix, isNormalized, t2name, FWHM, resliceMM] = validatePreprocSub(p)
%check all inputs
if exist('spm','file')~=2; fprintf('%s requires SPM\n',which(mfilename)); return; end;
isSPM12orNewerSub; %check recent SPM
prefix = ''; %originally '', if 'swa' it means we have smoothed, warped, aligned data
if ~isfield(p,'fmriname'), error('%s requires field "fmriname"',mfilename); end;
if ~isfield(p,'t1name'), error('%s requires field "t1name"',mfilename); end;
if ~isfield(p,'TRsec'), p.TRsec = 0; end;
if ~isfield(p,'slice_order'), p.slice_order = 0; end;
if ~isfield(p,'phase'), p.phase = ''; end;
if ~isfield(p,'magn'), p.magn = ''; end;
if isfield(p,'resliceMM') && isnumeric(p.resliceMM)
    resliceMM = p.resliceMM;
else
    resliceMM = 2; %resolution for reslicing data
end;
if isfield(p,'FWHM') && isnumeric(p.FWHM)
    FWHM = p.FWHM;
else
    FWHM = 6; 

end;
if ~isfield(p,'t2name') || isempty(p.t2name), 
    t2name = [];
else
    t2name = findImgSub(p.t2name, '');
    if ~exist(t2name, 'file')
        fprintf('T2 will not be used for fMRI normalization (unable to find registered %s\n', t2name);
        t2name = '';
    else
        disp(t2name)
        [pth, nm, ext] = fileparts(t2name);
        t2name = fullfile(pth, ['r', nm, ext]);
    end;
end;

isNormalized = false;
if (~ischar(p.t1name)) %we use p.t1name = -1 to signify skipping T1
    t1name = [];
else
    if isempty(p.t1name)
        t1name = spm_select(1,'image','Select T1 image volume');
        isNormalized = isNormalizedExists(t1name);
    end
    t1name = findImgSub(p.t1name, '');
    isNormalized = isNormalizedExists(t1name);
end;
p.fmriname = findImgSub(p.fmriname, p.t1name);
p.phase =  findImgSub(p.phase, p.t1name);
p.magn = findImgSub(p.magn, p.t1name);
fmriname = p.fmriname;
slice_order = p.slice_order;
phase = p.phase;
magn = p.magn;
%X spm('Defaults','fMRI');
%X spm_jobman('initcfg'); % useful in SPM8 only
clear matlabbatch
fmriCell = getsesvolsSubHier(prefix, fmriname);
nSessions = numel(fmriCell);
nVol = sum(cellfun('prodofsize',fmriCell));
fprintf('fMRI has %d sessions for a total of %d volumes\n',nSessions, nVol);
if (nVol < 12)
    error('Too few volumes: this script expects 4D fMRI images');
end
if (nSessions > 5)
    error('Too many sessions: provide the FIRST volume from each 4D image');
end
%determine TR for fMRI data (in secounds)
TRsec = p.TRsec;
TRfmri = getTRSub(deblank (fmriname(1,:)), TRsec);
if  (TRsec == 0)
    if (TRfmri ==0)
        answer = inputdlg('TR (sec)', 'Input required',1,{'2'});
        TRfmri = str2double(answer{1});
    else
        fprintf('Auto-detected repetition time (TR) as %gsec\n', TRfmri); 
    end
    TRsec = TRfmri;
else
    if (TRfmri > 0) && (abs(TRsec - TRfmri) > 0.001)
       warning('Please check TR, header reports %g (not %g)\n', TRfmri, TRsec);
    end
end
%determine slice order
if ~exist('slice_order','var')  || isempty(slice_order) || (slice_order == 0)
    slice_order = getSliceOrderSub(deblank (fmriname(1,:)));
end
if (slice_order ==0)
    answer = inputdlg('slice order (1=[1234],2=[4321],3=[1324],4=[4231], 5=[2413],6=[3142])', 'Input required',1,{'1'});
    slice_order = str2double(answer{1});
end
%end validateInputsSub()

function isNormalized = isNormalizedExists(t1name)
%has the T1 scan been already normalized?
% returns true if brain-extracted T1 ('b' prefix) and T1 deformation exists
[pth,nam,ext] = spm_fileparts(t1name);
defname = fullfile(pth,['y_' nam ext]); %y_ deformation file
if ~exist(defname,'file')
    defname = fullfile(pth,['y_e' nam ext]); %y_e deformation for enantiomorphic normalization
end
bname = fullfile(pth,['b' nam ext]); %'b' = brain extracted image
if (exist(defname,'file')) && (exist(bname,'file'))
    isNormalized = true;
    fprintf('Will use prior normalization estimate for %s\n', t1name)
else
    isNormalized = false;
end
%end()

function slice_order =  getSliceOrderSub(fmriname)
%detect whether slices were acquired ascending, descending, interleaved
[pth,nam,ext] = spm_fileparts( deblank(fmriname(1,:)));
fMRIname1 = fullfile(pth,[ nam, ext]); %'img.nii,1' -> 'img.nii'
fid = fopen(fMRIname1);
fseek(fid,122,'bof');
slice_order = fread(fid,1,'uint8');
fclose(fid);
fMRIname1 = fullfile(pth,[ nam, ext, ',1']); %'img.nii' -> 'img.nii,1'
hdr = spm_vol(fMRIname1);
if (hdr.dim(1) == 90) && (hdr.dim(2) == 90) && (hdr.dim(3) == 50)
    slice_order = -1; %skip: multiband!
    return;
end
if ~isempty(findstr(hdr(1).private.aux_file,'_MB'))
    slice_order = -1; %skip: multiband!
    fprintf('Multi-band detected: slice time correction skipped\n');
    return;    
end
if (slice_order > 0) && (slice_order <= 7)
    fprintf('Auto-detected slice order as %d\n',slice_order);
else
    if (hdr.dim(1) == 80) && (hdr.dim(2) == 80) && (hdr.dim(3) == 42)
        kNIFTI_SLICE_SEQ_DEC = 2; %4,3,2,1
        slice_order = kNIFTI_SLICE_SEQ_DEC;
        fprintf('Assuming UC Irvine descending slice order\n');
    elseif (hdr.dim(1) == 64) && (hdr.dim(2) == 64) && (hdr.dim(3) == 34)
        kNIFTI_SLICE_SEQ_DEC = 2; %4,3,2,1
        slice_order = kNIFTI_SLICE_SEQ_DEC;
        fprintf('Assuming descending slice order\n');
    else
        fprintf('Unable to detect slice order. Please manually specify slice order.\n');
    end;    
end;
%end getSliceOrderSub()


function fimg = findImgSub(img, guess)
fimg = [];
for i = 1: size(img,1)
   img1 = deblank(img(i,:));
   fimg = strvcat(fimg, findImg1Sub(img1, guess)); %#ok<DSTRVCT>
end
%end findImgSub()

function img = findImg1Sub(img, guess)
%find full name for image, validate path
if isempty(img), return; end; %file empty
if exist(img,'file') == 2, return; end; %file found
[pth, nam, ext] = spm_fileparts(img);
img = fullfile(pth, [nam, ext]);
if exist(img,'file') == 2, return; end; %file found - SPM volume number confused us
img = fullfile(pwd, [nam, ext]);
if exist(img,'file') == 2, return; end; %file found - current directory
if exist('guess','var') && ~isempty(guess)
    if exist(guess, 'file') == 7 %guess is a directory
        pth = guess;
    else
        pth = spm_fileparts(guess); %assume guess is a file
    end
    img = fullfile(pth, [nam, ext]);
    if exist(img,'file') == 2, return; end; %file found - guess folder
end
img =  spm_select(1,'image',['Please find image ' nam]);
%end findImg1Sub()

function [tr] =  getTRSub(fmriname, trInput)
%returns Repeat Time in seconds for volume fMRIname
% n.b. for original images from dcm2nii - SPM will strip this information
hdr = spm_vol(fmriname);
hdr = hdr(1);
if (sum(ismember(fieldnames(hdr(1,1).private),'timing')) > 0)  && isfield(hdr(1,1).private.timing,'tspace')
  tr = hdr(1,1).private.timing.tspace;
else %unable to auto-detect TR
  tr = 0;
  if (hdr.dim(1) == 90) && (hdr.dim(2) == 90) && (hdr.dim(3) == 50)
    tr = 1.65;
  end;
  if (hdr.dim(1) == 64) && (hdr.dim(2) == 64) && (hdr.dim(3) == 34)
    tr = 1.85;
  end;
  if tr == 0
    if exist('trInput','var')
        fprintf('%s error: unable to determine TR for image %s (assuming %gs, perhaps SPM stripped this information)\n',mfilename, fmriname, trInput);
    else
        fprintf('%s error: unable to determine TR for image %s (perhaps SPM stripped this information)\n',mfilename,fmriname);
    end
  else
    fprintf('%s warning: guessing TR (%g) based on image dimensions %s\n',mfilename, tr, fmriname);
  end
  fprintf(' Please uncomment "try, N.timing = V.private.timing; end" from spm_create_vol.m\n');
end
%end getTRSub()

function meanname = mocoSub(prefix, fmriname)
%motion correct fMRIdata
fMRIses = getsesvolsSubHier(prefix, fmriname);
fprintf('Motion correction\n');
matlabbatch{1}.spm.spatial.realign.estwrite.data = fMRIses;
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.quality = 0.9;
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.sep = 4;
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.fwhm = 5;
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.rtm = 1;
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.interp = 2;
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.wrap = [0 0 0];
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.weight = '';
matlabbatch{1}.spm.spatial.realign.estwrite.roptions.which = [0 1]; %0 <- do not reslice data!
matlabbatch{1}.spm.spatial.realign.estwrite.roptions.interp = 4;
matlabbatch{1}.spm.spatial.realign.estwrite.roptions.wrap = [0 0 0];
matlabbatch{1}.spm.spatial.realign.estwrite.roptions.mask = 1;
matlabbatch{1}.spm.spatial.realign.estwrite.roptions.prefix = 'r';
spm_jobman('run',matlabbatch);
[pth,nam,ext] = fileparts(deblank(fmriname(1,:)));
meanname = fullfile(pth,['mean', nam, ext]); %moco creates mean image, used for subsequent processing
%end mocoSub()

function prefix = slicetimeSub(prefix, fmriname, TRsec, slice_order)
if slice_order < 1, return; end; %skip time correction
%slice time correct data
kNIFTI_SLICE_SEQ_INC = 1; %1,2,3,4
kNIFTI_SLICE_SEQ_DEC = 2; %4,3,2,1
%kNIFTI_SLICE_ALT_INC = 3; %1,3,2,4 Siemens: interleaved with odd number of slices, interleaved for other vendors
%kNIFTI_SLICE_ALT_DEC = 4; %4,2,3,1 descending interleaved
kNIFTI_SLICE_ALT_INC2 = 5; %2,4,1,3 Siemens interleaved with even number of slices
kNIFTI_SLICE_ALT_DEC2 = 6; %3,1,4,2 Siemens interleaved descending with even number of
[pth,nam,ext] = spm_fileparts( deblank(fmriname(1,:)));
fMRIname1 = fullfile(pth,[ nam, ext]); %'img.nii,1' -> 'img.nii'
hdr = spm_vol([fMRIname1 ',1']);
nslices = hdr.dim(3);
if nslices <= 1 %automatically detect TR
    error('Fatal Error: image %s does not have multiple slices per volume - slice time correction inappropriate. Please edit m-file named %s\n',fMRIname,which(mfilename));
end;
if (slice_order == kNIFTI_SLICE_ALT_INC2) || (slice_order == kNIFTI_SLICE_ALT_DEC2) %sequential
    isSiemens = true;
    if mod(nslices,2) ~= 0, warning('Potential error: Siemens only acquires ALT interleaves with an even number of slice (perhaps multiband you do not want to slice time correct)'); end;
else
    isSiemens = false;
end;
if (slice_order == kNIFTI_SLICE_SEQ_INC) || (slice_order == kNIFTI_SLICE_SEQ_DEC) %sequential
	so = 1:1:nslices;
else % if sequential else Interleaved
	if (mod(nslices,2) == 0) && (isSiemens) %even number of slices, Siemens
		so =[2:2:nslices 1:2:nslices ];
	else
		so =[1:2:nslices 2:2:nslices];
	end
end
if (mod(slice_order,2) == 0) %isDescending
	so = (nslices+1)-so;
end; %isDescending
TA = (TRsec/nslices)*(nslices-1);
fprintf('  Slice order=%d, slices=%d, TR= %0.3fsec, TA= %fsec, referenced to 1st slice.\n', slice_order,nslices, TRsec,TA);
if (TRsec < 0.1) || (TRsec > 5.0)
    fprintf('  Aborting: strange Repeat Time (TR). Please edit the m-file.\n');
    if  (TRsec > 5.0)
          fprintf('  Long TR often used with sparse imaging: if this is a sparse design please set the TA manually.\n');
    end;
    if  (TRsec < 0.1)
          fprintf('  Short TR may be due to DICOM-to-NIfTI conversion. Perhaps use dcm2nii.\n');
    end;
    return;
end; %unusual TR
fMRIses = getsesvolsSubHier(prefix, fmriname);
matlabbatch{1}.spm.temporal.st.scans = fMRIses;
matlabbatch{1}.spm.temporal.st.nslices = nslices;
matlabbatch{1}.spm.temporal.st.tr = TRsec;
matlabbatch{1}.spm.temporal.st.ta = TA;
matlabbatch{1}.spm.temporal.st.so = so;
matlabbatch{1}.spm.temporal.st.refslice = so(1); %set slice order to the first slice http://www.alivelearn.net/?p=1037
%so(1) is first acquired slice, set it as the refernece
%  to see how this works add the line "fprintf('spm_slice_timing slice %d shift %f\n',k, shift amount);' if the 'for k =' loop of spm_slice_timing.m
fprintf('Setting reference slice as %d\n',so(1));
matlabbatch{1}.spm.temporal.st.prefix = 'a';
spm_jobman('run',matlabbatch);
prefix = ['a' prefix]; %new files will have 'a' appended to filename
%end slicetimeSub()

function [fMRIses] = getsesvolsSubHier(prefix, fmriname)
%load all images from all sessions... AS MULTIPLE SESSIONS (e.g. moco)
nsessions = length(fmriname(:,1));
fMRIses = cell(nsessions,1);
for s = 1 : nsessions
	[pth,nam,ext] = spm_fileparts( deblank (fmriname(s,:)));
	sesname = fullfile(pth,[prefix, nam, ext]);
    fMRIses(s,1) = {getsesvolsSub(sesname)};
end;
%end getsesvolsSubHier()

function [fMRIses] = getsesvolsSubFlat(prefix, fmriname)
%load all images from all sessions... AS SINGLE SESSION (e.g. norm writing)
nsessions = length(fmriname(:,1));
fMRIses = '';
for s = 1 : nsessions
	[pth,nam,ext] = spm_fileparts( deblank (fmriname(s,:)));
	sesname = fullfile(pth,[prefix, nam, ext]);
    fMRIses = [fMRIses; getsesvolsSub(sesname)]; %#ok<AGROW>
end;
% end getsesvolsSubFlat()

function [sesvols] = getsesvolsSub(sesvol1)
% input: single volume from 4D volume, output: volume list
%  example 4D file with 3 volumes input= 'img.nii', output= {'img.nii,1';'img.nii,2';'img.nii,3'}
[pth,nam,ext] = spm_fileparts( deblank (sesvol1));
sesname = fullfile(pth,[nam, ext]);
hdr = spm_vol(sesname);
nvol = length(hdr);
sesvols = cellstr([sesname,',1']);
if nvol < 2, return; end;
for vol = 2 : nvol
    sesvols = [sesvols; [sesname,',',int2str(vol)]  ]; %#ok<AGROW>
end;
%end getsesvolsSub()

function namBet = betSub(nam)
%apply brain extraction tool on an image
normIntensitySub(nam);
hdr = spm_vol(nam);
if numel(hdr) > 1, hdr = hdr(1); error('betSub is slow with 4D images - are you sure?'); end;
   %flags.template=fullfile(spm('Dir'),'toolbox','FieldMap','T1.nii');
   flags.template=fullfile(spm('Dir'),'toolbox','OldNorm','EPI.nii');
   flags.fwhm=5;
   flags.nerode=2;
   flags.ndilate=4;
   flags.thresh=0.5;
   flags.reg = 0.05;%0.02; %default can generate "Input to SCHUR must not contain NaN or Inf."
   flags.graphics=0;
pm_brain_mask(hdr, flags);
bnam = prefixSub('bmask', nam);
bimg = spm_read_vols(spm_vol(bnam));
if (min(bimg(:)) == max(bimg(:)) ), error('No variability in %s', bnam); end;
hdr = spm_vol(nam);
img = spm_read_vols(hdr);
thresh = max(bimg(:))/2;
img(bimg < thresh) = min(img(:));
namBet = prefixSub('b', nam);
hdr.fname = namBet;
spm_write_vol(hdr,img);
%end betSub()

function normIntensitySub(nam)
%adjust NIfTI image so brightness ranges from 0..1
% the wild intensity of some Philips scans disrupts pm_brain_mask()
hdr = spm_vol(nam);
img = spm_read_vols(hdr);
if max(img(:)) == min(img(:)), return; end;
img = img-min(img(:)); %translate to 0..max
img = img/max(img(:)); %scale to 0..1
hdr.pinfo = [0; 0; 352];
spm_write_vol(hdr,img);
%end normIntensitySub()

function [longname] = addpth(shortname)
%adds path if not specified
[pth,nam,ext,~] = spm_fileparts(shortname);
if isempty(pth), pth = pwd;end;
longname = fullfile(pth,[nam, ext]);
if exist(longname,'file')~=2; fprintf('Warning: unable to find image %s - cd to approrpiate working directory?\n',longname); end;
longname = [longname, ',1'];
%end addpth()

function stat_1st_levelSub (prefix, basefmriname, kTR, s)
% first level statistics
%  basefmriname : names of 4D fMRI data (one per session)
%  kTR : repeat time in seconds
%  s : structure with statistics
%
%Examples
% stat_1st_blockSub('swa','fMRI.nii', 2)
% stat_1st_blockSub('swa',strvcat('fMRI1.nii','fMRI2.nii'), 2)

if ~isfield(s,'mocoRegress'), s.mocoRegress = false; end;
nCond = numel(s.names);
if nCond ~= size(s.onsets,2)
    error('"s.names" specifies %d conditions while "s.onsets" specifies %d',  nCond, size(s.onsets,2));
end
nSessions = size(s.onsets,1);
fprintf('Experiment has %d sessions with %d conditions\n',nSessions,nCond);
if nSessions ~= size(basefmriname,1)
    error('There must be %d sessions of fMRI data', nSessions);
end
%prepare SPM
if exist('spm','file')~=2; fprintf('%s requires SPM\n',which(mfilename)); return; end;
%spm('Defaults','fMRI');
%spm_jobman('initcfg'); % useful in SPM8 only
clear matlabbatch
%get files if not specified....
if ~exist('kTR','var') || (kTR <= 0)
    error('%s requires the repeat-time (TR) in seconds', mfilename);
end
%next make sure each image has its full path
fmriname = [];
for ses = 1:nSessions
    [pth,nam,ext] = spm_fileparts( deblank (basefmriname(ses,:)));
    if isempty(pth)
        pth = pwd;
    end;
    fmriname = strvcat(fmriname, fullfile(pth, [nam ext]));  %#ok<DSTRVCT>
end
%next - generate output folder
if ~exist('statdirname','var') %no input: select fMRI file[s]
 [pth,nam] = spm_fileparts(deblank(fmriname(1,:)));
 statdirname = nam;
 if s.mocoRegress
    statdirname = [nam '_mc'];
 end
 fprintf('Directory for SPM.mat file not specified - using folder named %s\n',statdirname);
end
predir = pwd;
%create new stat directory
if isempty(pth); pth=pwd; end;
statpth = fullfile(pth, statdirname);
if exist(statpth, 'file') ~= 7; mkdir(statpth); end;
fprintf(' SPM.mat file saved in %s\n',statdirname);
hpf = 128;
if (min([s.duration{:}]) > 2) && (max([s.duration{:}]) < 32);
    if (min([s.duration{:}]) > 5)
        hpf = mean([s.duration{:}]) * 4;
        fprintf('Block design : using %.1fs high pass filter with no temporal derivative.\n',hpf);
    end
	temporalderiv = false;
    if (isfield(s,'forceTemporalDeriv')) && (s.forceTemporalDeriv == true)
        temporalderiv = true;
        warning('Using temporal derivatives');
    end
else
    temporalderiv = true;

	fprintf('Event-related design : using %.1fs high pass filter with a temporal derivative.\n',hpf);
end;
% MODEL SPECIFICATION
%--------------------------------------------------------------------------
clear matlabbatch
matlabbatch{1}.spm.stats.fmri_spec.dir = {statpth};
matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'secs';

matlabbatch{1}.spm.stats.fmri_spec.timing.RT = kTR;
%666 microscale time adjustment
if kTR > 4
    %https://www.jiscmail.ac.uk/cgi-bin/webadmin?A2=spm;27ab2eb2.1010
    matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = round(10*kTR);    
else
    matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = 16;
end
%matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = 16;
matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 1;
warning('fmri_t = %g, fmri_t0 = %g', matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t, matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0);
for ses = 1:nSessions
    %sesFiles = fmriname{ses};%getsesvolsSubSingle(fmriname, ses);
    sesFiles = getsesvolsSubSingle(prefix, fmriname, ses);
    fprintf('Session %d has %d volumes\n',ses, length(sesFiles) );
    %sesFiles = getsesvolsSubSingle(basefmriname, ses);
    matlabbatch{1}.spm.stats.fmri_spec.sess(ses).scans = sesFiles;
    for c = 1:nCond
        matlabbatch{1}.spm.stats.fmri_spec.sess(ses).cond(c).name = deblank(char(s.names{c}));
        matlabbatch{1}.spm.stats.fmri_spec.sess(ses).cond(c).onset = cell2mat(s.onsets(ses, c));
        if numel(s.duration) == 1
            matlabbatch{1}.spm.stats.fmri_spec.sess(ses).cond(c).duration = s.duration{1};
        else
            matlabbatch{1}.spm.stats.fmri_spec.sess(ses).cond(c).duration = cell2mat(s.duration(ses, c));
        end
        matlabbatch{1}.spm.stats.fmri_spec.sess(ses).cond(c).tmod = 0;
        matlabbatch{1}.spm.stats.fmri_spec.sess(ses).cond(c).pmod = struct('name', {}, 'param', {}, 'poly', {});
        matlabbatch{1}.spm.stats.fmri_spec.sess(ses).cond(c).orth = 1; %SPM12
    end;
    matlabbatch{1}.spm.stats.fmri_spec.sess(ses).multi = {''};
    matlabbatch{1}.spm.stats.fmri_spec.sess(ses).regress = struct('name', {}, 'val', {});
    if s.mocoRegress
        [p,n] = spm_fileparts(deblank(fmriname(ses,:)));
        motionFile = fullfile(p, ['rp_', n, '.txt']);
        if ~exist(motionFile, 'file')
            error('Unable to find realign parameters for motion correction %s', motionFile);
        end
        matlabbatch{1}.spm.stats.fmri_spec.sess(ses).multi_reg = {motionFile};
    else
        matlabbatch{1}.spm.stats.fmri_spec.sess(ses).multi_reg = {''};
    end
    matlabbatch{1}.spm.stats.fmri_spec.sess(ses).hpf = hpf;
end
matlabbatch{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
if temporalderiv
	matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [1 0];
else
	matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
end;
matlabbatch{1}.spm.stats.fmri_spec.volt = 1;
matlabbatch{1}.spm.stats.fmri_spec.global = 'None';
matlabbatch{1}.spm.stats.fmri_spec.mthresh = 0.8; %SPM12
matlabbatch{1}.spm.stats.fmri_spec.mask = {''};
matlabbatch{1}.spm.stats.fmri_spec.cvi = 'AR(1)';
% MODEL ESTIMATION
%--------------------------------------------------------------------------
matlabbatch{2}.spm.stats.fmri_est.spmmat = cellstr(fullfile(statpth,'SPM.mat'));
% INFERENCE
%--------------------------------------------------------------------------
matlabbatch{3}.spm.stats.con.spmmat = cellstr(fullfile(statpth,'SPM.mat'));
matlabbatch{3}.spm.stats.con.consess{1}.tcon.name    = 'Task>Rest';
matlabbatch{3}.spm.stats.con.consess{1}.tcon.convec = ones(1,nCond);
nContrast = 1;
if nCond > 1
    for pos = 1: nCond
        for neg = 1: nCond
            if pos ~= neg
                nContrast = nContrast + 1;
                c = zeros(1,nCond);
                c(pos) = 1;
                c(neg) = -1;
                matlabbatch{3}.spm.stats.con.consess{nContrast}.tcon.convec = c;
                matlabbatch{3}.spm.stats.con.consess{nContrast}.tcon.name = [char(s.names{pos}) '>' char(s.names{neg})];
            end % j ~= i
        end %for j
    end %for i
end % > 1 conditions
fprintf('Contrasts:\n');
for c = 1: nContrast
    fprintf(' %d\t%s\n', c, matlabbatch{3}.spm.stats.con.consess{c}.tcon.name); 
end
if temporalderiv %zero pad temporal derivations
    for c = 1 : nContrast
        for cond = 1 : nCond
            v = matlabbatch{3}.spm.stats.con.consess{c}.tcon.convec;
            c2 = cond * 2;
            v = [v(1:(c2-1)) 0 v(c2:end)];
            matlabbatch{3}.spm.stats.con.consess{c}.tcon.convec = v;

        end
    end
end
if s.mocoRegress %add 6 nuisance regressors for motion paramets (rotation + translation]
    for c = 1 : nContrast
    	 matlabbatch{3}.spm.stats.con.consess{c}.tcon.convec = [ matlabbatch{3}.spm.stats.con.consess{c}.tcon.convec  0 0 0 0 0 0];
    end
end
if (nSessions > 1) %replicate contrasts for each session
    for c = 1 : nContrast
        v = matlabbatch{3}.spm.stats.con.consess{c}.tcon.convec;
        for s = 2: nSessions
         matlabbatch{3}.spm.stats.con.consess{c}.tcon.convec = [matlabbatch{3}.spm.stats.con.consess{c}.tcon.convec v];
        end
    end
end
spm_jobman('run',matlabbatch);
cd(predir); %return to starting directory...
%end stat_1st_levelSub()

function [sesvols] = getsesvolsSubSingle(prefix, fmriname, session)
%* Load all fMRI images from single sessions
[pth,nam,ext] = spm_fileparts( deblank (fmriname(session,:))); %#ok<*NASGU>
sesname = fullfile(pth,[prefix, nam, ext]);
hdr = spm_vol(sesname);
nvol = length(hdr);
if (nvol < 2), fprintf('Error 4D fMRI data required %s\n', sesname); return; end;
sesvols = cellstr([sesname,',1']);
for vol = 2 : nvol
    sesvols = [sesvols; [sesname,',',int2str(vol)]  ]; %#ok<AGROW>
end;
%end getsesvolsSubSingle()

function coivox = setOriginSub(vols, modality)
%Align images so that origin and alignment roughly match MNI space
%  vols : cell string of image name(s) - first image used for estimate, others yoked
%  modality : modality of first image 1=T1, 2=T2, 3=EPI
%Example
%  setOrigin('T1.nii',1); %align T1 scan
%  setOrigin({'T1s005.nii', 'fmriblocks009.nii'},1); %use T1 to align T1 and fMRI data
%  setOrigin %use graphical interface
%Chris Rorden 12/2014 (now supports SPM12)
if ~exist('vols','var') %no files specified
 vols = spm_select(inf,'image','Reset origin for selected image(s) (estimated from 1st)');
end
if ischar(vols)
    vols = cellstr(vols);
end
if ~exist('modality','var') %no files specified
 modality = 1;
 fprintf('%s Modality not specified, assuming T1\n', mfilename);
end
coivox = ones(4,1);
%extract filename
[pth,nam,ext, ~] = spm_fileparts(deblank(vols{1}));
fname = fullfile(pth,[nam ext]); %strip volume label
%report if filename does not exist...
if (exist(fname, 'file') ~= 2)
 	fprintf('%s error: unable to find image %s.\n',mfilename,fname);
	return;
end;
hdr = spm_vol([fname,',1']); %load header
img = spm_read_vols(hdr); %load image data
img = img - min(img(:));
img(isnan(img)) = 0;
%find center of mass in each dimension (total mass divided by weighted location of mass
% img = [1 2 1; 3 4 3];
sumTotal = sum(img(:));
coivox(1) = sum(sum(sum(img,3),2)'.*(1:size(img,1)))/sumTotal; %dimension 1
coivox(2) = sum(sum(sum(img,3),1).*(1:size(img,2)))/sumTotal; %dimension 2
coivox(3) = sum(squeeze(sum(sum(img,2),1))'.*(1:size(img,3)))/sumTotal; %dimension 3
XYZ_mm = hdr.mat * coivox; %convert from voxels to millimeters
fprintf('%s center of brightness differs from current origin by %.0fx%.0fx%.0fmm in X Y Z dimensions\n',fname,XYZ_mm(1),XYZ_mm(2),XYZ_mm(3));
for v = 1:   numel(vols)
    fname = deblank(vols{v});
    if ~isempty(fname)
        [pth,nam,ext, ~] = spm_fileparts(fname);
        fname = fullfile(pth,[nam ext]);
        hdr = spm_vol([fname ',1']); %load header of first volume
        fname = fullfile(pth,[nam '.mat']);
        if exist(fname,'file')
            destname = fullfile(pth,[nam '_old.mat']);
            copyfile(fname,destname);
            fprintf('%s is renaming %s to %s\n',mfilename,fname,destname);
        end
        hdr.mat(1,4) =  hdr.mat(1,4) - XYZ_mm(1);
        hdr.mat(2,4) =  hdr.mat(2,4) - XYZ_mm(2);
        hdr.mat(3,4) =  hdr.mat(3,4) - XYZ_mm(3);
        spm_create_vol(hdr);
        if exist(fname,'file')
            delete(fname);
        end
    end
end%for each volume
coregEstTemplateSub(vols, modality);
for v = 1:   numel(vols)
    [pth, nam, ~, ~] = spm_fileparts(deblank(vols{v}));
    fname = fullfile(pth,[nam '.mat']);
    if exist(fname,'file')
        delete(fname);
    end
end %for each volume
%end setOriginSub()

function coregEstTemplateSub(vols, modality)
%vols: images to coregister - first used for estimate
if modality == 2
   template = fullfile(spm('Dir'),'canonical','avg152T2.nii');
elseif modality == 3
    template  = fullfile(spm('Dir'),'toolbox','OldNorm','EPI.nii');
else
    template = fullfile(spm('Dir'),'canonical','avg152T1.nii');
end
if ~exist(template,'file')
    error('Unable to find template named %s\n', template);
end
if ischar(vols)
    vols = cellstr(vols);
end
matlabbatch{1}.spm.spatial.coreg.estimate.ref = {template};
matlabbatch{1}.spm.spatial.coreg.estimate.source = {[deblank(vols{1}),',1']};%{'/Users/rorden/Desktop/3D.nii,1'};
if  numel(vols) > 1
    matlabbatch{1}.spm.spatial.coreg.estimate.other = vols(2:end);% {''};
else
    matlabbatch{1}.spm.spatial.coreg.estimate.other = {''};
end
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];
spm_jobman('run',matlabbatch);
%end coregEstTemplateSub()


function namBet = betSubFSL(nam)
%apply brain extraction tool on an image
[p,n,x] = spm_fileparts(nam);
namBet = fullfile(p, ['bmask', n, '.nii.gz']);
command=sprintf('bet "%s" "%s" -R -f 0.3',nam, namBet);
doFslCmd(command);
namBet = unGzCSub(namBet);
%end betSubFSL()

function [status,cmdout]  = doFslCmd (command, verbose)
if ~exist('verbose', 'var'),
    verbose = true;
end;
fslEnvSub;
cmd = fslCmdSub(command);
if verbose
    fprintf('Running \n %s\n', cmd);
    [status,cmdout]  = system(cmd,'-echo');
else
  [status,cmdout]  = system(cmd);
end
%end doFslCmd()

function fsldir = fslDirSub
if (getenv("FSLDIR"))
    fsldir = strcat(getenv("FSLDIR"),"/");
    if exist(fsldir,"dir")
        return;
    end
end
fsldir= '/usr/local/fsl/'; %CentOS intall location
if ~exist(fsldir,'dir')
    fsldir = '/usr/share/fsl/5.0/'; %Debian install location (from neuro debian)
    if ~exist(fsldir,'dir')
        error('FSL is not in the standard CentOS or Debian locations, please check your installation!');
    end
end
%end fslDirSub()

function fslEnvSub
fsldir = fslDirSub;
curpath = getenv('PATH');
k = strfind(curpath, fsldir);
if isempty(k)
	setenv('PATH',sprintf('%s:%s',fullfile(fsldir,'bin'),curpath));
end
%end fslEnvSub()

function command  = fslCmdSub (command)
fsldir = fslDirSub;
cmd=sprintf('sh -c ". %setc/fslconf/fsl.sh; ',fsldir);
command = [cmd command '"'];
%fslCmdSub

function [fnm, isGz] = unGzCSub (fnm)
if isempty(fnm) || ~exist(fnm, 'file'), return; end;
fnm = deblank(fnm);
[pth,nam,ext] = spm_fileparts(fnm);
if ~strcmpi(ext,'.gz'), return; end;
ofnm = fnm;
fnm = char(gunzip(fnm));
delete(ofnm);
%end unGzSub()