function matName = nii_preprocess(imgs, matName, checkForUpdates, hideInteractiveGraphs)
%preprocess data from multiple modalities3
% imgs.T1: filename of T1 scan - ONLY REQUIRED INPUT: ALL OTHERS OPTIONAL
% imgs.T2: filename used to draw lesion, if absent lesion drawn on T1
% imgs.Lesion : lesion map
% imgs.ASL : pCASL or PASL sequence
% imgs.ASLrev : phase reversed sequence
% imgs.Rest : Resting state sequence
% DTI : Diffusion scan
% DTIrev : Diffusion scan with reverse phase encoding of 'DTI'
% fMRI : fMRI scan
%Examples
% imgs.T1 = 'T1.nii'; imgs.ASL = 'ASL.nii';
% nii_preprocess(imgs);

%check dependencies
fprintf('%s version 27May2021\n', mfilename);   
%warning('Not checking for updates 6666');
if ~exist('checkForUpdates','var') || checkForUpdates
    checkForUpdate(fileparts(mfilename('fullpath')));
end
nii_check_dependencies;
if nargin < 1 %, error('Please use nii_preprocess_gui to select images'); 
    %imgs = 'T1_limegui.mat';
    [f, p] = uigetfile('*limegui.mat', 'Select a mat file');
    imgs = fullfile(p,f);
end;
if ~isstruct(imgs) %user passes folder, e.g. "nii_preprocess(pwd)"
    if iscell(imgs), imgs = imgs{1}; end;
    pth = imgs;
    if ~exist(pth, 'dir'), error('Unable to find folder %s', pth); end;
    imgs = findImgs(pth);
    if isempty(imgs) 
        error('Unable to find NIfTI images (e.g. "T1_.nii") in folder %s\n', pth);
    end
end 
if isempty(spm_figure('FindWin','Graphics')), 
    spm fmri; 
    %spm_get_defaults('cmdline',true); %enable command line mode in scripts
end; %launch SPM if it is not running
%f = spm_figure('FindWin','Graphics'); clf(f.Number); %clear SPM window
if ~exist('hideInteractiveGraphs','var') || hideInteractiveGraphs
    fg = spm_figure('FindWin','Interactive');
    if ~isempty(fg), close(fg); end
end
if ischar(imgs), imgs = load(imgs); end;
%set structures
if ~isfield(imgs,'T1') || isempty(imgs.T1), error('T1 scan is required'); end;
if ~isfield(imgs,'T2'), imgs.T2 = []; end;
if ~isfield(imgs,'Lesion'), imgs.Lesion = []; end;
if ~isfield(imgs,'ASL'), imgs.ASL = []; end;
if ~isfield(imgs,'ASLrev'), imgs.ASLrev = []; end;
if ~isfield(imgs,'Rest'), imgs.Rest = []; end;
if ~isfield(imgs,'DTI'), imgs.DTI = []; end;
if ~isfield(imgs,'DTIrev'), imgs.DTIrev = []; end;
if ~isfield(imgs,'fMRI'), imgs.fMRI = []; end;

if ~exist('matName','var') || isempty(matName)
    [p,n] = filepartsSub(imgs.T1);
    if ~exist(p,'dir')
        error('Please provide a matName: files no longer located in %s', p);
    end
    if endsWith(n,'_'), n = n(1:end-1); end
    matName = fullfile(p, [n, '_lime.mat']);
end
if true
    diary([matName, '.log.txt'])
    imgs = unGzAllSub(imgs); %all except DTI - fsl is OK with nii.gz
    tStart = tic;
    imgs = doT1Sub(imgs, matName); %normalize T1
    imgs = doI3MSub(imgs, matName);
    tStart = timeSub(tStart,'T1');
    imgs = doAslSub(imgs, matName);
    tStart = timeSub(tStart,'ASL');   
    %666x -> recommented in by Roger for Jill's Data
    imgs = doRestSub(imgs, matName); %TR= 1.850 sec, descending; %doRestSub(imgs, matName, 2.05, 5); %Souvik study
    %666x <-
    
    %tStart = timeSub(tStart,'REST');

    %666x -> recommented in by Roger for Data with fMRI
    imgs = dofMRISub(imgs, matName);
    %666x <-
    
    %666 VBM
    imgs = doVBMSub(imgs, matName);   
    
    tStart = timeSub(tStart,'fMRI');
    %warning('Skipping DTI');
    if true
        if ~isempty(imgs.DTI)
            if ~nii_check_dims({imgs.DTI; imgs.DTIrev}), error('Fix DTI'); end;
            imgs = removeDotDtiSub(imgs);
            dtiDir = fileparts(imgs.DTI);
            doDtiSub(imgs);
            doFaMdSub(imgs, matName);
            %-->(un)comment next line for JHU tractography
            doDtiTractSub(imgs, matName, dtiDir, 'jhu');
            %-->(un)comment next line for AICHA tractography
            doDtiTractSub(imgs, matName, dtiDir, 'AICHA'); % uncommented by GY, Aug10_2016
            %-->compute scalar DTI metrics
            
            doTractographySub(imgs);
            %666 doDkiSub(imgs, matName);
            tStart = timeSub(tStart,'DTI');
        end
        %666 doDkiSub(imgs, matName, true);
    end
    tStart = timeSub(tStart,'DKI');
    %matName
    %addLimeVersionSub(matName); %update versioning
end
%print output
pth = '~/Desktop';
if ~exist(pth,'file'), pth = ''; end;
%printDTISub(imgs, fullfile(pth,'MasterDTI')); %show results - DTI
nii_mat2ortho(matName, fullfile(pth,'MasterNormalized')); %do after printDTI (spm_clf) show results - except DTI
diary off
%nii_preprocess()

function imgs = findImgs(pth);
%find images in path
imgs = [];
ms = {'T1', 'T2', 'Lesion', 'ASL', 'ASLrev', 'Rest', 'DTI','DTIrev', 'fMRI', 'DKI','DKIrev'};
for m = 1 : numel(ms)
    fnm = fullfile(pth,[ms{m}, '_*.n*']);
    fnms = nii_dir(fnm, false);
    if isempty(fnms), continue; end;
    if numel(fnms) > 1 
       error('Multiple NIfTI files have the name %s\n', fnm); 
    end
    imgs.(ms{m}) = fnms{1};
end

%end findImgs
function addLimeVersionSub(matName)
%add 'timestamp' to file allowing user to autodetect if there mat files are current
% e.g. after running
%  m = load(matName);
%  fprintf('LIME version (YYYY.MMDD): %.4f\n', m.T1.lime);
v = nii_matver;
stat.T1.lime = v.lime;
fprintf('adding lime version to %s\n', matName);
if exist(matName,'file')
    old = load(matName);
    stat = nii_mergestruct(stat,old); %#ok<NASGU>
end
save(matName,'-struct', 'stat');
%end addLimeVersionSub()

function checkForUpdate(repoPath)
prevPath = pwd;
cd(repoPath);
if exist('.git','dir') %only check for updates if program was installed with "git clone"
    [s, r] = system('git fetch origin','-echo');
    if strfind(r,'fatal')
        warning('Unabe to check for updates. Network issue?');
        return;
    end
    [~, r] = system('git status','-echo');
    if strfind(r,'behind')
        if askToUpdate
            [~, r] = system('git pull','-echo');
            showRestartMsg
        end
    end
else %do nothing for now
    warning(sprintf('To enable updates run "!git clone git@github.com:neurolabusc/%s.git"',mfilename));
end
cd(prevPath);
%end checkForUpdate()

function a = askToUpdate
% Construct a questdlg
choice = questdlg(sprintf('An update for %s is available. Would you like to update?',mfilename), ...
    'Auto update', ...
    'Yes','No','Yes');
% Handle response
switch choice
    case 'Yes'
        a = true;
    case 'No'
        a = false;
end
%end askToUpdate()

function showRestartMsg
uiwait(msgbox('The program must be restarted for changes to take effect. Click "OK" to quit the program. You will need to restart it just as you normally would','Restart Program'))
exit;
%end showRestartMsg()

function tStart = timeSub(tStart, timeComment)
elapsed =  toc(tStart);
if elapsed > 1.0
    fprintf('Stage %s required\t%g\tseconds\n', timeComment, elapsed);
end
tStart = tic;
%timeSub

function printDTISub(imgs, matName)
if isempty(spm_figure('FindWin','Graphics')), spm fmri; end; %launch SPM if it is not running
if  isempty(imgs.DTI) , return; end; %required
spm_clf;
spm_figure('Clear', 'Graphics');
spm_orthviews('Reset');
f = spm_figure('FindWin','Graphics'); clf(f.Number); %clear SPM window
FA = prepostfixSub('n', 'd_FA', imgs.DTI);
ROI = prepostfixSub('', '_roi', imgs.DTI);
if ~exist(FA,'file') || ~exist(ROI,'file') , return; end; %required
FA = unGzSub (FA);
ROI = unGzSub (ROI);
spm_clf;
spm_figure('Clear', 'Graphics');
spm_orthviews('Reset');
fig = spm_figure('FindWin','Graphics');
ax  = axes('Position',[0.1 0.5 0.35 0.3],'Visible','off','Parent',fig);
orthSub(FA,[0.01 0.33 .48 .49]); %spm_orthviews('Image',FA,[0.01 0.33 .48 .49]);
orthSub(ROI,[.51 0.33 .48 .49]); %spm_orthviews('Image',ROI,[.51 0.33 .48 .49]);
%[~, nam] = fileparts(FA);
%text(0,0.95, sprintf('DTI %s',nam),'Parent',ax, 'FontSize',10, 'fontn','Arial');
spm_print(matName);
%end printDTISub()

function orthSub(fnm, pos)
% e.g. mat2niiSub(mat, 'i3mT1', [0.51 0.01 .4 .49]
if ~exist(fnm,'file'), return; end;
spm_orthviews('Image',fnm,pos);
fig = spm_figure('FindWin','Graphics');
ax  = axes('Position',[0.0 0.1 1.0 1.0],'Visible','off','Parent',fig);
[~,n] = filepartsSub(fnm);
text(pos(1)+0.0,pos(2)+0.01, n,'Parent',ax, 'FontSize',8, 'fontn','Arial', 'color','red');
%end orthSub()


%adding RNN 666 
function imgs = doVBMSub(imgs, matName)

%these lines confirm presence of eT1 (stolen from doDTISub)
eT1 = prefixSub('e',imgs.T1); %enantimorphic image
if ~exist(eT1,'file'), eT1 = imgs.T1; end; %if no lesion, use raw T1
if ~exist(eT1,'file'), fprintf('doVBM unable to find %s\n', eT1); return; end; %required

global ForceVBM; %e.g. user can call "global ForceVBM;  ForceVBM = true;"

%call nii_VBM, only needs eT1 as parameter
%nii_VBM(eT1);

%add results to matName

%end doVBMSub()



function imgs = dofMRISub(imgs, matName)
if isempty(imgs.T1) || isempty(imgs.fMRI), return; end; %required
imgs.fMRI = removeDotSub (imgs.fMRI);
[~,n] = fileparts(imgs.fMRI);
[p] = fileparts(imgs.fMRI);
cstat = fullfile(p,[n], 'con_0002.nii');
bstat = fullfile(p,[n], 'beta_0001.nii');
global ForcefMRI; %e.g. user can call "global ForcefMRI;  ForcefMRI = true;"
%do not check cstat files anymore, as folder names may have changed...
%if isempty(ForcefMRI) && isFieldSub(matName, 'fmri') && exist(cstat, 'file') && exist(bstat,'file'), fprintf('Skipping fMRI (already done) %s\n',imgs.fMRI); return;  end;
if isempty(ForcefMRI) && isFieldSub(matName, 'fmri'), fprintf('Skipping fMRI (already done) %s\n',imgs.fMRI); return;  end;
if ~exist(imgs.fMRI,'file'), warning('Unable to find %s', imgs.fMRI); return; end;
%if ~isempty(ForcefMRI)
    d = fullfile(p,n);
    if exist(d,'file'), rmdir(d,'s'); end; %delete statistics directory
    delImgs('sw', imgs.fMRI);
    delMat(imgs.fMRI)
    %delImgs('swa', imgs.fMRI); %we never slice-time correct sparse data!
%end
if ~exist('nii_fmri60.m','file')
    fnm = fullfile(fileparts(which(mfilename)), 'nii_fmri');
    if ~exist(fnm,'file')
        error('Unable to find %s', fnm);
    end
    addpath(fnm);
end

%GY, Oct 5, 2017
%nii_fmri60(imgs.fMRI, imgs.T1, imgs.T2); %use fMRI for normalization
nii_fmri60(imgs.fMRI, imgs.T1, imgs.T2, imgs.Lesion);

if ~exist(cstat, 'file') || ~exist(bstat,'file')
    error('fMRI analysis failed : %s\n  %s', bstat, cstat);
end
nii_nii2mat(cstat, 'fmri' , matName); %12
nii_nii2mat(bstat, 'fmrib', matName); %13
vox2mat(prefixSub('wmean',imgs.fMRI), 'fMRIave', matName);
vox2mat(prefixSub('wbmean',imgs.fMRI), 'fMRIave', matName);
%end dofMRISub()

function XYZmm = getCenterOfIntensitySub(vols)
XYZmm = ones(3,1);
if ischar(vols), vols = cellstr(vols); end;
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
coivox = ones(4,1);
coivox(1) = sum(sum(sum(img,3),2)'.*(1:size(img,1)))/sumTotal; %dimension 1
coivox(2) = sum(sum(sum(img,3),1).*(1:size(img,2)))/sumTotal; %dimension 2
coivox(3) = sum(squeeze(sum(sum(img,2),1))'.*(1:size(img,3)))/sumTotal; %dimension 3
XYZmm = hdr.mat * coivox; %convert from voxels to millimeters
XYZmm = XYZmm(1:3);
%end setCenterOfIntensitySub()

function doTractographySub(imgs)
if isempty(imgs.DTI), return; end; %required
FAimg = prepostfixSub('', 'd_FA', imgs.DTI);
if ~exist(FAimg,'file') , return; end; %required
[p,n] = fsl_filepartsSub(imgs.DTI);
basename = fullfile(p,n);
vtkname = [basename, '.vtk'];
%next: check for bedpost
BPimg = fullfile(p, 'bedpost.bedpostX', 'dyads1.nii');
if ~exist(BPimg,'file'), BPimg = fullfile(p, 'bedpost.bedpostX', 'dyads1.nii.gz'); end;
if exist(BPimg,'file'), basename = BPimg; end;
global ForceDTI;
if isempty(ForceDTI) && exist(vtkname,'file') , fprintf('Skipping tractography (already done): %s\n', vtkname); return; end; %required
if ismac
    TrakExe = fullfile(fileparts(which(mfilename)), 'tracktionOSX');
elseif isunix % Code to run on Linux plaform
    TrakExe = fullfile(fileparts(which(mfilename)), 'tracktionLX');
end
PathTrakExe = findExeSub(TrakExe);
if isempty(PathTrakExe), fprintf('Tractography skipped: unable to find %s', TrakExe); return; end;
cmd = sprintf('%s -o "%s" "%s"',TrakExe, vtkname, basename);
[status, fullnam]  = system(cmd,'-echo');
%end doTractographySub()

%old version uses SingleTensorFT that does not preserve NIfTI orientation
%and can not use BEDPOSTX
% function doTractographySub(imgs)
% if isempty(imgs.DTI), return; end; %required
% Mask = prepostfixSub('', '_FA_thr', imgs.DTI);
% if ~exist(Mask,'file') , return; end; %required
% [p,n] = fsl_filepartsSub(imgs.DTI);
% basename = fullfile(p,n);
% vtkname = [basename, '.vtk'];
% if exist(vtkname,'file') , fprintf('Skipping tractography (alredy done): %s\n', vtkname); return; end; %required
% imgCnv = [basename, '_dtitk.nii.gz'];
% ConvExe = 'TVFromEigenSystem';
% PathConvExe = findExeSub(ConvExe);
% if isempty(PathConvExe), fprintf('Tractography skipped: unable to find %s', ConvExe); return; end;
% TrakExe = 'SingleTensorFT';
% PathTrakExe = findExeSub(TrakExe);
% if isempty(PathTrakExe), fprintf('Tractography skipped: unable to find %s', TrakExe); return; end;
% cmd = sprintf('%s -basename "%s" -out %s -type FSL',ConvExe,basename, imgCnv);
% [status, fullnam]  = system(cmd,'-echo');
% %e.g. TVFromEigenSystem -basename 199_99_AP_7 -type FSL
% cmd = sprintf('%s -in "%s" -seed "%s" -out "%s"',TrakExe, imgCnv, Mask, vtkname );
% %e.g. SingleTensorFT -in 199_99_AP_7.nii.gz -seed 199_99_AP_7_FA_thr.nii.gz -out t.vtk
% [status, fullnam]  = system(cmd,'-echo');
% %end doTractographySub()

function fullnam = findExeSub(nam)
%findExeSub('nano') returns path of executable ("/usr/bin/nano") or empty if not found
cmd=sprintf('which "%s"',nam);
[status, fullnam]  = system(cmd);
if (fullnam(1) ~= filesep), fullnam = []; end;
%end findExeSub


function nam = prepostfixSub (pre, post, nam)
if isempty(nam), return; end;
[p, n, x] = filepartsSub(nam);
nam = fullfile(p, [pre, n, post, x]);
if ~exist(nam, 'file') && strcmpi(x,'.nii')
   nam = fullfile(p, [pre, n, post, '.nii.gz']);
   if ~exist(nam, 'file')
        nam = fullfile(p, [pre, n, post, x]);
   end
end
if ~exist(nam, 'file') && strcmpi(x,'.nii.gz')
   nam = fullfile(p, [pre, n, post, '.nii']);
   if ~exist(nam, 'file')
        nam = fullfile(p, [pre, n, post, x]);
   end
end
%end prepostfixSub()

function fname = rescaleSub(fname)
hdr = spm_vol(fname);
img = spm_read_vols(hdr);
if min(img(:)) < 0, error('Image intensity must be positive'); end;
if min(img(:)) == max(img(:)), error('Image has no variability'); end;
img = img - min(img(:)); %scale from zero
img = img/max(img(:)); %scale from 0..1
[pth, nm, ext] = spm_fileparts(fname);
img = power(img, 0.5);
fname = fullfile(pth, ['n' nm ext]);
hdr.fname = fname;
spm_write_vol(hdr,img);
%end rescaleSub()

function doDkiSub(imgs, matName, isDki)
if exist('isDki','var') && (isDki)
    if ~isfield(imgs,'DKI'), return; end;
    if isempty(imgs.DKI), return; end; %nothing to do!
    global ForceDTI;
    if isempty(ForceDTI) && isFieldSub(matName, 'mk')
        fprintf('Skipping MK estimates: already computed\n');
        return;
    end; %stats already exist
    warning('If you are using fsl 5.0.9, ensure you have patched version of eddy ("--repol" option) https://fsl.fmrib.ox.ac.uk/fsldownloads/patches/eddy-patch-fsl-5.0.9/centos6/');
    %2017: dti_1_eddy_cuda now auto-detects if cuda is installed
    %if isEddyCuda7Sub()
    if HalfSphere(imgs.DKI)
        command= [fileparts(which(mfilename)) filesep 'dti_1_eddy_cuda_half.sh'];
    else
        command= [fileparts(which(mfilename)) filesep 'dti_1_eddy_cuda.sh'];
        %command= [fileparts(which(mfilename)) filesep 'dti_1_eddy_correct.sh'];
    end
    %else
    %    command= [fileparts(which(mfilename)) filesep 'dti_1_eddy.sh'];
    %end
    command=sprintf('%s "%s"',command, imgs.DKI);
    doFslCmd (command);
    doDkiCoreSub(imgs.T1, imgs.DKI, matName)
else
    doDkiCoreSub(imgs.T1, imgs.DTI, matName);
end
%end doDkiSub()

function doDkiCoreSub(T1, DTI, matName)
if isempty(T1) || isempty(DTI), return; end; %required
if isFieldSub(matName, 'mk'), fprintf('skipping DKI: already computed\n'); return; end; %stats already exist
wbT1 = prefixSub('wb',T1); %warped brain extracted image
if ~exist('dkifx2','file'),  error('skipping DKI: requires dkifx2 script\n'); return; end;
mask=prepostfixSub('', 'db_mask', DTI);
mask = unGzSub (mask);
if ~exist(mask,'file'),  error('unable to find %s\n', mask); return; end;
MKfa=prepostfixSub('', 'd_FA', DTI);
if ~exist(MKfa,'file'),  error('unable to find %s\n', MKfa); return; end;
dti_u=prepostfixSub('', 'du', DTI);
[pth,nam] = filepartsSub(DTI);
bvalnm = fullfile(pth, [nam, 'both.bval']); %assume topup
if ~exist(bvalnm, 'file')
    bvalnm = fullfile(pth, [nam, '.bval']);
end
if ~exist(mask, 'file') || ~exist(dti_u, 'file') || ~exist(bvalnm, 'file')|| ~exist(wbT1, 'file')
    fprintf('Skipping DKI preprocessing. Missing file(s): %s %s %s %s\n', mask, dti_u, bvalnm, wbT1);
end;
bval=load(bvalnm);
if max(bval(:)) < 1500
    fprintf('Skipping DKI preprocessing. Require B-values >1500. Max b-value: %d\n',max(bval(:)));
    return;
end
MK=dkifx2(dti_u, bvalnm, mask, true);
MKmask=mask;
%MKmask = prepostfixSub('', 'u_ldfDKI_MASK', DTI);
%MK=prepostfixSub('', 'u_ldfDKI_MK', DTI);
if ~exist(MK, 'file') || ~exist(MKmask, 'file')
    fprintf('Serious error: no kurtosis images named %s %s\n', MK, MKmask);
    return;
end;
if ~exist(wbT1,'file'), error('unable to find %s',wbT1); end;
%normalize mean kurtosis to match normalized, brain extracted T1
wMK = prepostfixSub('w', '', MK);
oldNormSub( {MKfa, MK}, wbT1, 8, 8 );
nii_nii2mat(wMK, 'mk', matName);
%save note
fid = fopen('dki.txt', 'a+');
fprintf(fid, '%s\n', matName);
fclose(fid);
%end doDkiCoreSub()

function doFaMdSub(imgs, matName)
if isempty(imgs.T1) || isempty(imgs.DTI), return; end; %required
T1 = prefixSub('wb',imgs.T1); %warped brain extracted image
FA = prepostfixSub('', 'd_FA', imgs.DTI);
MD = prepostfixSub('', 'd_MD', imgs.DTI);
if ~exist(T1,'file') || ~exist(FA,'file') || ~exist(MD,'file'), return; end; %required
global ForceDTI;
if isempty(ForceDTI) && isFieldSub(matName, 'fa')
    fprintf('Skipping MD/FA estimates: already computed\n');
    return;
end; %skip: previously computed
FA = unGzSub (FA);
nii_famask(FA, true); %8/2016: remove speckles at rim of cortex
MD = unGzSub (MD);
wFA = prepostfixSub('w', '', FA);
wMD = prepostfixSub('w', '', MD);
if ~exist(wFA,'file') || ~exist(wMD,'file')
    nFA = rescaleSub(FA);
    %nii_setOrigin12({nFA, FA, MD}, 1,false); %rFA looks like T1
    oldNormSub( {nFA, FA,MD}, T1, 8, 10 );
end
nii_nii2mat(wFA, 'fa', matName); %6
nii_nii2mat(wMD, 'md', matName); %8
%end doFaMdSub()

function doDtiSub(imgs)
if isempty(imgs.T1) || isempty(imgs.DTI), return; end; %required
betT1 = prefixSub('b',imgs.T1); %brain extracted image
if ~exist(betT1,'file'), fprintf('doDti unable to find %s\n', betT1); return; end; %required
eT1 = prefixSub('e',imgs.T1); %enantimorphic image
if ~exist(eT1,'file'), eT1 = imgs.T1; end; %if no lesion, use raw T1
if ~exist(eT1,'file'), fprintf('doDti unable to find %s\n', eT1); return; end; %required
global ForceDTI;
if isempty(ForceDTI) && isDtiDoneBedpost(imgs), fprintf('Skipping DTI processing (bedpost done)\n'); return; end;
n = bvalCountSub(imgs.DTI);
if (n < 1)
    fprintf('UNABLE TO FIND BVECS/BVALS FOR %s\n', imgs.DTI);
    return
end
if (n < 12)
    fprintf('INSUFFICIENT BVECS/BVALS FOR %s\n', imgs.DTI);
    return
end
%preprocess - denoise
dti_d=prepostfixSub('', 'd', imgs.DTI);
if exist(imgs.DTIrev)
    dti_dr=prepostfixSub('', 'd', imgs.DTI);
end;
if isempty(ForceDTI) && exist(dti_d, 'file')
   fprintf('Skipping DTI denoising: found %s\n', dti_d);
else
    mm = imgMM(imgs.DTI);
    degibbs = (mm > 1.9); %partial fourier used for HCP 1.5mm isotropic images
    dti_d = nii_dwidenoise (imgs.DTI, degibbs);
    if exist(imgs.DTIrev)
        dti_dr = nii_dwidenoise (imgs.DTIrev, degibbs);
    end;
end;
%preprocess - eddy
dti_u=prepostfixSub('', 'du', imgs.DTI);
if isempty(ForceDTI) && exist(dti_u, 'file')
    fprintf('Skipping DTI preprocessing: found %s\n', dti_u);
else
    clipSub (imgs); %topup requires images with even dimensions
    %1 - eddy current correct
    %2017: dti_1_eddy_cuda now auto-detects if cuda is installed
    %if isEddyCuda7Sub()
    %9/2017: always use Eddy
    %if isFullSphereSub(imgs.DTI)
    if HalfSphere(imgs.DTI)
        command= [fileparts(which(mfilename)) filesep 'dti_1_eddy_cuda_half.sh'];
    
    else
        command= [fileparts(which(mfilename)) filesep 'dti_1_eddy_cuda.sh'];
    end;
    %else
    %    command= [fileparts(which(mfilename)) filesep 'dti_1_eddy_correct.sh'];
    %end
    %else
    %    command= [fileparts(which(mfilename)) filesep 'dti_1_eddy.sh'];
    %end
    if isempty(imgs.DTIrev)
        command=sprintf('%s "%s"',command, dti_d);
    else
        nr = bvalCountSub(dti_dr);
        if (nr ~= n)
            fprintf('WARNING: BVECS/BVALS DO NOT MATCH %s %s\n', dti_d, dti_dr);
            %return
        end
        command=sprintf('%s "%s" "%s"',command, dti_d, dti_dr);
    end
    cleanupDtiDir(imgs.DTI);
    doFslCmd (command);
end
doDtiBedpostSub(imgs.DTI);
%end doDtiSub()

function mm = imgMM(fnm)
%mean voxel size, 2x2x2=2, 2x2x4=2.52
h = spm_vol([fnm ',1']);
mm = h.mat(1:3, 1:3) * [1 1 1]';
mm = nthroot(prod(abs(mm)), 3);
%end imgMM

function isEddyCuda = isEddyCuda7Sub
%eddy_cuda7.0 is a beta version supplied by Jesper Andersson that supports the GTX970 and has new features
isEddyCuda = false;
if ~isGpuInstalledSub(), return; end;
eddyName = '/usr/local/fsl/bin/eddy_cuda7.0';
isEddyCuda = exist(eddyName,'file') > 0;
if ~isEddyCuda, printf('Hint: Eddy will run faster if you install %s', eddyName); end;
%end isEddyCuda7Sub

function done = isDtiDoneBedpost(imgs)
done = false;
if isempty(imgs.DTI), return; end;
pth = fileparts( imgs.DTI );
bed_dirX=fullfile(pth, 'bedpost.bedpostX');
bed_done=fullfile(bed_dirX, 'xfms', 'eye.mat');
if exist(bed_done, 'file')
    done = true;
end;
%end isDtiDone()


%{ 
function done = isDtiDone(imgs)
done = false;
if isempty(imgs.DTI), return; end;
p = fileparts( imgs.DTI );
pDir = fullfile(p,'probtrackx');
if exist(pDir, 'file')
    done = true;
end;
%end isDtiDone()
%}

function clipSub (imgs)
if isempty(imgs.DTIrev)
    nii_clipeven(imgs.DTI, true);
else
    nii_clipeven({imgs.DTI, imgs.DTIrev}, true);
end
%end clipSub()

% function doDtiBedpostSub(dti) %warp atlas to DTI
% t_start=tic;
% pth = fileparts(dti);
% bed_dir=fullfile(pth, 'bedpost');
% if ~exist(bed_dir, 'file'), mkdir(bed_dir); end;
% bed_dirX=fullfile(pth, 'bedpost.bedpostX');
% if ~exist(bed_dirX, 'file'), mkdir(bed_dirX); end;
% bed_done=fullfile(bed_dirX, 'xfms', 'eye.mat');
% if exist(bed_done,'file'), fprintf('Skipping bedpost (already done)\n'), return, end;
% dti_u=prepostfixSub('', 'u', dti);
% dti_x=fullfile(bed_dir, 'data.nii.gz');
% if ~exist(dti_u,'file'), error('Bedpost unable to find %s', dti_u); end;
% copyfile(dti_u, dti_x);
% [bvec, bval] = getBVec(dti);
% dti_x=fullfile(bed_dir, 'bvecs');
% copyfile(bvec, dti_x);
% dti_x=fullfile(bed_dir, 'bvals');
% copyfile(bval, dti_x);
% dti_faThr=prepostfixSub('', '_FA_thr', dti);
% dti_x=fullfile(bed_dir, 'nodif_brain_mask.nii.gz');
% copyfile(dti_faThr, dti_x);
% command = fullfile(fslDirSub, 'bin', 'bedpostx1');
% if ~exist(command,'file'), error('You need the special file %s', command); end;
% command=sprintf('bedpostx1 "%s" ', bed_dir);
% doFslCmd (command);
% bed_cmd=fullfile(bed_dirX, 'commands.txt');
% if ~exist(bed_cmd, 'file'), error('Unable to find %s', bed_cmd); end;
% maxThreads = feature('numCores');
% if maxThreads > 15, maxThreads = maxThreads - 1; end;
% nPerm = 5000;
% %cmds = textread(bed_cmd, '%s', 'delimiter', '\n'); %deprecated - use textscan
% cmds = dataread('file', bed_cmd, '%s', 'delimiter', '\n'); %deprecated - use textscan
% %fid = fopen(bed_cmd);
% %cmds = textscan(fid,'%s','delimiter','\n');
% %fclose(fid);
% threads = {};
% for i = 1: numel(cmds)
%     if mod(i,maxThreads) == 0, fprintf('bedpost %d/%d (%d threads)\n', i, numel(cmds), maxThreads); end;
%     doFslCmd ([cmds{i} ' &'], false);
%     isDone =fullfile(bed_dirX, 'diff_slices', ['data_slice_' sprintf('%4.4d',i-1)],'dyads1.nii.gz'); %file created when thread completes
%     threads = [threads; {isDone}]; %#ok<AGROW>
%     threads = threadSub(threads, maxThreads);
% end
% threadSub(threads, 0);
% pause(2);
% command = fullfile(fslDirSub, 'bin', 'bedpostx_postproc.sh');
% if ~exist(command,'file'), error('You need the special file %s', command); end;
% command=[command ' ', bed_dir];
% doFslCmd (command);
% fprintf ('Bedpost took %f seconds to run.\n', toc(t_start) );
%end doDtiBedpostSub()

function doDtiBedpostSub(dti) %warp atlas to DTI
t_start=tic;
pth = fileparts(dti);
bed_dirX=fullfile(pth, 'bedpost.bedpostX');
%if exist(bed_dirX, 'file'), rmdir(bed_dirX, 's'); end; %666 ForceBedpost
bed_done=fullfile(bed_dirX, 'xfms', 'eye.mat');
global ForceDTI;
if isempty(ForceDTI) && exist(bed_done,'file'), fprintf('Skipping bedpost (already done)\n'), return, end;
if ~exist(bed_dirX, 'file'), mkdir(bed_dirX); end;
bed_dir=fullfile(pth, 'bedpost');
%if exist(bed_dir, 'file'), rmdir(bed_dir, 's'); end; %666 ForceBedpost
if ~exist(bed_dir, 'file'), mkdir(bed_dir); end;
dti_u=prepostfixSub('', 'du', dti);
dti_x=fullfile(bed_dir, 'data.nii.gz');
if ~exist(dti_u,'file'), error('Bedpost unable to find %s', dti_u); end;
copyfile(dti_u, dti_x);
[bvec, bval] = getBVec(dti);
dti_x=fullfile(bed_dir, 'bvecs');
copyfile(bvec, dti_x);
dti_x=fullfile(bed_dir, 'bvals');
copyfile(bval, dti_x);
dti_faThr=prepostfixSub('', 'd_FA_thr', dti);
dti_x=fullfile(bed_dir, 'nodif_brain_mask.nii.gz');
copyfile(dti_faThr, dti_x);
if isGpuInstalledSub
    command=sprintf('bedpostx_gpu "%s" ', bed_dir);
else
    command=sprintf('bedpostx "%s" ', bed_dir);
end
fslParallelSub;
doFslCmd (command);
if ~exist(bed_done,'file')
    error('Fatal error running bedpostx');
end
%while ~exist(bed_done,'file')
%    pause(1.0);
%end
%end doDtiBedpostSub()

function [bvec, bval] = getBVec(dti)
%both.bval
dti = unGzNameSub(deblank(dti));
[p,n] = fileparts(dti);
dti = fullfile(p,n);
bvec = [dti 'du.eddy_rotated_bvecs'];
if ~exist(bvec,'file')
    bvec = [dti 'both.bvec'];
end
bval = [dti 'dboth.bval'];
if exist(bvec,'file') && exist(bval,'file'), return; end; %combined AP/PA bvecs
%fall back to original bvecs and bvals
if ~exist(bvec,'file')
    bvec = [dti '.bvec'];
end
bval = [dti '.bval'];
if ~exist(bvec,'file') || ~exist(bval,'file'), error('Can not find files %s %s', bvec, bval); end;
%end getBVec()

function doDtiTractSub(imgs, matName, dtiDir, atlas)
dti = imgs.DTI;
global ForceDTI;
if ~exist('atlas','var'),
    atlas = 'jhu';
    if isempty(ForceDTI) &&  isFieldSub(matName, 'dti'), fprintf('Skipping tractography (JHU DTI already computed) %s\n', imgs.ASL); return; end;
end;
if isempty(ForceDTI) &&   isFieldSub(matName, ['dti_' atlas]), fprintf('Skipping tractography (DTI already computed) %s\n', imgs.ASL); return; end;
doDtiWarpSub(imgs, atlas); %warp atlas to DTI
pth = fileparts(dti);
bed_dir=fullfile(pth, 'bedpost');
bed_dirX=fullfile(pth, 'bedpost.bedpostX');
bed_merged=fullfile(bed_dirX, 'merged');
bed_mask=fullfile(bed_dirX, 'nodif_brain_mask');
if ~exist(bed_dir,'file') || ~exist(bed_dirX,'file')
    fprintf('Please run bedpost to create files %s %s\n',bed_dir, bed_dirX);
    return;
end
dti_u=prepostfixSub('', 'du', dti);
if ~exist(dti_u,'file')
    fprintf('Can not find undistorted DTI %s\n',dti_u);
    return;
end
template_roiW=prepostfixSub('', '_roi', dti);
if exist(template_roiW,'file') && strcmpi(atlas,'jhu')
    atlasext = '';
else
   atlasext = ['_' atlas];
end
template_roiW=prepostfixSub('', ['_roi', atlasext], dti);
dti_faThr=prepostfixSub('', 'd_FA_thr', dti);
mask_dir=fullfile(pth, ['masks', atlasext]);
if ~exist(template_roiW,'file') || ~exist(dti_u,'file') || ~exist(dti_faThr,'file')
    error('doDtiTractSub Can not find %s or %s or %s\n',template_roiW, dti_u, dti_faThr);
    return;
end
template_roiWThr=prepostfixSub('', ['_roi_thr', atlasext], dti);
command=sprintf('fslmaths "%s" -mas "%s" "%s"',template_roiW, dti_faThr, template_roiWThr);
fprintf('Creating thresholded image %s\n', template_roiWThr);
doFslCmd (command);
template_roiWThr=prepostfixSub('', ['_roi_thr', atlasext], dti);
fprintf('PROBTRACKX: Create seed data\n');
if ~exist(mask_dir, 'file'), mkdir(mask_dir); end;
nROI = nRoiSub(template_roiWThr);
%now run probtrackx
prob_dir=fullfile(pth, ['probtrackx' atlasext]);
if ~exist(prob_dir, 'file'), mkdir(prob_dir); end;
nPerm = 5000; % <- arbitrary: choose your number
t_start=tic;
commands = [];
[hdr,img] = loadSub(template_roiWThr);
hdr.dt = [2 0];
hdr.pinfo = [1; 0; 0];
for i = 1: nROI
    %maski=fullfile(mask_dir, [num2str(i),'.nii.gz']);
    %command=sprintf('fslmaths "%s" -thr %d -uthr %d -bin "%s" -odt char', template_roiWThr, i,i, maski);
    %doFslCmd (command, i == 1); %only show text for 1st region
    imgi = zeros(size(img));
    imgi(img == i) = 1;
    if max(imgi) == 0, continue; end;
    maski=fullfile(mask_dir, [num2str(i),'.nii']);
    hdr.fname = maski;
    spm_write_vol(hdr,imgi);

    %command=sprintf('fslstats "%s" -M',  maski);
    %[~,cmdout] = doFslCmd (command, i == 1); %only show text for 1st region
    %if str2num(cmdout) < 1.0  %#ok<ST2NM>
    %    %delete(maski);
    %else
        prob_diri=fullfile(prob_dir, num2str(i));
        if exist(prob_diri, 'file'), rmdir(prob_diri, 's'); end;
        mkdir(prob_diri);
        if isGpuInstalledSub()
            exeName=sprintf('probtrackx2_gpu');
        else
            exeName=sprintf('probtrackx2');
        end
        command=sprintf('%s -x "%s" --dir="%s" --forcedir  -P %d -s "%s" -m "%s" --opd --pd -l -c 0.2 --distthresh=0', ...
            exeName, maski, prob_diri, nPerm ,bed_merged, bed_mask );
        commands = [commands {command}]; %#ok<AGROW>
    fprintf("%s\n", command);
    %end %if voxels survive
end %for each region
if numel(commands) < 1
   fprintf('No regions survive thresholding with FA (poor normalization?) %s', dti);
   return;
end
fprintf ('computing probtrackx for %d regions (this may take a while)\n',numel(commands) );
doThreads(commands, prob_dir);
fprintf ('probtrackx2 took %f seconds to run.\n', toc(t_start) ); %t_start=tic;
if isempty(ForceDTI)
    nii_fiber_quantify(matName, dtiDir, atlas);
else
    nii_fiber_quantify(matName, dtiDir, atlas, true);
end
%sum(nOK(:))

function [hd,im] = loadSub(fnm);
[p,n,x] = spm_fileparts(fnm);
if (length(x)==3)  && min((x=='.gz')==1)
    fnm = char(gunzip(fnm));
    delnam = fnm;
    [p,n,x] = spm_fileparts(char(fnm));
else
    delnam = '';
end
hd = spm_vol(fnm); %input header
im = spm_read_vols(hd);%Input image
hd = hd(1);
if ~isempty(delnam), delete(delnam); end;
%end loadSub()

function nROI = nRoiSub(fnm)
[hdr,img] = loadSub(fnm);
nROI = max(img(:));
%if ~spm_type(hdr.dt,'intt'), fprintf('WARNING: expected integer datatype %s\n', fnm); end;
%end nRoiSub()


function doThreads(commands, out_dir)
fslParallelSub;
command_file=fullfile(out_dir, 'commands.txt');
log_dir=fullfile(out_dir, 'logs');
if exist(log_dir, 'file'), rmdir(log_dir, 's'); end;
mkdir(log_dir);
fid = fopen(command_file, 'wt');
for i = 1 : numel(commands)
    fprintf(fid,'%s\n',commands{i});
end
fclose(fid);
command=sprintf('fsl_sub -l %s -N probtrackx -T 1 -t %s', log_dir, command_file );
doFslCmd (command, false);
%end doThreads

function fslParallelSub (maxThreads)
if ~exist('maxThreads', 'var')
    maxThreads = feature('numCores');
    if maxThreads > 15, maxThreads = maxThreads - 1; end;
end
setenv('FSLPARALLEL', num2str(maxThreads));
%end fslParallelSub()


function doDtiWarpSub(imgs, atlas)

if isempty(imgs.T1) || isempty(imgs.DTI), return; end; %required
if ~exist('atlas','var'), atlas = 'jhu'; end;
if strcmpi(atlas,'jhu')
    atlasext = '_roi';
else
   atlasext = ['_roi_' atlas];
end
T1 = prefixSub('wb',imgs.T1); %warped brain extracted image
FA = prepostfixSub('', 'd_FA', imgs.DTI);
MD = prepostfixSub('', 'd_MD', imgs.DTI);
if ~exist(T1,'file') fprintf('Unable to find image: %s\n',T1); return; end; %required
if ~exist(FA,'file') error('Unable to find image: %s\n',FA); return; end; %required
if ~exist(MD,'file') error('Unable to find image: %s\n',MD); return; end; %required
FA = unGzSub (FA);
MD = unGzSub (MD);
nFA = rescaleSub(FA);
atlasImg = fullfile(fileparts(which('NiiStat')), 'roi' , [atlas '.nii']);
if ~exist(atlasImg,'file')
    error('Unable to find template %s', atlasImg);
end
[outhdr, outimg] = nii_reslice_target(atlasImg, '', T1, 0) ;
roiname = prepostfixSub('', atlasext, imgs.DTI);
global ForceDTI;
if isempty(ForceDTI) && ~strcmpi(atlas,'jhu') &&  exist(roiname,'file')
    fprintf('Skipping doDtiWarpSub: file exists %s\n', roiname);
    return;
end
if exist(roiname,'file') %e.g. FSL made .nii.gz version
    delete(roiname);
    roiname = prepostfixSub('', atlasext, imgs.DTI);
end
roiname = unGzNameSub(roiname);
outhdr.fname = roiname;
spm_write_vol(outhdr,outimg);
oldNormSub({T1, roiname}, nFA, 8, 10, 0 );
delete(roiname);
wroiname = prepostfixSub('w', atlasext, imgs.DTI);
movefile(wroiname, roiname);
%end doFaMdSub()

function [p,n,x] = fsl_filepartsSub(imgname)
[p, n, x] = fileparts(imgname);
if strcmpi(deblank(x),'.gz') %.nii.gz
    [p, n, x] = fileparts(fullfile(p,n));
    x = [x, '.gz'];
end
%end fsl_filepartsSub()

function imgname = unGzNameSub(imgname)
[p, n, x] = fileparts(imgname);
if strcmpi(deblank(x),'.gz') %.nii.gz
    imgname = fullfile(p,n);
end
%end unGzNameSub()

function oldNormSub(src, tar, smoref, reg, interp)
%coregister T2 to match T1 image, apply to lesion
if isempty(src) || isempty(tar), return; end;
if ~exist('smoref','var'), smoref = 0; end;
if ~exist('reg','var'), reg = 1; end;
if ~exist('interp','var'), interp = 1; end;
matlabbatch{1}.spm.tools.oldnorm.estwrite.subj.source = src(1);
matlabbatch{1}.spm.tools.oldnorm.estwrite.subj.wtsrc = '';
matlabbatch{1}.spm.tools.oldnorm.estwrite.subj.resample = src(:);
matlabbatch{1}.spm.tools.oldnorm.estwrite.eoptions.template = {tar};
matlabbatch{1}.spm.tools.oldnorm.estwrite.eoptions.weight = '';
matlabbatch{1}.spm.tools.oldnorm.estwrite.eoptions.smosrc = 8;
matlabbatch{1}.spm.tools.oldnorm.estwrite.eoptions.smoref = smoref; % <-- !!!
matlabbatch{1}.spm.tools.oldnorm.estwrite.eoptions.regtype = 'mni';
matlabbatch{1}.spm.tools.oldnorm.estwrite.eoptions.cutoff = 25;
matlabbatch{1}.spm.tools.oldnorm.estwrite.eoptions.nits = 16;
matlabbatch{1}.spm.tools.oldnorm.estwrite.eoptions.reg = reg;
matlabbatch{1}.spm.tools.oldnorm.estwrite.roptions.preserve = 0;
matlabbatch{1}.spm.tools.oldnorm.estwrite.roptions.bb = [nan nan nan; nan nan nan];%[-78 -112 -70; 78 76 85];
matlabbatch{1}.spm.tools.oldnorm.estwrite.roptions.vox = [nan nan nan];%[1 1 1];
matlabbatch{1}.spm.tools.oldnorm.estwrite.roptions.interp = interp;
matlabbatch{1}.spm.tools.oldnorm.estwrite.roptions.wrap = [0 0 0];
matlabbatch{1}.spm.tools.oldnorm.estwrite.roptions.prefix = 'w';
spm_jobman('run',matlabbatch);
%end oldNormSub()


function cleanupDtiDir(imgname)
FA = prepostfixSub('', '_FA', imgname);
MD = prepostfixSub('', '_MD', imgname);
deleteImg(FA);
deleteImg(MD);
[p, n, x] = fileparts(imgname);
bDir = fullfile(p,'bedpost.bedpostX');
if exist(bDir, 'file'), rmdir(bDir, 's'); end;
bDir = fullfile(p,'bedpost');
if exist(bDir, 'file'), rmdir(bDir, 's'); end;
bDir = fullfile(p,'masks');
if exist(bDir, 'file'), rmdir(bDir, 's'); end;
bDir = fullfile(p,'probtrackx');
if exist(bDir, 'file'), rmdir(bDir, 's'); end;
%end rmBedpostDir()

function nii_check_dependencies
%make sure we can run nii_preprocess optimally
if isempty(which('NiiStat')), error('NiiStat required (https://github.com/neurolabusc/NiiStat)'); end;
%we no longer use ASLtbx:
% if isempty(which('asl_perf_subtract')), error('%s requires ASLtbx (asl_perf_subtract)\n',which(mfilename)); return; end;
% if isempty(which('spm_realign_asl')), error('%s requires ASLtbx (spm_realign_asl)\n',which(mfilename)); return; end;
if isempty(which('spm')) || ~strcmp(spm('Ver'),'SPM12'), error('SPM12 required'); end;
if isempty(which('nii_batch12')),
    p = fileparts(which('nii_preprocess'));
    p = fullfile(p,'nii_fmri');
    addpath(p);
    if isempty(which('nii_batch12'))
        error('Unable to find nii_batch12');
    end;
end;
if ~exist('spm_create_vol','file'), error('SPM12 required'); end;
if ~exist('prctile' , 'file'), error('Resting state requires "prctile" from Statistics and Machine Learning Toolbox (https://www.mathworks.com/matlabcentral/mlc-downloads/downloads/submissions/53545/versions/5/previews/tools/prctile.m/index.html).'); end;
%make sure user uncomments timing line in spm
fnm = which('spm_create_vol');
txt = fileread(fnm);
badStr = '%try, N.timing = V.private.timing; end';
if ~isempty(strfind(txt,badStr))
   error('Please uncomment "%s" from "%s"', badStr, fnm);
end
%make sure the user had the correct version of fsl_sub
fsldir = fslDirSub;
fnm = fullfile(fsldir, 'bin', 'fsl_sub');
if ~exist(fnm,'file'), error('%s required', fnm); end;
txt = fileread(fnm);
goodStr = '_gpu';
if isempty(strfind(txt,goodStr))
   error('Unless you are using SLURM, Please overwrite "%s" with file from https://github.com/neurolabusc/fsl_sub', fnm);
end
%check MRTrix installed: required for DTI
if isempty(nii_mrtrix_pth() )
   error('Please install MRTrix3'); 
end
%check FSL version
fnm = fullfile(fsldir, 'etc', 'fslversion');
fileExistsOrErrorSub(fnm);
fid = fopen(fnm);
vers = strsplit(fgetl(fid), ':'); %6.0.2:a4f562d9
fclose(fid);
if numel(vers) < 1, error('FSL version not detected: %s', fnm);  end
vers = vers{1};
versReq = '6.0.4';
if ~strcmpi(vers,versReq), error('nii_preprocess and FSL version mismatch: Expected "%s", found "%s": %s', versReq, vers, fnm); end 
%make sure eddy supports repol
cmd = fullfile(fsldir, 'bin', 'eddy_openmp');
fslEnvSub;
cmd = fslCmdSub(cmd);
[~,cmdout]  = system(cmd);
goodStr = '--repol';
if isempty(strfind(cmdout,goodStr))
   error('Update eddy files to support "repol" https://fsl.fmrib.ox.ac.uk/fsldownloads/patches/eddy-patch-fsl-5.0.9/centos6/');
end
fnm = fullfile(fsldir, 'bin', 'eddy_cuda');
if ~exist(fnm,'file')
    error('Run FSL command "configure_eddy" and test installation (https://github.com/neurolabusc/gpu_test): %s required', fnm); 
end
%requirements for BASIL, ASL
if isempty(which('nii_tool'))
    error('%s requires nii_tool for ASL (https://github.com/xiangruili/dicm2nii)\n',which(mfilename));
end
fileExistsOrErrorSub(fullfile(fsldir, 'bin', 'fsl_anat'));
fileExistsOrErrorSub(fullfile(fsldir, 'bin', 'oxford_asl'));
%end nii_check_dependencies()

function fsldir = fslDirSub;
fsldir= '/usr/local/fsl/'; %CentOS intall location
if ~exist(fsldir,'dir')
    fsldir = '/usr/share/fsl/5.0/'; %Debian install location (from neuro debian)
    if ~exist(fsldir,'dir')
        error('FSL is not in the standard CentOS or Debian locations, please check you installation!');
    end
end
%end fslDirSub()

function fileExistsOrErrorSub(fnm)
if exist(fnm,'file'), return; end;
error('Update FSL, missing required file: %s', fnm); 
%end fileExistsOrErrorSub()

% function fslEnvSub
% fsldir = fslDirSub;
% curpath = getenv('PATH');
% k = strfind(curpath, fsldir);
% if isempty(k)
% 	setenv('PATH',sprintf('%s:%s',fullfile(fsldir,'bin'),curpath));
% end
% %end fslEnvSub()

function fslEnvSub
fsldir = fslDirSub;
curpath = getenv('PATH');
k = strfind(curpath, fsldir);
if isempty(k)
	setenv('PATH',sprintf('%s:%s',fullfile(fsldir,'bin'),curpath));
end
%next code for GPUs...
ldpath = getenv('LD_LIBRARY_PATH');
if contains(ldpath,'cuda') return; end %cuda already in LD path
%look for all /usr/local/cuda* folders
pth = '/usr/local/';
d = dir([pth, '*']);
isub = [d(:).isdir];
d = {d(isub).name}';
dcuda = d(contains(d,'cuda'));
if isempty(dcuda), fprintf('GPU tools may fail: Unable to find cuda folders in %s', pth); end 
for i = 1 : numel(dcuda)
    cudalib=[[pth, dcuda{i}], filesep, 'lib64'];
    if ~exist(cudalib), continue; end
    fprintf('LD_LIBRARY_PATH now includes %s\n', cudalib); 
    setenv('LD_LIBRARY_PATH',sprintf('%s:%s',cudalib,ldpath));
end
%end fslEnvSub()

function command  = fslCmdSub (command)
fsldir = fslDirSub;
cmd=sprintf('sh -c ". %setc/fslconf/fsl.sh; ',fsldir);
command = [cmd command '"'];
%fslCmdSub

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
if status ~= 0
   warning('doFslCmd error %s\n', command); 
end
%end doFslCmd()

% function [status,cmdout]  = doFslCmd (command, verbose)
% if ~exist('verbose', 'var'),
%     verbose = true;
% end;
% fsldir= '/usr/local/fsl/';
% setenv('FSLDIR', fsldir);
% curpath = getenv('PATH');
% setenv('PATH',sprintf('%s:%s',fullfile(fsldir,'bin'),curpath));
% cmd=sprintf('sh -c ". %setc/fslconf/fsl.sh; ',fsldir);
% cmd = [cmd command '"'];
% if verbose
%     fprintf('Running \n %s\n', cmd);
% end
% [status,cmdout]  = system(cmd);
% setenv('PATH',sprintf('%s',curpath));
%end doFslCmd()


% function dtiSub(dtia)
% fsldir= '/usr/local/fsl/';
% if ~exist(fsldir,'dir'), error('%s: fsldir (%s) not found',mfilename,fsldir); end
% flirt = [fsldir 'bin/flirt'];
% if ~exist(flirt,'file')
% 	error('%s: fsl not installed (%s)',mfilename,flirt);
% end
% setenv('FSLDIR', fsldir);
% curpath = getenv('PATH');
% setenv('PATH',sprintf('%s:%s',fullfile(fsldir,'bin'),curpath));
% %./dti_travis.sh "DTIA_LM1001""
% command= [fileparts(which(mfilename)) filesep 'dti_travis2.sh'];
% command=sprintf('sh -c ". %setc/fslconf/fsl.sh; %s "%s""\n',fsldir,command, dtia);
% system(command);
% %end dtiSub()

function n = bvalCountSub(fnm)
[pth,nam, ext] = filepartsSub(fnm);
bnm = fullfile(pth, [nam, '.bval']);
vnm = fullfile(pth, [nam, '.bvec']);
if ~exist(bnm, 'file') || ~exist(vnm, 'file')
    n = 0;
    return;
end
fileID = fopen(bnm,'r');
[A, n] = fscanf(fileID,'%g'); %#ok<ASGLU>
fclose(fileID);
%check that number of elements in bval matches nifti!
if n < 2, return; end;
if strcmpi(ext,'.nii')
   h = spm_vol(fnm);
   if numel(h) ~= n, error('Number of volumes does not match "%s" "%s"', fnm, bnm); end;
end
%end bvalCountSub()

function imgs = doRestSub(imgs, matName)

if isempty(imgs.T1) || isempty(imgs.Rest), return; end; %we need these images
imgs.Rest = removeDotSub (imgs.Rest);
global ForceRest; %e.g. user can call "global ForceRest;  ForceRest = true;"
if isempty(ForceRest) && isFieldSub(matName, 'alf'), fprintf('Skipping Rest (already computed) %s\n', imgs.Rest); return; end;
delImgs('fdsw', imgs.Rest);
delImgs('fdswa', imgs.Rest);
delImgs('fdw', imgs.Rest); %old routines combine 's'mooth and 'd'etrend
delImgs('fdwa', imgs.Rest);
delMat(imgs.Rest);
%%% commented out by GY, March 4, 2017
%nii_rest(imgs);
%%% instead: (GY)
rest_prefix = nii_rest (imgs)

%7/2016 "dsw" nof "dw" as smoothing is now prior to detrending (for Chinese-style ALFF)
prefix = 'a'; %assume slice time 'a'ligned
restName = prefixSub(['dsw', prefix ],imgs.Rest);
if ~exist(restName,'file'),
    prefix = ''; %unaligned: multi-band
    restName = prefixSub(['dsw', prefix ],imgs.Rest);
    if ~exist(restName,'file')
        error('Catastrophic unspecified resting state error %s', restName);
    end
end; %required
%nii_nii2mat (prefixSub(['fdsw', prefix ],imgs.Rest), 'rest', matName)
nii_nii2mat (prefixSub(rest_prefix,imgs.Rest), 'rest', matName) % slightly modified by GY, March 3
nii_nii2mat (prefixSub(['palf_dsw', prefix ],imgs.Rest), 'alf', matName) %detrended
nii_nii2mat (prefixSub(['palf_sw', prefix ],imgs.Rest), 'palf', matName) %conventional linear trends only
vox2mat(prefixSub(['mean' ],imgs.Rest), 'NativeRestAve', matName); %no prefix: prior to slice time
vox2mat(prefixSub(['wmean' ],imgs.Rest), 'RestAve', matName); %no prefix: prior to slice time
%vox2mat(prefixSub(['wbmean' ],imgs.Rest), 'RestAve', matName);
vox2mat(prefixSub(['wbmaskmean' ],imgs.Rest), 'RestAve', matName); % Grigori added "mask"

%end doRestSub()

function delImgs(prefix, fnm)
%e.g. delImgs('sw', 'X.nii') would delete wX.nii and swX.nii
[pth,nam] = spm_fileparts(fnm);
fnmMat = fullfile(pth,[nam,'.mat']);
for i = 1: numel(prefix),
    p = prefix(end-i+1:end);
    d = prefixSub(p, fnm);
    if exist(d,'file'), delete(d); end; %delete image
    d = prefixSub(p, fnmMat);
    if exist(d,'file'), delete(d); end; %delete lmat
end;
%end delImgs()

function delMat(fnm)
[pth,nam] = spm_fileparts(fnm);
fnmMat = fullfile(pth,[nam,'.mat']);
if exist(fnmMat,'file'), delete(fnmMat); end; %delete lmat
%end delMat()

function idx = roiIndexSub(roiName)
[~, ~, idx] = nii_roi_list(roiName);
if idx < 1, error('Invalid roi name %s', roiName); end;
%end roiIndexSub()

function imgs = doAslSubOld(imgs, matName)
if isempty(imgs.T1) || isempty(imgs.ASL), return; end; %we need these images
imgs.ASL = removeDotSub (imgs.ASL);
global ForceASL; %e.g. user can call "global ForceASL;  ForceASL = true;"
if isempty(ForceASL) && isFieldSub(matName, 'cbf'), fprintf('Skipping ASL (CBF already computed) %s\n', imgs.ASL); return; end;
delImgs('src', imgs.ASL); %remove previously preprocessed data
delMat(imgs.ASL);
[nV, nSlices] = nVolSub (imgs.ASL) ;
[mx, ind] = max(nV);
if (mod(mx,2) == 0) && ( (nSlices ~= 17) && (nSlices ~= 16))
	fprintf('Error: nii_pasl12 only designed for pCASL sequences with 17 or 16 slices not %d: %s\n', nSlices);
	return;
end
if mx < 60,
    fprintf('not enough ASL volumes (only %d volumes with %d slices) for %s\n', nV, nSlices, imgs.ASL);
    imgs.ASL =[];
    return;
end;
asl = imgs.ASL(ind,:);
if ~exist(prefixSub('b',imgs.T1),'file'), return; end; %required
if isempty(ForceASL) && exist(prefixSub('wmeanCBF_0_src',asl),'file'), fprintf('Skipping ASL (wmeanCBF already computed) %s\n', imgs.ASL); return; end; %already computed
[cbf, c1L, c1R, c2L, c2R] = nii_pasl12(asl, imgs.T1);
if isempty(cbf), warning('Catastrophic error with nii_pasl12'); return; end;
nii_nii2mat(cbf, 'cbf', matName); %2
stat = load(matName);
stat.cbf.nV = nV;
stat.cbf.c1L = c1L;
stat.cbf.c1R = c1R;
stat.cbf.c2L = c2L;
stat.cbf.c2R = c2R;
save(matName,'-struct', 'stat');
%end doAslSubOld()


function imgs = doAslSub(imgs, matName)
if isempty(imgs.T1) || isempty(imgs.ASL), return; end;
imgs.ASL = removeDotSub (imgs.ASL);
global ForceASL;
inCalScan = []; %M0 scan
anatDir = '';
dryRun=false;
basilDir = fullfile(fileparts(imgs.T1), 'basil');
if isempty(ForceASL) && isFieldSub(matName, 'basilStd')
    fprintf('Skipping ASL (Basil already computed) %s\n', imgs.ASL);
    return;
end;
if exist(basilDir,'file'), rmdir(basilDir,'s'); end; %delete ASL directory
if ~exist('basilDir', 'dir'), mkdir(basilDir); end;
[p,n] = fileparts(imgs.ASL);
inJSON = fullfile(p, [n, '.json']);
if ~exist(inJSON, 'file')
   error('Unable to find %s', inJSON); 
end
T1 = imgs.T1;
[p,n,x] = fileparts(T1);
healedT1 = fullfile(p,['e',n,x]);
if exist(healedT1, 'file') 
    T1 = healedT1;
else
    error('Unable to find %s\n', healedT1);
end
exitCode = nii_basil(imgs.ASL, T1, imgs.ASLrev, inJSON, inCalScan, anatDir, dryRun, basilDir);
if exitCode ~= 0, return; end;
basil2mat(basilDir, matName);
%end doAslSub()

function imgs = doAslSubOldish(imgs, matName)
if isempty(imgs.T1) || isempty(imgs.ASL), return; end;
imgs.ASL = removeDotSub (imgs.ASL);
global ForceASL;
ASLrev = []; %reverse phase encoded image
jsonFile = [];   %name of json file that stores parameters for this scan
[numVolASLImg, ~] = nVolSub (imgs.ASL) ;
switch numVolASLImg
    case 101 %handles POLAR and legacy pasl sequences
        jsonFile = which('POLAR_dummy.json');
    case 74
        jsonFile = which('LEGACY_PCASL_dummy.json');
    case 97
        jsonFile = which('LARC_dummy.json'); 
        [pth,nam,ext] = spm_fileparts(imgs.ASL);
        ASLrev = [pth '/' nam 'rev',ext];
    case 60
        jsonFile = which('SEN_dummy.json'); 
    otherwise
        error('Can''t guess sequence type');
end
pth = spm_fileparts(imgs.ASL);
basilDir = fullfile(pth, 'BASIL');
if ~exist('basilDir', 'dir'), mkdir(basilDir); end;
[filepath,name,ext] = fileparts(imgs.T1)
healedT1 = [filepath,'/e',name,ext];
%function exitCode = nii_basil(asl, t1, aslRev, inJSON, inCalScan, anatDir, dryRun, overwrite, outDir)
exitCode = nii_basil(imgs.ASL,healedT1, ASLrev, jsonFile, '', '', false, basilDir)
if exitCode ~= 0, error('BASIL failed\n'); end;
basil2mat(basilDir, matName);
%end doAslSub()


function [nVol, nSlices] = nVolSub (fnm)
%Report number of volumes
% v= nVols('img.nii');
% v= nVolS(strvcat('ASL_P005.nii', 'ASL_P005_1.nii'))
% v= nVolS({'ASL_P005.nii', 'ASL_P005_1.nii'})
nVol = [];
nSlices = 0;
fnm = cellstr(fnm);
for v = 1 : numel(fnm)
    hdr = spm_vol(deblank(char(fnm(v,:))));
    nVol = [nVol numel(hdr)]; %#ok<AGROW>
    nSlices = hdr(1).dim(3);
end
%end nVol()

function imgs = doI3MSub(imgs, matName)
if isempty(imgs.T1), return; end; %we need these images
imgs.T1 = removeDotSub (imgs.T1);
if isFieldSub(matName, 'i3mT1'), fprintf('Skipping i3m (already computed) %s\n', imgs.ASL); return; end;
w = prefixSub('w',imgs.T1); %warped T1
if  ~exist(w, 'file')  return; end; %exit: we require warped T1
i3m = prefixSub('zw',imgs.T1);
if isFieldSub(matName, 'i3mT1'), return; end;
nii_i3m(w,'',0.5,10,0.25,1); %i3m T1 image
nii_nii2mat(i3m, 'i3mT1', matName); %4
%end doI3MSub()

function imgs = doT1Sub(imgs, matName)
if isempty(imgs.T1), return; end;
if isFieldSub(matName, 'T1'), fprintf('Skipping T1 normalization (already done): %s\n',imgs.T1); return; end;
imgs.T1 = removeDotSub (imgs.T1);
if size(imgs.T1,1) > 1 || size(imgs.T2,1) > 1 || size(imgs.Lesion,1) > 1
    error('Require no more than one image for these modalities: T1, T2, lesion');
end;
%we have a separate version of nii_enat_norm in SPM scripts as well
nii_enat_norm_npp(imgs.T1,imgs.Lesion,imgs.T2);
if ~isFieldSub(matName, 'T1')
    bT1 = prefixSub('wb',imgs.T1);
    vox2mat(bT1,'T1',matName);
    %nii_nii2mat(bT1,'T1',matName);
end
if isempty(imgs.Lesion), return; end;
if isFieldSub(matName, 'lesion'), return; end;
wr = prefixSub('wr',imgs.Lesion);
if ~exist(wr,'file'), wr = prefixSub('ws',imgs.Lesion); end; %T1 but no T2
if ~exist(wr,'file'), wr = prefixSub('wsr',imgs.Lesion); end; %T1 and T2, smoothed
if ~exist(wr,'file'), error('Unable to find %s', wr); end;
nii_nii2mat(wr, 'lesion', matName); %1
%end doT1Sub()

function vox2mat(imgName, fieldName, matName)
%save voxelwise data as field 'fieldName' in the .mat file matName
%  similar to nii_nii2mat, but does not compute ROI values
%Example
% vox2mat('wbT1.nii', 'T1', 'M2012.mat');
if ~exist(imgName, 'file'), fprintf('vox2mat Unable to find %s\n', imgName); return; end;
hdr = spm_vol(imgName);
img = spm_read_vols(hdr);
if ndims(img) ~= 3, fprintf('Not a 3D volume %s\n', imgName); return; end;
stat = [];
stat.(fieldName).hdr = hdr;
stat.(fieldName).dat = img;
if exist(matName,'file')
    old = load(matName);
    stat = nii_mergestruct(stat,old); %#ok<NASGU>
end
save(matName,'-struct', 'stat');
fprintf('Saving %s in %s\n', fieldName, matName);
%end vox2mat()

function bt1 = normSub(t1)
template = fullfile(spm('Dir'),'tpm','TPM.nii');
if ~exist(template,'file')
    error('Unable to find template named %s',template);
end
if iscell(t1), t1 = char(t1); end;
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
bt1 = extractSub(0.005, t1, prefixSub('c1', t1), prefixSub('c2', t1));
%end normSub();

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


function is = isFieldSub(matname, fieldname)
is = false;
if ~exist(matname, 'file'), return; end;
m = load(matname);
is = isfield(m, fieldname);
%end isFieldSub

function nam = prefixSub (pre, nam)
if isempty(nam), return; end;
[p, n, x] = spm_fileparts(nam);
nam = fullfile(p, [pre, n, x]);
%end prefixSub()

function imgs = unGzAllSub(imgs)
imgs.ASL = unGzSub(imgs.ASL);
%imgs.DTI = unGzSub(imgs.DTI); %fsl is fine with gz
%imgs.DTIrev = unGzSub(imgs.DTIrev); %fsl is fine with gz
imgs.fMRI = unGzSub(imgs.fMRI);
imgs.Lesion = unGzSub(imgs.Lesion);
imgs.Rest = unGzSub(imgs.Rest);
imgs.T1 = unGzSub(imgs.T1);
imgs.T2 = unGzSub(imgs.T2);
%end subjStructSub()

function imgs = removeDotDtiSub(imgs)
imgs.DTI = removeDotSub(imgs.DTI); %fsl is fine with gz
imgs.DTIrev = removeDotSub(imgs.DTIrev); %fsl is fine with gz
%end subjStructSub()

function fnm = removeDotSub (fnm)
%topup will make 'my.img.nii.gz -> my.topup and my.img.topup.nii.gz
if isempty(fnm), return; end; %CR210116 - e.g. if DTI exists but DTIrev is empty
[p,n, x] = filepartsSub(fnm);
if isempty(strfind(n,'.')), return; end;
nn = strrep(n, '.', '_');
oldfnm = fnm;
fnm = fullfile(p, [nn, x]);
if exist(oldfnm, 'file') %if a script is re-run, return the corrected name
    movefile(oldfnm, fnm);
end;
%Taylor Hanayik added for .bvec and .bval 'un'dotting
%if isDTI
    if exist(fullfile(p,[n, '.bvec']), 'file')
        movefile(fullfile(p,[n, '.bvec']), fullfile(p,[nn, '.bvec']));
    end
    if exist(fullfile(p,[n, '.bval']), 'file')
        movefile(fullfile(p,[n, '.bval']), fullfile(p,[nn, '.bval']));
    end
%end
%end removeDotSub

function fnm = unGzSub (fnm)
if isempty(fnm), return; end;
innames = fnm;
fnm = [];
for i = 1: size(innames,1)
    f = deblank(innames(i,:));
    f = unGzCSub(f);
    fnm = strvcat(fnm, f);
end
%end unGzSub()

function deleteImg(fnm)
if isempty(fnm), return; end;
fnm = deblank(fnm);
[pth,nam,ext] = spm_fileparts(fnm);
if strcmpi(ext,'.gz') %.nii.gz
    fnm = fullfile(pth, nam);
    [pth,nam,ext] = spm_fileparts(fnm);
end
if strcmpi(ext,'.nii') %.nii.gz
    fnm = fullfile(pth, nam);
end
fnmN = [fnm, '.nii'];
if exist(fnmN,'file')
  delete(fnmN);
end
fnmZ = [fnm, '.nii.gz'];
if exist(fnmZ,'file')
  delete(fnmZ);
end
%
function [fnm, isGz] = unGzCSub (fnm)
if isempty(fnm), return; end;
fnm = deblank(fnm);
[pth,nam,ext] = spm_fileparts(fnm);
isGz = false;
if strcmpi(ext,'.nii') % fsl can not handle .nii coexisting with .nii.gz
    fnmz = fullfile(pth, [nam ext '.gz']);
    if exist(fnmz, 'file') && exist(fnm, 'file')
        delete(fnmz); %remove .nii.gz if .nii exists
    elseif exist(fnmz, 'file') %img.nii.gz exists but img.nii does not!
        fnm = fnmz;
        [pth,nam,ext] = spm_fileparts(fnm);
    end
end
if ~exist(fnm, 'file'), return; end;
if strcmpi(ext,'.gz') %.nii.gz
    ofnm = fnm;
    fnm = char(gunzip(fnm));
    if isunix && ~ismac
        if ~exist(fnm,'file') %current working directory?
            fprintf('#######  Taylor''s special code executed "%s" "%s"\n', pth, fnm);
            fnm = fullfile(pth,fnm); %Taylor Hanayik added on 08Sep2015 to fix error in processing
        end;
        %fnm = fullfile(pth,fnm); %Taylor Hanayik added on 08Sep2015 to fix error in processing
    end
    isGz = true;
    delete(ofnm);
elseif strcmpi(ext,'.voi') %.voi ->
    onam = char(gunzip(fnm));
    fnm = fullfile(pth, [nam '.nii']);
    movefile(onam,fnm);
    isGz = true;
end;
%end unGzSub()


function imgs = subjStructSub(subjName)
imgs.name = subjName;
imgs.ASL = '';
imgs.DTI = '';
imgs.DTIrev = '';
imgs.fMRI = '';
imgs.Lesion = '';
imgs.Rest = '';
imgs.T1 = '';
imgs.T2 = '';
%end subjStructSub()

function imgName = imgfindSub (imgName, imgKey, inDir)
%look for a filename that includes imgKey in folder inDir or subfolders
% for example if imgKey is 'T1' then T1 must be in both folder and image name myFolder\T1\T1.nii
%if ~isempty(imgName), return; end;
[pth, nam] = fileparts(inDir); %#ok<ASGLU> %e.g. 'T1folder' for /test/T1folder
nameFiles = subFileSub(inDir);
nameFiles = sort(nameFiles); %take first file for multiimage sequences, e.g. ASL
for i=1:size(nameFiles,1)
    pos = isStringInKey (nameFiles(i), imgKey);
    if pos == 1 && isImgSub(char(nameFiles(i)))
        imgName = strvcat(imgName, [inDir, filesep, char(nameFiles(i))]);

    end; %do not worry about bvec/bval
end
if isempty(imgName), fprintf('WARNING: unable to find any "%s" images in folder %s\n',deblank(imgKey(1,:)), inDir); end;
if size(imgName) > 1, imgName = deblank(imgName(1,:)); end;
%end imgfindSub()

function isKey = isStringInKey (str, imgKey)
isKey = true;
for k = 1 : size(imgKey,1)
    key = deblank(imgKey(k,:));
    pos = strfind(lower(char(str)),lower(key));
    if ~isempty(pos), isKey = pos(1); return; end;
end
isKey = false;
%isStringInKey()

% function imgName = imgfindSub (imgName, imgKey, inDir, imgKey2)
% %look for a filename that includes imgKey in folder inDir or subfolders
% % for example if imgKey is 'T1' then T1 must be in both folder and image name myFolder\T1\T1.nii
% if ~exist('imgKey2','var'), imgKey2 = imgKey; end;
% if ~isempty(imgName), return; end;
% [pth, nam] = fileparts(inDir); %#ok<ASGLU> %e.g. 'T1folder' for /test/T1folder
% if isempty(strfind(lower(char(nam)), lower(imgKey))), return; end;
% if exist([inDir,filesep, 'Native'],'file')
%     inDir = [inDir,filesep, 'Native'];
%     %fprintf('xxxx %s\n', inDir);
% end
% nameFiles = subFileSub(inDir);
% nameFiles = sort(nameFiles); %take first file for multiimage sequences, e.g. ASL
% for i=1:size(nameFiles,1)
%     pos = strfind(lower(char(nameFiles(i))),lower(imgKey));
%     if isempty(pos)
%         pos = strfind(lower(char(nameFiles(i))),lower(imgKey2));
%     end
%     if ~isempty(pos) && isImgSub(char(nameFiles(i)))
%         imgName = [inDir, filesep, char(nameFiles(i))];
%         return
%     end; %do not worry about bvec/bval
% end
% fprintf('WARNING: unable to find any "%s" images in folder %s\n',imgKey, inDir);
% %end imgfindSub

function nameFiles=subFileSub(pathFolder)
d = dir(pathFolder);
isub = ~[d(:).isdir];
nameFiles = {d(isub).name}';
%end subFileSub()

function nameFolds=subFolderSub(pathFolder)
d = dir(pathFolder);
isub = [d(:).isdir];
nameFolds = {d(isub).name}';
nameFolds(ismember(nameFolds,{'.','..'})) = [];
%end subFolderSub()

function isImg = isImgSub (fnm)
[pth,nam,ext] = spm_fileparts(fnm); %#ok<ASGLU>
isImg = false;
if strcmpi(ext,'.gz') || strcmpi(ext,'.voi') || strcmpi(ext,'.hdr') || strcmpi(ext,'.nii')
    isImg = true;
end;
%end isImgSub()

function [pth,nam,ext,num] = filepartsSub(fname)
% extends John Ashburner's spm_fileparts.m to include '.nii.gz' as ext
num = '';
if ~ispc, fname = strrep(fname,'\',filesep); end
[pth,nam,ext] = fileparts(deblank(fname));
ind = find(ext==',');
if ~isempty(ind)
    num = ext(ind(1):end);
    ext = ext(1:(ind(1)-1));
end
if strcmpi(ext,'.gz')
   [pth nam ext] = fileparts(fullfile(pth, nam));
   ext = [ext, '.gz'];
end
%end filepartsSub()

function v = isGpuInstalledSub()
v = ~isempty(strfind(getenv('PATH'),'cuda'));
%v = 0; %default to 0
%if ~exist('/usr/local/cuda/bin/nvcc','file')
%    return;
%end
%if nargin < 1
%    vstr = 'release 6.5';
%end
%[s, thisVer] = system(['/usr/local/cuda/bin/nvcc --version | grep -o "' vstr '"']);
%v = strcmpi(deblank(thisVer),vstr);
%end isGpuInstalledSub()

function isHalfSphere = HalfSphere(bvec_nam)
%return true if DWI vectors sample half sphere
% bvec_nam : filename of b-vector file
%Example
% HalfSphere('dki.bvec')
if ~exist('bvec_nam','var') %file not specified
   fprintf('Select b-vec file\n');
   [A,Apth] = uigetfile({'*.bvec'},'Select b-vec file');
   if isnumeric(A), return; end; %nothing selected
   bvec_nam = strcat(Apth,char(A));
end;
% %handle different extensions bvec/nii/nii.gz
[p,n,x]=fileparts(bvec_nam);
if ~strcmpi(x,'.bvec')
    fnm = fullfile(p,[n,'.bvec']);
    if ~exist(fnm,'file') %assume img.nii.gz
        [~,n,x] = fileparts(n);
        fnm = fullfile(p,[n,'.bvec']);
    end
    if ~exist(fnm,'file')
        error('Unable to find bvec file for %s', vNam);
    end
    bvec_nam = fnm;
end
isHalfSphere = false;
bvecs = load(bvec_nam); 
bvecs = unitLengthSub(bvecs)';
bvecs(~any(~isnan(bvecs')),:) = [];
if isempty(bvecs), return; end;
mn = unitLengthSub(mean(bvecs)')';
if isnan(mn), return; end;
mn = repmat(mn,size(bvecs,1),1);
minV = min(dot(bvecs',mn'));
thetaInDegrees = acosd(minV);
if thetaInDegrees < 110
    isHalfSphere = true;
    fprintf('Sampling appears to be half sphere (%g-degrees): %s\n', thetaInDegrees, bvec_nam);
end;
%end HalfSphere()

function Xu = unitLengthSub(X) 
%Use pythagorean theorem to set all vectors to have length 1
%does not require statistical toolbox 
if size(X,1) ~= 3, error('expected each column to have three rows (X,Y,Z coordinates)'); end;
Dx =  sqrt(sum((X.^2)));
Dx = [Dx; Dx; Dx];
Xu = X./Dx;
%end unitLengthSub()

% function [isFullSphere, meanLength] = isFullSphereSub(vNam) %does bvec sample whole sphere?
% %given FSL style bvec, determine if sampling over full or half sphere
% %useful for determining if one should run eddy (requires full sphere) or eddy_correct
% %see nii_plotBvec for visualization
% % vNam : name of bvec file
% %Output
% % isFullSphere: returns true if mean vector near zero (vectors cancel each other)
% % meanLength: length of mean vector, near zero for full sphere
% %Examples
% % isFull = nii_meanBvec('a.bvec');
% % isFull = nii_meanBvec(); %use GUI
% if ~exist('vNam','var') %no files
%  vNam = spm_select(1,'^.*\.(bvec)$','Select bvec file to view');
% end;
% %handle different extensions bvec/nii/nii.gz
% [p,n,x]=fileparts(vNam);
% if ~strcmpi(x,'.bvec')
%     fnm = fullfile(p,[n,'.bvec']);
%     if ~exist(fnm,'file') %assume img.nii.gz
%         [~,n,x] = fileparts(n);
%         fnm = fullfile(p,[n,'.bvec']);
%     end
%     if ~exist(fnm,'file')
%         error('Unable to find bvec file for %s', vNam);
%     end
%     vNam = fnm;
% end
% if ~exist(vNam, 'file'), error('Unable to find file %s',vNam); end;
% %read input
% v = textread(vNam);
% v( :, all( ~any( v ), 1 ) ) = []; %delete vectors with all zeros (e.g. B=0)
% v = mean(v,2); %mean vector
% meanLength = sqrt(sum((v) .^ 2)); %mean vector length
% %if sperical then vectors cancel out and mean length near zero
% %if half-sphere than the "center of mass" for vectors is biased
% isFullSphere = (meanLength < 0.25);
% %end isFullSphereSub()




















