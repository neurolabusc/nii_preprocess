function matName = nii_preprocess(imgs, matName)
%preprocess data from multiple modalities3
% imgs.T1: filename of T1 scan - ONLY REQUIRED INPUT: ALL OTHERS OPTIONAL
% imgs.T2: filename used to draw lesion, if absent lesion drawn on T1
% imgs.lesion : lesion map
% imgs.ASL : pCASL or PASL sequence
% imgs.Rest : Resting state sequence
% DTI : Diffusion scan
% DTIrev : Diffusion scan with reverse phase encoding of 'DTI'
% fMRI : fMRI scan
%Examples
% imgs.T1 = 'T1.nii'; imgs.ASL = 'ASL.nii';
% nii_preprocess(imgs);
%check dependencies

fprintf('%s version 27July2016\n', mfilename);
checkForUpdate(fileparts(mfilename('fullpath')));
if nargin < 1, error('Please use nii_preprocess_gui to select images'); end;
if isempty(which('NiiStat')), error('NiiStat required'); end;
if isempty(which('spm')) || ~strcmp(spm('Ver'),'SPM12'), error('SPM12 required'); end;
if isempty(spm_figure('FindWin','Graphics')), spm fmri; end; %launch SPM if it is not running
%f = spm_figure('FindWin','Graphics'); clf(f.Number); %clear SPM window

%set structures
if ~isfield(imgs,'T1') || isempty(imgs.T1), error('T1 scan is required'); end;
if ~isfield(imgs,'T2'), imgs.T2 = []; end;
if ~isfield(imgs,'Lesion'), imgs.Lesion = []; end;
if ~isfield(imgs,'ASL'), imgs.ASL = []; end;
if ~isfield(imgs,'Rest'), imgs.Rest = []; end;
if ~isfield(imgs,'DTI'), imgs.DTI = []; end;
if ~isfield(imgs,'DTIrev'), imgs.DTIrev = []; end;
if ~isfield(imgs,'fMRI'), imgs.fMRI = []; end;

if ~exist('matName','var')
    [p,n] = filepartsSub(imgs.T1);
    matName = fullfile(p, [n, '_lime.mat']);
end

if true
    diary([matName, '.log.txt'])
    imgs = unGzAllSub(imgs); %all except DTI - fsl is OK with nii.gz
    tStart = tic;
    imgs = doT1Sub(imgs, matName); %normalize T1
    imgs = doI3MSub(imgs, matName);
    tStart = timeSub(tStart,'T1');
    imgs = doRestSub(imgs, matName); %TR= 1.850 sec, descending; %doRestSub(imgs, matName, 2.05, 5); %Souvik study
    tStart = timeSub(tStart,'REST');
    imgs = doAslSub(imgs, matName);
    tStart = timeSub(tStart,'ASL');
    imgs = dofMRISub(imgs, matName);
    tStart = timeSub(tStart,'fMRI');
    if ~isempty(imgs.DTI)
        imgs = removeDotDtiSub(imgs);
        dtiDir = fileparts(imgs.DTI);
        doDtiSub(imgs);
        %-->(un)comment next line for JHU tractography
        doDtiTractSub(imgs, matName, dtiDir, 'jhu');
        %-->(un)comment next line for AICHA tractography
        %doDtiTractSub(imgs, matName, dtiDir, 'AICHA')
        %-->compute scalar DTI metrics
        doFaMdSub(imgs, matName);
        doTractographySub(imgs);
        doDkiSub(imgs, matName);
        tStart = timeSub(tStart,'DTI');
    end
    doDkiSub(imgs, matName, true);
    tStart = timeSub(tStart,'DKI');
end
%print output
pth = '/home/crlab/Desktop';
if ~exist(pth,'file'), pth = ''; end;
printDTISub(imgs, fullfile(pth,'MasterDTI')); %show results - DTI
nii_mat2ortho(matName, fullfile(pth,'MasterNormalized')); %do after printDTI (spm_clf) show results - except DTI
diary off
%nii_preprocess()

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
FA = prepostfixSub('n', '_FA', imgs.DTI);
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

function imgs = dofMRISub(imgs, matName)
if isempty(imgs.T1) || isempty(imgs.fMRI), return; end; %required
imgs.fMRI = removeDotSub (imgs.fMRI);
[~,n] = fileparts(imgs.fMRI);
[p] = fileparts(imgs.fMRI);
cstat = fullfile(p,[n], 'con_0002.nii');
bstat = fullfile(p,[n], 'beta_0001.nii');
global ForcefMRI; %e.g. user can call "global ForcefMRI;  ForcefMRI = true;"
if isempty(ForcefMRI) && isFieldSub(matName, 'fmri') && exist(cstat, 'file') && exist(bstat,'file'), fprintf('Skipping fMRI (already done) %s\n',imgs.fMRI); return;  end;
%if ~isempty(ForcefMRI)
    d = fullfile(p,n);
    if exist(d,'file'), rmdir(d,'s'); end; %delete statistics directory
    delImgs('sw', imgs.fMRI);
    %delImgs('swa', imgs.fMRI); %we never slice-time correct sparse data!
%end
if ~exist('nii_fmri60.m','file')
    fnm = fullfile(fileparts(which(mfilename)), 'nii_fmri');
    if ~exist(fnm,'file')
        error('Unable to find %s', fnm);
    end
    addpath(fnm);
end
nii_fmri60(imgs.fMRI, imgs.T1, imgs.T2); %use fMRI for normalization
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
Mask = prepostfixSub('', '_FA_thr', imgs.DTI);
if ~exist(Mask,'file') , return; end; %required
[p,n] = fsl_filepartsSub(imgs.DTI);
basename = fullfile(p,n);
vtkname = [basename, '.vtk'];
if exist(vtkname,'file') , fprintf('Skipping tractography (alredy done): %s\n', vtkname); return; end; %required
imgCnv = [basename, '_dtitk.nii.gz'];
ConvExe = 'TVFromEigenSystem';
PathConvExe = findExeSub(ConvExe);
if isempty(PathConvExe), fprintf('Tractography skipped: unable to find %s', ConvExe); return; end;
TrakExe = 'SingleTensorFT';
PathTrakExe = findExeSub(TrakExe);
if isempty(PathTrakExe), fprintf('Tractography skipped: unable to find %s', TrakExe); return; end;
cmd = sprintf('%s -basename "%s" -out %s -type FSL',ConvExe,basename, imgCnv);
[status, fullnam]  = system(cmd,'-echo');
%e.g. TVFromEigenSystem -basename 199_99_AP_7 -type FSL
cmd = sprintf('%s -in "%s" -seed "%s" -out "%s"',TrakExe, imgCnv, Mask, vtkname );
%e.g. SingleTensorFT -in 199_99_AP_7.nii.gz -seed 199_99_AP_7_FA_thr.nii.gz -out t.vtk
[status, fullnam]  = system(cmd,'-echo');
%end doTractographySub()

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
    if isFieldSub(matName, 'mk'), return; end; %stats already exist
    if isEddyCuda7Sub()
         command= [fileparts(which(mfilename)) filesep 'dti_1_eddy_cuda.sh'];
    else
        command= [fileparts(which(mfilename)) filesep 'dti_1_eddy.sh'];
    end
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
if ~exist('dkifx','file'),  fprintf('skipping DKI: requires dkifx script\n'); return; end;
mask=prepostfixSub('', 'b_mask', DTI);
dti_u=prepostfixSub('', 'u', DTI);
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
dkifx(dti_u, bvalnm, mask);
MKmask = prepostfixSub('', 'u_ldfDKI_MASK', DTI);
MK=prepostfixSub('', 'u_ldfDKI_MK', DTI);
if ~exist(MK, 'file') || ~exist(MKmask, 'file')
    fprintf('Serious error: no kurtosis images named %s %s\n', MK, MKmask);
    return;
end;
%normalize mean kurtosis to match normalized, brain extracted T1
wMK = prepostfixSub('w', '', MK);
oldNormSub( {MKmask, MK}, wbT1, 8, 8 );
nii_nii2mat(wMK, 'mk', matName);
%save note
fid = fopen('dki.txt', 'a+');
fprintf(fid, '%s\n', matName);
fclose(fid);
%end doDkiCoreSub()


function doFaMdSub(imgs, matName)
if isempty(imgs.T1) || isempty(imgs.DTI), return; end; %required
T1 = prefixSub('wb',imgs.T1); %warped brain extracted image
FA = prepostfixSub('', '_FA', imgs.DTI);
MD = prepostfixSub('', '_MD', imgs.DTI);
if ~exist(T1,'file') || ~exist(FA,'file') || ~exist(MD,'file') , return; end; %required
FA = unGzSub (FA);
MD = unGzSub (MD);
wFA = prepostfixSub('w', '', FA);
wMD = prepostfixSub('w', '', MD);
if ~exist(wFA,'file') || ~exist(wMD,'file')
    nFA = rescaleSub(FA);
    %nii_setOrigin12({nFA, FA, MD}, 1,false); %rFA looks like T1
    oldNormSub( {nFA, FA,MD}, T1, 8, 10 );
end
if isFieldSub(matName, 'fa'), return; end; %stats already exist
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
if isDtiDone(imgs), fprintf('Skipping DTI processing (probtrackx done)\n'); return; end;
n = bvalCountSub(imgs.DTI);
if (n < 1)
    fprintf('UNABLE TO FIND BVECS/BVALS FOR %s\n', imgs.DTI);
    return
end
if (n < 12)
    fprintf('INSUFFICIENT BVECS/BVALS FOR %s\n', imgs.DTI);
    return
end
dti_u=prepostfixSub('', 'u', imgs.DTI);
if exist(dti_u, 'file')
    fprintf('Skipping DTI preprocessing: found %s\n', dti_u);
else
    clipSub (imgs); %topup requires images with even dimensions
    %1 - eddy current correct
    if isEddyCuda7Sub()
         command= [fileparts(which(mfilename)) filesep 'dti_1_eddy_cuda.sh'];
    else
        command= [fileparts(which(mfilename)) filesep 'dti_1_eddy.sh'];
    end
    if isempty(imgs.DTIrev)
        command=sprintf('%s "%s"',command, imgs.DTI);
    else
        nr = bvalCountSub(imgs.DTIrev);
        if (nr ~= n)
            fprintf('WARNING: BVECS/BVALS DO NOT MATCH %s %s\n', imgs.DTI, imgs.DTIrev);
            %return
        end
        command=sprintf('%s "%s" "%s"',command, imgs.DTI, imgs.DTIrev);
    end
    cleanupDtiDir(imgs.DTI);
    doFslCmd (command);
end
doDtiBedpostSub(imgs.DTI);
%end doDtiSub()

function isEddyCuda = isEddyCuda7Sub
%eddy_cuda7.0 is a beta version supplied by Jesper Andersson that supports the GTX970 and has new features
isEddyCuda = false;
if ~isGpuInstalledSub(), return; end;
eddyName = '/usr/local/fsl/bin/eddy_cuda7.0';
isEddyCuda = exist(eddyName,'file') > 0;
if ~isEddyCuda, printf('Hint: Eddy will run faster if you install %s', eddyName); end;
%end isEddyCuda7Sub

function done = isDtiDone(imgs)
done = false;
if isempty(imgs.DTI), return; end;
p = fileparts( imgs.DTI );
pDir = fullfile(p,'probtrackx');
if exist(pDir, 'file')
    done = true;
end;
%end isDtiDone()

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
%if exist(bed_dirX, 'file'), rmdir(bed_dirX, 's'); end; %666
bed_done=fullfile(bed_dirX, 'xfms', 'eye.mat');
if exist(bed_done,'file'), fprintf('Skipping bedpost (already done)\n'), return, end;
if ~exist(bed_dirX, 'file'), mkdir(bed_dirX); end;
bed_dir=fullfile(pth, 'bedpost');
%if exist(bed_dir, 'file'), rmdir(bed_dir, 's'); end; %666
if ~exist(bed_dir, 'file'), mkdir(bed_dir); end;
dti_u=prepostfixSub('', 'u', dti);
dti_x=fullfile(bed_dir, 'data.nii.gz');
if ~exist(dti_u,'file'), error('Bedpost unable to find %s', dti_u); end;
copyfile(dti_u, dti_x);
[bvec, bval] = getBVec(dti);
dti_x=fullfile(bed_dir, 'bvecs');
copyfile(bvec, dti_x);
dti_x=fullfile(bed_dir, 'bvals');
copyfile(bval, dti_x);
dti_faThr=prepostfixSub('', '_FA_thr', dti);
dti_x=fullfile(bed_dir, 'nodif_brain_mask.nii.gz');
copyfile(dti_faThr, dti_x);
if isGpuInstalledSub
    command=sprintf('bedpostx_gpu "%s" ', bed_dir);
else
    command=sprintf('bedpostx "%s" ', bed_dir);
end
fslParallelSub;
doFslCmd (command);
while ~exist(bed_done,'file')
    pause(1.0);
end
fprintf ('Bedpost took %f seconds to run.\n', toc(t_start) );
%end doDtiBedpostSub()

function [bvec, bval] = getBVec(dti)
%both.bval
dti = unGzNameSub(deblank(dti));
[p,n] = fileparts(dti);
dti = fullfile(p,n);
bvec = [dti 'both.bvec'];
bval = [dti 'both.bval'];
if exist(bvec,'file') && exist(bval,'file'), return; end;
bvec = [dti '.bvec'];
bval = [dti '.bval'];
if ~exist(bvec,'file') || ~exist(bval,'file'), error('Can not find files %s %s', bvec, bval); end;
%end getBVec()

function doDtiTractSub(imgs, matName, dtiDir, atlas)
dti = imgs.DTI;
if ~exist('atlas','var'), atlas = 'jhu'; end;
if isFieldSub(matName, ['dti_' atlas]), fprintf('Skipping tractography (DTI already computed) %s\n', imgs.ASL); return; end;
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
dti_u=prepostfixSub('', 'u', dti);
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
dti_faThr=prepostfixSub('', '_FA_thr', dti);
mask_dir=fullfile(pth, ['masks', atlasext]);
if ~exist(template_roiW,'file') || ~exist(dti_u,'file') || ~exist(dti_faThr,'file')
    fprintf('Can not find %s or %s or %s\n',template_roiW, dti_u, dti_faThr);
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
nPerm = 5000; %666
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
    %end %if voxels survive
end %for each region
if numel(commands) < 1
   fprintf('No regions survive thresholding with FA (poor normalization?) %s', dti);
   return;
end
fprintf ('computing probtrackx for %d regions (this may take a while)\n',numel(commands) );
doThreads(commands, prob_dir);
fprintf ('probtrackx2 took %f seconds to run.\n', toc(t_start) ); %t_start=tic;
nii_fiber_quantify(matName, dtiDir, atlas);
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
hdr = spm_vol(fnm);
img = spm_read_vols(hdr);
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
FA = prepostfixSub('', '_FA', imgs.DTI);
MD = prepostfixSub('', '_MD', imgs.DTI);
if ~exist(T1,'file') fprintf('Unable to find image: %s\n',T1); return; end; %required
if ~exist(FA,'file') fprintf('Unable to find image: %s\n',FA); return; end; %required
if ~exist(MD,'file') fprintf('Unable to find image: %s\n',MD); return; end; %required
FA = unGzSub (FA);
MD = unGzSub (MD);
nFA = rescaleSub(FA);
atlasImg = fullfile(fileparts(which('NiiStat')), 'roi' , [atlas '.nii']);
if ~exist(atlasImg,'file')
    error('Unable to find template %s', atlasImg);
end
[outhdr, outimg] = nii_reslice_target(atlasImg, '', T1, 0) ;
roiname = prepostfixSub('', atlasext, imgs.DTI);
if ~strcmpi(atlas,'jhu') &&  exist(roiname,'file')
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
%end unGzNameSub()

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

function fsldir = fslDirSub;
fsldir= '/usr/local/fsl/'; %CentOS intall location
if ~exist(fsldir,'dir')
    fsldir = '/usr/share/fsl/5.0/'; %Debian install location (from neuro debian)
    if ~exist(fsldir,'dir')
        error('FSL is not in the standard CentOS or Debian locations, please check you installation!');
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
[pth,nam] = filepartsSub(fnm);
bnm = fullfile(pth, [nam, '.bval']);
vnm = fullfile(pth, [nam, '.bvec']);
if ~exist(bnm, 'file') || ~exist(vnm, 'file')
    n = 0;
    return;
end
fileID = fopen(bnm,'r');
[A, n] = fscanf(fileID,'%g'); %#ok<ASGLU>
fclose(fileID);
%end bvalCountSub()

function imgs = doRestSub(imgs, matName)
if isempty(imgs.T1) || isempty(imgs.Rest), return; end; %we need these images
imgs.Rest = removeDotSub (imgs.Rest);
global ForceRest; %e.g. user can call "global ForceRest;  ForceRest = true;"
if isempty(ForceRest) && isFieldSub(matName, 'alf'), fprintf('Skipping Rest (already computed) %s\n', imgs.Rest); return; end;
delImgs('fdsw', imgs.Rest);
delImgs('fdswa', imgs.Rest);
nii_rest(imgs);
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
nii_nii2mat (prefixSub(['fdsw', prefix ],imgs.Rest), 'rest', matName)
nii_nii2mat (prefixSub(['palf_dsw', prefix ],imgs.Rest), 'alf', matName) %detrended 
nii_nii2mat (prefixSub(['palf_sw', prefix ],imgs.Rest), 'palf', matName) %conventional linear trends only
vox2mat(prefixSub(['wmean', prefix ],imgs.Rest), 'RestAve', matName);
vox2mat(prefixSub(['wbmean', prefix ],imgs.Rest), 'RestAve', matName);
%end doRestSub()

function delImgs(prefix, fnm) 
%e.g. delImgs('sw', 'X.nii') would delete wX.nii and swX.nii 
for i = 1: numel(prefix), 
    p = prefix(end-i+1:end); 
    d = prefixSub(p, fnm);
    if exist(d,'file'), delete(d); end; %delete last attempt
end;
%end delImgs() 

function idx = roiIndexSub(roiName)
[~, ~, idx] = nii_roi_list(roiName);
if idx < 1, error('Invalid roi name %s', roiName); end;
%end roiIndexSub()

function imgs = doAslSub(imgs, matName)
if isempty(imgs.T1) || isempty(imgs.ASL), return; end; %we need these images
imgs.ASL = removeDotSub (imgs.ASL);
if isFieldSub(matName, 'cbf'), fprintf('Skipping ASL (CBF already computed) %s\n', imgs.ASL); return; end;
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
if exist(prefixSub('wmeanCBF_0_src',asl),'file'), return; end; %already computed
[cbf, c1L, c1R, c2L, c2R] = nii_pasl12(asl, imgs.T1);
nii_nii2mat(cbf, 'cbf', matName); %2
stat = load(matName);
stat.cbf.nV = nV;
stat.cbf.c1L = c1L;
stat.cbf.c1R = c1R;
stat.cbf.c2L = c2L;
stat.cbf.c2R = c2R;
save(matName,'-struct', 'stat');
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
nii_enat_norm(imgs.T1,imgs.Lesion,imgs.T2);
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






















