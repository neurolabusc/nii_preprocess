function nii_bold_undistort(vols, meanbold, se, seRev)
%Align images so that origin and alignment roughly match MNI space
%  vols : cell string of BOLD image name(s) to undistort, assumes vols have been motion corrected
%  meanbold : average BOLD image (output of motion correct)
%  se : spin echo image with same parameters as vols
%  seRev : as se, but with reversed phase encoding
%Example
% nii_bold_undistort('fmri.nii', 'se.nii', 'seRev.nii');
% nii_bold_undistort({'fmri.nii.gz','fmri2.nii.gz'}, 'se.nii', 'seRev.nii');
% Taylor Hanayik 03/2016

if ~exist('vols','var') || isempty(vols) %no files specified
    vols = spm_select(inf,'image','Select distorted BOLD image(s)');
end
if ~exist('meanbold','var') || isempty(meanbold) %no files specified
    meanbold = spm_select(1,'image','Select mean BOLD image');
end
if ~exist('se','var') || isempty(se) %no files specified
    se = spm_select(inf,'image','Select Spin Echo image(s)');
end
if ~exist('seRev','var') || isempty(seRev) %no files specified
    seRev = spm_select(inf,'image','Select Reverse Phase Spin Echo image(s)');
end
if ischar(vols), vols = cellstr(vols); end %make cell if only one file passed
fprintf('%s version 03/2016\n', mfilename);
if isempty(which('spm')) || ~strcmp(spm('Ver'),'SPM12'), error('SPM12 required'); end;
if isempty(spm_figure('FindWin','Graphics')), spm fmri; end; %launch SPM if it is not running

fsldir = fslEnvSub; % set up fsl 
meanSE = realignSub(se); %make average image
meanSErev = realignSub(seRev); %make average image (phase reversed)

%Taylor - the coregEstimate does NOTHING for FSL - do you want to reslice?
%coregEstSub(meanbold,meanSE, prefix, meanSErev);
[meanSE, meanSErev] = coregEstWriteSub(meanbold, meanSE, meanSErev);
p.meanspinechoAP = meanSE;
p.meanspinechoPA = meanSErev;
p.boldname = vols;
p.spinechoAP = se;
p.spinechoPA = seRev;
%1.) merge spin echo images (AP then PA) into one file
% fslmerge -t bothSpinEchoFieldMap SpinEchoFieldMap_AP_TBV002.nii SpinEchoFieldMap_PA_TBV002.nii
p = fslMergeSpinEchoSub(fsldir, p);
%2.) run topup on combined SpinEcho image (~30min)
% topup --imain=bothSpinEchoFieldMap --datain=acqparams.txt --config=b02b0.cnf --out=topupResults --fout=myField --iout=ubothSpinEchoFieldMap --verbose
p = fslRunTopupSub(fsldir, p);
%5.) apply topup (linear interp may be faster, but spline (used here) is better)
%applytopup --imain=ztRest_AP_TBV002,Rest_PA_TBV002 --datain=acqparams.txt --inindex=1,4 --method=jac --topup=topupResults --out=uuRest_PA_TBV002
imgToFix = [p.boldname, meanbold];
fslApplyTopupSub(fsldir, imgToFix, p);

%========================================================================
%              Sub functions below
%========================================================================

function  [coregName, coregShadow] = coregEstWriteSub(ref, imgToCoreg, ShadowCoreg)
if nargin < 4
    otherImgs = [];
end
[p, n, x] = spm_fileparts(imgToCoreg);
coregName = fullfile(p, ['r', n, x]);
coregShadow = '';
%coregister fmri data to match T1 image
fprintf('Coregistering %s to match %s\n',imgToCoreg,ref);
%fMRIses = getsesvolsSubHier(prefix, fmriname);
matlabbatch{1}.spm.spatial.coreg.estwrite.ref = {ref};
matlabbatch{1}.spm.spatial.coreg.estwrite.source = {imgToCoreg};
if ~exist('ShadowCoreg','var') || isempty(ShadowCoreg) 
    matlabbatch{1}.spm.spatial.coreg.estwrite.other = {''};
else
    matlabbatch{1}.spm.spatial.coreg.estwrite.other = {ShadowCoreg};
    [p, n, x] = spm_fileparts(ShadowCoreg);
    coregShadow = fullfile(p, ['r', n, x]); 
end
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'nmi';
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.sep = [4 2];
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.fwhm = [7 7];
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.interp = 4;
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.mask = 0;
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.prefix = 'r';

spm_jobman('run',matlabbatch);
%end coregEstSub()



function fsldir = fslEnvSub
fsldir = fslDirSub;
curpath = getenv('PATH');
k = strfind(curpath, fsldir);
if isempty(k)
	setenv('PATH',sprintf('%s:%s',fullfile(fsldir,'bin'),curpath));
end
%end fslEnvSub()

function fsldir = fslDirSub
fsldir= '/usr/local/fsl/'; %CentOS intall location
if ~exist(fsldir,'dir')
    fsldir = '/usr/share/fsl/5.0/'; %Debian install location (from neuro debian)
    if ~exist(fsldir,'dir')
        error('FSL is not in the standard CentOS or Debian locations, please check you installation!');
    end
end
%end fslDirSub()

% function coregEstSub(ref, imgToCoreg, prefix, otherImgs)
% if nargin < 4
%     otherImgs = [];
% end
% %coregister fmri data to match T1 image
% fprintf('Coregistering %s to match %s\n',imgToCoreg,ref);
% %fMRIses = getsesvolsSubHier(prefix, fmriname);
% matlabbatch{1}.spm.spatial.coreg.estimate.ref = {ref};
% matlabbatch{1}.spm.spatial.coreg.estimate.source = {imgToCoreg};
% if ~exist('otherImgs','var') || isempty(otherImgs) 
%     matlabbatch{1}.spm.spatial.coreg.estimate.other = {''};
% else
%     otherImgsSes = getsesvolsSubFlat(prefix, otherImgs);
%     matlabbatch{1}.spm.spatial.coreg.estimate.other = otherImgsSes;
% end;
% matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
% matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
% matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
% matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];
% spm_jobman('run',matlabbatch);
% %end coregEstSub()


function meanImgOut = realignSub(inName)
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
    matlabbatch{1}.spm.spatial.realign.estwrite.roptions.which = [0 1];
    matlabbatch{1}.spm.spatial.realign.estwrite.roptions.interp = 4;
    matlabbatch{1}.spm.spatial.realign.estwrite.roptions.wrap = [0 0 0];
    matlabbatch{1}.spm.spatial.realign.estwrite.roptions.mask = 1;
    matlabbatch{1}.spm.spatial.realign.estwrite.roptions.prefix = 'r';
    spm_jobman('run',matlabbatch);
else
    disp(['File already exists: ' meanImgOut]);
end
%end  realignSub()


function p = fslMergeSpinEchoSub(fsldir,p)
[pth,~,ext] = spm_fileparts(p.spinechoAP); % get path of spinecho images
p.spinechoMerged = fullfile(pth,['mSpinEcho' ext]); % construct merged spinecho filename. Order in new file will be AP then PA
cmdstr = ['fslmerge -t ' p.spinechoMerged ' ' p.meanspinechoAP ' ' p.meanspinechoPA];
disp(['Running: ' cmdstr]);
if ~exist(p.spinechoMerged,'file') % only run if this file doesn't exist already
    doFslCmd(cmdstr);
else
    disp(['File already exists: ' p.spinechoMerged])
end
p = fslWriteTopUpParamFile(fsldir,p);
%end fslMergeSpinEchoSub


function p = fslWriteTopUpParamFile(fsldir,p)
if ~exist(p.spinechoMerged,'file'); error(['Cant find file: ' p.spinechoMerged]); end;
[pth,nm,~] = spm_fileparts(p.spinechoMerged);
p.acqparamsfile = fullfile(pth,['acqparams_' nm '.txt']);
fid = fopen(p.acqparamsfile,'w'); %open file for writing, discard existing contents if any
%read number of volumes in merged spinecho file
cmdstr = ['fslnvols ' p.spinechoMerged];
disp(cmdstr);
[~,r] = doFslCmd( cmdstr);
nvols = str2double(r);
p.inindex1 = 1;
p.inindex2 = (nvols/2)+1;
disp([nm ' has ' num2str(nvols) ' volumes'])
disp('Writing params file for topup');
for i = 1:nvols
    if i <= nvols/2
        % write AP params
        fprintf(fid,'0 -1 0 0.062\n');
    else
        %write PA params
        fprintf(fid,'0 1 0 0.062\n');
    end
end
fclose(fid);
%end fslWriteTopUpParamFile


function p = fslRunTopupSub(fsldir, p)
prevpth = pwd;
%topup --imain=bothSpinEchoFieldMap --datain=acqparams.txt --config=b02b0.cnf --out=topupResults --fout=spinechoField --iout=umSpinEchoFieldMap --verbose
[pth,nm,~] = spm_fileparts(p.spinechoMerged);
if ~isempty(pth); cd(pth); end;
p.topupoutfile = 'topupResults';
p.fieldfile = 'computedField';
p.umspinechoMerged = ['u' nm];
cmdstr = ['topup --imain=' p.spinechoMerged ' --datain=' p.acqparamsfile ' --config=b02b0.cnf --out=' p.topupoutfile ' --fout=' p.fieldfile ' --iout=' p.umspinechoMerged ' --verbose'];
disp(cmdstr);
if ~exist(fullfile(pth,[p.umspinechoMerged '.nii']),'file') % only run topup if it hasn't been run before since it take a long time to finish
    doFslCmd(cmdstr);
else
    disp('It looks like topup has been run on this data previously');
end
if ~isempty(prevpth); cd(prevpth); end;
%fslRunTopupSub


function [p, outname] = fslApplyTopupSub(fsldir, imgsToFix, p)
%apply topup (linear interp may be faster, but spline (used here) is better)
%applytopup --imain=ztRest_AP_TBV002,Rest_PA_TBV002 --datain=acqparams.txt --inindex=1,4 --method=jac --topup=topupResults --out=uuRest_PA_TBV002
%setFSLenv;
numFiles = size(imgsToFix,2);
for i = 1:numFiles
    img = imgsToFix{i};
    prevpth = pwd;
    [pth, nm, ext] = spm_fileparts(img);
    img = fullfile(pth,[nm ext]);
    outname = fullfile(pth,['u' nm ext]);
    if ~isempty(pth), cd(pth); end;
    %cmdstr = (['applytopup --imain=' APfilename ',' PAfilename ' --datain=' p.acqparamsfile ' --inindex=' num2str(p.inindex1) ',' num2str(p.inindex2) ' --method=jac ' '--topup=' p.topupoutfile ' --out=' outname]);
    %cmdstr = (['applytopup --imain=' img ' --datain=' p.acqparamsfile ' --inindex=' num2str(p.inindex2) ' --method=jac ' '--topup=' fullfile(pth,p.topupoutfile) ' --out=' outname]);
    % inindex=1 : bold images match SE, not SErev
    cmdstr = (['applytopup --imain=' img ' --datain=' p.acqparamsfile ' --inindex=1 --method=jac ' '--topup=' fullfile(pth,p.topupoutfile) ' --out=' outname]);
    
    if ~exist(outname,'file');
        disp(cmdstr);
        doFslCmd(cmdstr);
    else
        disp(['File already exists: ' outname]);
    end
    if ~isempty(prevpth), cd(prevpth); end;
end
%end fslApplyTopupSub

function [sesvols] = getsesvolsSub(sesvol1)
% input: single volume from 4D volume, output: volume list
%  example 4D file with 3 volumes input= 'img.nii', output= {'img.nii,1';'img.nii,2';'img.nii,3'}
[pth,nam,ext,vol] = spm_fileparts( deblank (sesvol1));
sesname = fullfile(pth,[nam, ext]);
hdr = spm_vol(sesname);
nvol = length(hdr);
sesvols = cellstr([sesname,',1']);
if nvol < 2, return; end;
for vol = 2 : nvol
    sesvols = [sesvols; [sesname,',',int2str(vol)]  ]; %#ok<AGROW>
end;
%end getsesvolsSub()

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

function command  = fslCmdSub (command)
fsldir = fslDirSub;
cmd=sprintf('sh -c ". %setc/fslconf/fsl.sh; export FSLOUTPUTTYPE=NIFTI; ',fsldir);
command = [cmd command '"'];
[status,output] = system(command);
%fslCmdSub
