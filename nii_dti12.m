function nii_dti12 (dti, dtiRev, lesion, t1, newPath)
%process DTI data
% dti : name of DTI image
% dtiRev : (optional) DTI image with reverse phase encoding
% lesion : name of lesion image
% t1 : name of T1-weighted image 
% newPath : optional - images processed in temporary folder...
%Examples 
% nii_dti12('1DTI_P084_1.nii.gz','', 'LS_P084.nii', 'T1_P084')
% nii_dti12('1DTI_P084_1.nii.gz','', 'LS_P084', 'T1_P084', [pwd filesep 'temp'])
% nii_dti12 %user selects files manually


if ~exist('dti','var') dti = spm_select(1,'image','Select DTI image'); end;
if ~exist('dtiRev','var') dtiRev = spm_select(1,'image','(Optional) Select reversed DTI'); end;
if ~exist('lesion','var') lesion = spm_select(1,'image','Select lesion map'); end;
if ~exist('t1','var') t1 = spm_select(1,'image','(Optional) Select T1 image'); end;
if ~exist('newPath','var') newPath = ''; end;
dti = stripExtSub (dti);
dtiRev = stripExtSub (dtiRev);
lesion = stripExtSub (lesion);
t1 = stripExtSub(t1);
%next: strip extension, move folder if requested
dti = cpImgSub (newPath,dti);
dtiRev = cpImgSub (newPath,dtiRev);
lesion = cpImgSub (newPath,lesion);
t1 = cpImgSub(newPath,t1);    
%next: process data
tic
do_dti_core(dti, dtiRev)
fprintf('DTI processing took %g seconds',toc);
%end nix()

function newName = cpImgSub(newPath,oldName)
newName = stripExtSub(oldName);
if isempty(newPath) || isempty(oldName)
   return;
end
if ~exist(newPath, 'dir')
    mkdir(newPath);
end
newName = [newName '.nii'];
[oldPath,nam] = fileparts(newName);
newName = fullfile(newPath,nam);
%doCpSub(oldPath,newPath,['m' nam ext]);
doCpSub(oldPath,newPath,[nam '.nii']);
doCpSub(oldPath,newPath,[nam '.nii.gz']);
doCpSub(oldPath,newPath,[nam '.bvec']);
doCpSub(oldPath,newPath,[nam '.bval']);
%end cpImgSub()

function doCpSub(oldPath,newPath,namext);
oldName = fullfile(oldPath,namext);
if exist(oldName, 'file') ~= 2 
    return;
end
newName = fullfile(newPath,namext);
copyfile(oldName,newName);
%end doCpSub

function do_dti_core(dti, dtiRev)
pth = fileparts(which(mfilename));
if isempty(pth), pth = pwd; end;
script = fullfile(pth, 'dti.sh');
if ~exist(script,'file'), error('%s: script (%s) not found',mfilename,script); end
fsldir= '/usr/local/fsl/';
if ~exist(fsldir,'dir'), error('%s: fsldir (%s) not found',mfilename,fsldir); end
fsldir= '/usr/local/fsl/';
if ~exist(fsldir,'dir'), error('%s: fsldir (%s) not found',mfilename,fsldir); end
command=sprintf('sh -c ". ${FSLDIR}/etc/fslconf/fsl.sh; %s "%s" "%s""\n', script, dti, dtiRev);
system(command);
%end do_dti()

% function do_dti(dti, dtiRev, lesion, t1)
% pth = fileparts(mfilename);
% if isempty(pth), pth = pwd; end;
% script = fullfile(pth, 'dti.sh');
% if ~exist(script,'file'), error('%s: script (%s) not found',mfilename,script); end
% fsldir= '/usr/local/fsl/';
% if ~exist(fsldir,'dir'), error('%s: fsldir (%s) not found',mfilename,fsldir); end
% fsldir= '/usr/local/fsl/';
% if ~exist(fsldir,'dir'), error('%s: fsldir (%s) not found',mfilename,fsldir); end
% command=sprintf('sh -c ". ${FSLDIR}/etc/fslconf/fsl.sh; %s "%s" "%s" "%s" "%s""\n', script, dti, dtiRev, lesion, t1);
% system(command);
%end do_dti()

function fnm = stripExtSub (fnm)
[p,n,x] = fileparts(fnm);
if strcmpi(x,'.gz') %.nii.gz
    [p,n,x] = fileparts(fullfile(p,n));
end;  
fnm = fullfile(p,n);
%end unGzSub()
