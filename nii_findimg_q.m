function nii_findimg_q (baseDir)
%find relevant images for patients
% baseDir: root folder that holds images
%Example Structure (%s is basefolder, with two subjects)
% ~\p001\T2_P001.nii
% ~\p001\LS_P001.nii
% ~\p001\T1_P001.nii
% ~\p001\fMRI_P001.nii
% ~\p002\T2_P002.nii
% ~\p002\LS_P002.nii
% ~\p002\T1_P002.nii
%Example
% nii_findimg_q %use gui
% nii_findimg_q(pwd)
% nii_findimg_q('/Users/rorden/Desktop/pre')

if ~exist('baseDir','var') || isempty(baseDir)
    baseDir = pwd; %uigetdir('','Pick folder that contains all subjects');
end
subjDirs = subFolderSub(baseDir);
subjDirs = sort(subjDirs);
subjDirs = {'P001'};
% 'LM1052', 'LM1054', 'LM1056RHD', 'P001', 'P011', 'P012', 'P021', 'P023', 'P024', 'P026', 'P027', 'P098', 'P106', 'P134', 'P171'

%subjDirs = {'LM1052', 'LM1054', 'LM1056RHD', 'P001', 'P011', 'P012', 'P021', 'P023', 'P024', 'P026', 'P027', 'P098', 'P106', 'P134', 'P171'
%subjDirs = {'LM1052', 'LM1054', 'LM1056RHD'}
fprintf('Found %d folders (subjects) in %s\n',numel(subjDirs), baseDir);
for s = 1:numel(subjDirs)% :-1: 1 %for each participant
    subjDir = [baseDir,filesep, deblank(subjDirs{s}) ]; %no filesep
    imgs = subjStructSub(deblank(subjDirs{s}));
    imgs.Lesion = imgfindSub(imgs.Lesion,strvcat('LS_'), subjDir); %#ok<REMFF1>
    imgs.T1 = imgfindSub(imgs.T1,'T1_', subjDir);
    imgs.T2 = imgfindSub(imgs.T2,'T2_', subjDir);
    imgs.ASL = imgfindSub(imgs.ASL,'ASL_', subjDir);
    imgs.DTI = imgfindSub(imgs.DTI,'APDTI_', subjDir);
    if isempty(imgs.DTI) %cannot find 'APDTI_' look for standard DTI
        imgs.DTI = imgfindSub(imgs.DTI,'DTI_', subjDir);
    else
        imgs.DTIrev = imgfindSub(imgs.DTIrev,'PADTI_', subjDir);
    end
    %nii_preprocess(imgs);
    imgs = unGzAllSub(imgs); %all except DTI - fsl is OK with nii.gz
    if ~isempty(imgs.T1) && ~isempty(imgs.Lesion)
        matName = [subjDirs{s} '.mat'];
        doT1Sub(imgs, matName); %normalize T1
        %doI3MSub(imgs, matName);%next i3m
        %doAslSub(imgs, matName);
        doDtiSub(imgs, matName);
        %doFaMdSub(imgs, matName);
        error('x');
    end
    checkDims(imgs);
end
%end nii_findimg()

function doFaMdSub(imgs, matName)
if isempty(imgs.T1) || isempty(imgs.DTI), return; end; %required
T1 = prefixSub('wb',imgs.T1); %warped brain extracted image
FA = prepostfixSub('', '_FA', imgs.DTI);
MD = prepostfixSub('', '_MD', imgs.DTI);
if ~exist(T1,'file') || ~exist(FA,'file') || ~exist(MD,'file') , return; end; %required
FA = unGzSub (FA);
MD = unGzSub (MD);
nFA = rescaleSub(FA);
nii_setOrigin12({nFA, FA, MD}, 1,false); %rFA looks like T1
oldNormSub( {nFA, FA,MD}, T1);
FA = prepostfixSub('w', '_FA', imgs.DTI);
MD = prepostfixSub('w', '_MD', imgs.DTI);
nii_nii2mat(FA, 6, matName);
nii_nii2mat(MD, 8, matName);
%end doFaMdSub()

function oldNormSub(src, tar)
%coregister T2 to match T1 image, apply to lesion
if isempty(src) || isempty(tar), return; end;
matlabbatch{1}.spm.tools.oldnorm.estwrite.subj.source = src(1);
matlabbatch{1}.spm.tools.oldnorm.estwrite.subj.wtsrc = '';
matlabbatch{1}.spm.tools.oldnorm.estwrite.subj.resample = src;
matlabbatch{1}.spm.tools.oldnorm.estwrite.eoptions.template = {tar};
matlabbatch{1}.spm.tools.oldnorm.estwrite.eoptions.weight = '';
matlabbatch{1}.spm.tools.oldnorm.estwrite.eoptions.smosrc = 8;
matlabbatch{1}.spm.tools.oldnorm.estwrite.eoptions.smoref = 0;
matlabbatch{1}.spm.tools.oldnorm.estwrite.eoptions.regtype = 'mni';
matlabbatch{1}.spm.tools.oldnorm.estwrite.eoptions.cutoff = 25;
matlabbatch{1}.spm.tools.oldnorm.estwrite.eoptions.nits = 16;
matlabbatch{1}.spm.tools.oldnorm.estwrite.eoptions.reg = 1;
matlabbatch{1}.spm.tools.oldnorm.estwrite.roptions.preserve = 0;
matlabbatch{1}.spm.tools.oldnorm.estwrite.roptions.bb = [-78 -112 -70; 78 76 85]; ;
matlabbatch{1}.spm.tools.oldnorm.estwrite.roptions.vox = [1 1 1];
matlabbatch{1}.spm.tools.oldnorm.estwrite.roptions.interp = 1;
matlabbatch{1}.spm.tools.oldnorm.estwrite.roptions.wrap = [0 0 0];
matlabbatch{1}.spm.tools.oldnorm.estwrite.roptions.prefix = 'w';
spm_jobman('run',matlabbatch);
%end oldNormSub()

function nam = prepostfixSub (pre, post, nam)
[p, n, x] = filepartsSub(nam);
nam = fullfile(p, [pre, n, post, x]);
if ~exist(nam) && strcmpi(x,'.nii')
   nam = fullfile(p, [pre, n, post, '.nii.gz']);
end
if ~exist(nam) && strcmpi(x,'.nii.gz')
   nam = fullfile(p, [pre, n, post, '.nii']);
end
%end prefixSub()

function [p,n,x] = filepartsSub (fnm)
if isempty(fnm), return; end;
fnm = char(deblank(fnm));
[p,n,x] = fileparts(fnm);
if strcmpi(x,'.gz') %.nii.gz
    [p2,n2,x2] = spm_fileparts(n);
    n = n2;
    x = [x2, x];
end;
%end filepartsSub()

function fname = rescaleSub(fname)
hdr = spm_vol(fname);
img = spm_read_vols(hdr);
if min(img(:)) < 0, error('Image intensity must be positive'); end;
if min(img(:)) == max(img(:)), error('Image has no variability'); end;
img = img - min(img(:)); %scale from zero
img = img/max(img(:)); %scale from 0..1
[pth nm ext] = spm_fileparts(fname);
img = power(img, 0.5);
fname = fullfile(pth, ['n' nm ext]);
hdr.fname = fname;
spm_write_vol(hdr,img);
%end rescaleSub()

%end


function checkDims(imgs)
if isempty(imgs.T1) || isempty(imgs.Lesion), return; end; %we need these images
hdr = spm_vol(imgs.Lesion);
[pth nm] = spm_fileparts(imgs.Lesion);
fprintf('%s\t%d\t%d\t%d\n', nm, hdr.dim(1), hdr.dim(2), hdr.dim(3) );

%end checkDims()


function doDtiSub(imgs, matName)
if isempty(imgs.T1) || isempty(imgs.DTI), return; end; %required
betT1 = prefixSub('b',imgs.T1); %brain extracted image
eT1 = prefixSub('e',imgs.T1); %enantimorphic image
if ~exist(betT1,'file') || ~exist(eT1,'file'), return; end; %required
n = bvalCountSub(imgs.DTI);
if (n < 1)
    fprintf('UNABLE TO FIND BVECS/BVALS FOR %s\n', imgs.DTI);
    return
end
if (n < 12)
    fprintf('INSUFFICIENT BVECS/BVALS FOR %s\n', imgs.DTI);
    return
end
%1 - eddy current correct
command= [fileparts(which(mfilename)) filesep 'dti_1_eddy.sh'];
if isempty(imgs.DTIrev)
    command=sprintf('%s "%s"',command, imgs.DTI);
else
    nr = bvalCountSub(imgs.DTIrev);
    if (nr ~= n)
        fprintf('BVECS/BVALS DO NOT MATCH %s %s\n', imgs.DTI, imgs.DTIrev);
        return
    end
    command=sprintf('%s "%s" "%s"',command, imgs.DTI, imgs.DTIrev);
end
%rmBedpostDir(imgs.DTI);
doFslCmd (command);
return;
%2 warp template
dti = deblank(imgs.DTI);
command= [fileparts(which(mfilename)) filesep 'dti_2_warp_template.sh'];
if isempty(eT1)
    command=sprintf('%s "%s" "%s"',command, dti, betT1);
else
    command=sprintf('%s "%s" "%s" "%s"',command, dti, betT1, eT1);
end
doFslCmd (command);
%3 tractography
command= [fileparts(which(mfilename)) filesep 'dti_3_tract.sh'];
command=sprintf('%s "%s" ',command, dti);
doFslCmd (command);
%end doDtiSub()

function rmBedpostDir(imgname);
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

function doFslCmd (command)
fsldir= '/usr/local/fsl/';
setenv('FSLDIR', fsldir);
curpath = getenv('PATH');
setenv('PATH',sprintf('%s:%s',fullfile(fsldir,'bin'),curpath));
cmd=sprintf('sh -c ". %setc/fslconf/fsl.sh; ',fsldir);
cmd = [cmd command '"'];
fprintf('Running \n %s\n', cmd);
system(cmd);
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
[pth,nam] = nii_filepartsSub(fnm);
bnm = fullfile(pth, [nam, '.bval']);
vnm = fullfile(pth, [nam, '.bvec']);
if ~exist(bnm, 'file') || ~exist(vnm, 'file')
    n = 0;
    return;
end
fileID = fopen(bnm,'r');
[A, n] = fscanf(fileID,'%g'); %#ok<ASGLU>
fclose(fileID);

%end


function doAslSub(imgs, matName)
if isempty(imgs.T1) || isempty(imgs.ASL), return; end; %we need these images
nV = nVolSub (imgs.ASL) ;
[mx, ind] = max(nV);
if mx < 73, fprintf('not enough ASL volumes for %s\n', matName); end;
asl = imgs.ASL(ind,:);
if ~exist(prefixSub('b',imgs.T1),'file'), return; end; %required
if exist(prefixSub('wmeanCBF_0_src',asl),'file'), return; end; %already computed
[cbf, c1L, c1R, c2L, c2R] = nii_pasl12(asl, imgs.T1);
nii_nii2mat(cbf, 2, matName);
stat = load(matName);
stat.cbf.nV = nV;
stat.cbf.c1L = c1L;
stat.cbf.c1R = c1R;
stat.cbf.c2L = c2L;
stat.cbf.c2R = c2R;
save(matName,'-struct', 'stat');
if (c1R < c2R)
    fid = fopen('errors.txt','a');
    fprintf(fid, 'ASL CBF higher in white matter\t%s\n', matName);
    fclose(fid);
end
%end doPaslSub()

function nVol = nVolSub (fnm)
%Report number of volumes
% v= nVols('img.nii');
% v= nVolS(strvcat('ASL_P005.nii', 'ASL_P005_1.nii'))
% v= nVolS({'ASL_P005.nii', 'ASL_P005_1.nii'})
nVol = [];
fnm = cellstr(fnm);
for v = 1 : numel(fnm)
    hdr = spm_vol(deblank(char(fnm(v,:))));
    nVol = [nVol numel(hdr)]; %#ok<AGROW>
end
%end nVol()

function doI3MSub(imgs, matName)
if isempty(imgs.T1), return; end; %we need these images
w = prefixSub('w',imgs.T1); %warped T1
if  ~exist(w, 'file')  return; end; %exit: we require warped T1
i3m = prefixSub('zw',imgs.T1);
if exist(i3m,'file'), return; end; %i3m already computed
nii_i3m(w,'',0.5,10,0.25,1); %i3m T1 image
nii_nii2mat(i3m, 4, matName);
%end doI3MSub()

function doT1Sub(imgs, matName)
if isempty(imgs.T1) || isempty(imgs.Lesion), return; end; %we need these images
if size(imgs.T1,1) > 1 || size(imgs.T2,1) > 1 || size(imgs.Lesion,1) > 1
    error('Require no more than one image for these modalities: T1, T2, lesion');
end;
b = prefixSub('b',imgs.T1);
if exist(b,'file'), return; end; %this stage was already run
imgs.T1
imgs.T2
imgs.Lesion
class(imgs.T1)
nii_enat_norm(imgs.T1,imgs.Lesion,imgs.T2);
nii_nii2mat(prefixSub('wr',imgs.Lesion), 1, matName);
%end doT1Sub()

function nam = prefixSub (pre, nam)
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

function [fnm, isGz] = unGzSub (fnm)
if isempty(fnm), return; end;
innames = fnm;
fnm = [];
for i = 1: size(innames,1)
    f = deblank(innames(i,:));
    f = unGzCSub(f);
    fnm = strvcat(fnm, f);
end
%end unGzSub()

function [fnm, isGz] = unGzCSub (fnm)
if isempty(fnm), return; end;
[pth,nam,ext] = spm_fileparts(fnm);
isGz = false;
if strcmpi(ext,'.gz') %.nii.gz
    ofnm = fnm;
    fnm = char(gunzip(fnm));
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
if size(imgName) > 1, imgName = imgName(1,:); end;
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

function [pth,nam,ext,num] = nii_filepartsSub(fname)
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
%end nii_filepartsSub()