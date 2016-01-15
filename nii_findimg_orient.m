function nii_findimg_orient (baseDir)
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
fprintf('Found %d folders (subjects) in %s\n',size(subjDirs,1), baseDir);
sv = {'P95','P96','P97','P99','P100','P106','LM104','LM1011', 'LM1019', 'LM1021', 'LM1022', 'LM1023', 'LM1029', 'LM1034', 'LM1035', 'LM1036', 'LM1037', 'LM1039', 'LM1047', 'LM1050', 'LM1052', 'LM1053', 'P001', 'P002', 'P004', 'P005', 'P006', 'P007', 'P009', 'P011', 'P012', 'P013', 'P014', 'P016', 'P017', 'P018', 'P019', 'P020', 'P021', 'P022', 'P023', 'P024', 'P025', 'P026', 'P027', 'P029', 'P030', 'P031', 'P036', 'P040', 'P042', 'P046', 'P059', 'P061', 'P069', 'P072', 'P074', 'P075', 'P076', 'P079', 'P083', 'P084', 'P086', 'P089', 'P092', 'P094', 'P095', 'P095', 'P125', 'P134', 'P137', 'P138', 'P141', 'P148', 'P150', 'P171', 'P175', 'P177', 'P180', 'P187', 'P188', 'P189', 'P192', 'P193'};

for s = 1 : size(subjDirs,1)% :-1: 1 %for each participant

    
    IndexC = strfind(sv,deblank(subjDirs{s}) ) ;
    
    Index = find(not(cellfun('isempty', IndexC)));
    if ~isempty(Index) continue; end
    
    subjDir = [baseDir,filesep, deblank(subjDirs{s}) ]; %no filesep
    imgs = subjStructSub( {deblank(subjDirs{s})});
    imgs.Lesion = imgfindSub(imgs.Lesion,strvcat('LS_'), subjDir); %#ok<REMFF1>
    imgs.T1 = imgfindSub(imgs.T1,'T1_', subjDir);
    imgs.fMRI = imgfindSub(imgs.fMRI,'fMRI_', subjDir);
    imgs.T2 = imgfindSub(imgs.T2,'T2_', subjDir); 
    imgs.ASL = imgfindSub(imgs.ASL,'ASL_', subjDir);
    imgs.DTI = imgfindSub(imgs.DTI,'APDTI_', subjDir);
    if isempty(imgs.DTI) %cannot find 'APDTI_' look for standard DTI
        imgs.DTI = imgfindSub(imgs.DTI,'DTI_', subjDir);
    else
        imgs.DTIrev = imgfindSub(imgs.DTIrev,'PADTI_', subjDir);
    end
    if isempty(imgs.T1), continue; end;
    %fprintf('%d %s -> %s\n', s, subjDir, imgs.T1);
    if isempty(imgs.Lesion), continue; end;
    newPath = fullfile('/Volumes/Flash/Pre',deblank(subjDirs{s}));
    if ~isempty(imgs.DTIrev)
        
        imgs.DTIrev = cpImgSub(newPath,imgs.DTIrev);
        imgs.DTI = cpImgSub(newPath,imgs.DTI);
        nii_setOrigin12({imgs.DTI, imgs.DTIrev}, 3, true); %EPI
        fprintf('-->\t%s\t%s\n',imgs.DTIrev, imgs.DTI);
    elseif ~isempty(imgs.DTI)
        imgs.DTI = cpImgSub(newPath,imgs.DTI);
        nii_setOrigin12({imgs.DTI}, 3, true); %EPI
        fprintf('-->\t%s\n',imgs.DTI);
            
    end
    
%     if ~strcmpi(deblank(subjDirs{s}), 'P081'), continue; end;
% 
%     if isempty(imgs.T2)
%         checkDimsSub(imgs.T1, imgs.Lesion);
%         nii_setOrigin12({imgs.T1, imgs.Lesion}, 1,true); %T1
%     else
%        nii_setOrigin12(imgs.T1, 1, true); %T1
%        checkDimsSub(imgs.T2, imgs.Lesion);
%        nii_setOrigin12({imgs.T2, imgs.Lesion}, 2, true); %T2
%     end
%     if ~isempty(imgs.ASL)
%         nii_setOrigin12(imgs.ASL, 3, true); %EPI
%     end
%     if ~isempty(imgs.fMRI)
%         nii_setOrigin12(imgs.fMRI, 3, true); %EPI
%     end 
%     if ~isempty(imgs.DTIrev)
%         nii_setOrigin12({imgs.DTI, imgs.DTIrev}, 3, true); %EPI
%     elseif ~isempty(imgs.DTI)
%         nii_setOrigin12(imgs.DTI, 3, true); %EPI 
%     end 
    
  %end
  
end
%end nii_findimg()

function checkDimsSub(anat, lesion)
if ~exist(anat,'file'),return; end;
if ~exist(lesion,'file'), return; end;
tarhdr = spm_vol(anat);
hdr = spm_vol(lesion);
img = spm_read_vols(hdr);
if prod(hdr.dim) ~= prod(hdr.dim)
    error('%s really does not match %s',anat, lesion);
end
if (hdr.dim(1) == tarhdr.dim(1)) && (hdr.dim(2) == tarhdr.dim(2)) && (hdr.dim(3) == tarhdr.dim(3))
    return; %need to reslice
end
[outhdr, outimg] = nii_reslice_target(hdr, img, tarhdr, 1) ;
delete(lesion);
outhdr.fname = hdr.fname;
spm_write_vol(outhdr,outimg);
%end checkDimsSub()

function newName = cpImgSub(newPath,oldName)
if size(oldName) > 1, oldName = oldName(1,:); end;
%oldName = deblank(oldName);
[newName, isGz] = stripExtSub(oldName);
%fprintf('%s\n', newName);
if isempty(newPath) || isempty(oldName)
   return;
end
if ~exist(newPath, 'dir')
    mkdir(newPath);
end

[oldPath,nam] = fileparts(newName);
if (isGz)
    newName = doCpSub(oldPath,newPath,[nam '.nii.gz']);
else
    newName = doCpSub(oldPath,newPath,[nam '.nii']);
end
%newName = fullfile(newPath,nam);
%doCpSub(oldPath,newPath,['m' nam ext]);


doCpSub(oldPath,newPath,[nam '.bvec']);
doCpSub(oldPath,newPath,[nam '.bval']);
%end cpImgSub()

function [fnm, isGz] = stripExtSub (fnm)
[p,n,x] = fileparts(fnm);
isGz = false;
if strcmpi(x,'.gz') %.nii.gz
    isGz = true;
    [p,n,x] = fileparts(fullfile(p,n));
end;  
fnm = fullfile(p,n);
%end unGzSub()

function newName = doCpSub(oldPath,newPath,namext)
oldName = fullfile(oldPath,namext);
if exist(oldName, 'file') ~= 2 
    return;
end
newName = fullfile(newPath,namext);
copyfile(oldName,newName);
%end doCpSub

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
doFslCmd (command);
%2 warp template
command= [fileparts(which(mfilename)) filesep 'dti_2_warp_template.sh'];
if isempty(et1)
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

function [pth,nam,ext,num] = nii_filepartsSub(fname)
% extends John Ashburner's spm_fileparts.m to include '.nii.gz' as ext
num = '';
if ~ispc, fname = strrep(fname,'\',filesep); end
[pth,nam,ext] = fileparts(fname);
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