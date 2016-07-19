function nii_findimg (baseDir)
%find relevant images for patients
% baseDir: root folder that holds images
%Example Structure (%s is basefolder, with two subjects)
% ~\p1\t1\patient1T1.nii
% ~\p1\lesion\patient1Lesion.nii
% ~\MUSC2\t1FOLDER\garbage.nii
% ~\MUSC2\lesionFOLD\garbage.nii
% ~\MUSC2\DTItrash\garbage.nii
%Example
% nii_findimg %use gui
% nii_findimg(pwd)
% nii_findimg('/Users/rorden/Desktop/merc')

if ~exist('baseDir','var') || isempty(baseDir)
    baseDir = uigetdir('','Pick folder that contains all subjects');
end
diary ON
subjDirs = subFolderSub(baseDir);
fprintf('Found %d folders (subjects) in %s\n',size(subjDirs,1), baseDir);
for s = 1:  size(subjDirs,1) %for each participant
    subjDir = [baseDir,filesep, deblank(subjDirs{s}) ]; %no filesep
    modalDirs = subFolderSub(subjDir);
    if size(modalDirs,1) < 1
        fprintf(' WARNING: No folders for %s\n',subjDir);
    else
        imgs = subjStructSub(deblank(subjDirs{s}));
        fprintf(' Found %d folders (modalities) in %s\n',size(modalDirs,1), subjDir);
        for m = 1: size(modalDirs,1) %for each participant
            modalDir = [subjDir,filesep, deblank(modalDirs{m}) ]; %no filesep
            imgs.fMRI = imgfindSub(imgs.fMRI,'fMRI', modalDir);
            imgs.ASL = imgfindSub(imgs.ASL,'ASL', modalDir);
            imgs.DTI = imgfindSub(imgs.DTI,strvcat('DTI','27AVE','DIFF'), modalDir);  %#ok<REMFF1>
             imgs.Lesion = imgfindSub(imgs.Lesion,strvcat('LESION','Lesion','LS'), modalDir); %#ok<REMFF1>
            imgs.T1 = imgfindSub(imgs.T1,'T1', modalDir);
            imgs.T2 = imgfindSub(imgs.T2,'T2', modalDir);
        end
        imgs.name = removeMusc(imgs.name);
        nii_preprocess(imgs);
    end
end
diary OFF
%end nii_findimg()

function str = removeMusc(str)
pos = strfind(lower(char(str)),lower('MUSC_'));
if isempty(pos), return; end
str = strrep(str, 'MUSC_', '');
%end removeMusc

function imgs = subjStructSub(subjName)
imgs.name = subjName;
imgs.ASL = '';
imgs.DTI = '';
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
if ~isStringInKey (nam, imgKey), return; end;
%if isempty(strfind(lower(char(nam)), lower(imgKey))), return; end;
if exist([inDir,filesep, 'Native'],'file')
     imgName = imgfindXSub (imgName, imgKey, [inDir,filesep, 'Native']);
end
if exist([inDir,filesep, 'RUN1'],'file')
    imgName = imgfindXSub (imgName, imgKey, [inDir,filesep, 'RUN1']);
end
if exist([inDir,filesep, 'RUN2'],'file')
    imgName = imgfindXSub (imgName, imgKey, [inDir,filesep, 'RUN1']);
end
imgName = imgfindXSub (imgName, imgKey, inDir);
if isempty(imgName)
    fprintf('WARNING: unable to find any "%s" images in folder %s\n',deblank(imgKey(1,:)), inDir);
end
%end imgfindSub()

function imgName = imgfindXSub (imgName, imgKey, inDir)
nameFiles = subFileSub(inDir);
nameFiles = sort(nameFiles); %take first file for multiimage sequences, e.g. ASL
for i=1:size(nameFiles,1)
    if isStringInKey (char(nameFiles(i)), imgKey) && isImgSub(char(nameFiles(i)))
        imgName = strvcat(imgName, [inDir, filesep, char(nameFiles(i))]);
        return
    end; %do not worry about bvec/bval
end
%end imgfindXSub()

function isKey = isStringInKey (str, imgKey)
isKey = false;
if (isempty(str)) || (str(1) == '.'), return; end;
isKey = true;
for k = 1 : size(imgKey,1)
    key = deblank(imgKey(k,:));
    pos = strfind(lower(char(str)),lower(key));
    if ~isempty(pos), return; end;
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