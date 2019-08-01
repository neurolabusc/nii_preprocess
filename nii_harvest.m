function nii_harvest (baseDir)

%baseDir = '/home/crlab/Desktop/testDB';
%outDir = '/home/crlab/Desktop/testIn';

%ENTIRE MASTER DB!
%outDir = '/media/research/MasterChief/Master_In';
%baseDir = '/media/research/MasterChief/Master_DB'; %'/Root'

%POLAR STUDY ONLY
%outDir = '/media/research/POLAREXP/POLAR_Master_In';
%baseDir = '/media/research/POLAREXP/POLAR_Master_Db'; %'/Root'
%outDir = '/media/research/POLAREXP/POLAR_Master_In';
%baseDir = '/media/research/POLAREXP/POLAR_Master_Db';



%Big Data Master Chief Backup
outDir = '/media/coffee/BigData1/MasterChief_12-12-2018/Master_In';
baseDir = '/media/coffee/BigData1/MasterChief_12-12-2018/Master_DB';

%BRIE INTERPERSONAL CORRELATION ONLY
%outDir = '/home/research/Desktop/BRIE_PILOT/Master_IN';
%baseDir = '/home/research/Desktop/BRIE_PILOT/Master_DB';


% outDir = '/home/research/In';
% baseDir = '/home/research/DB';

isExitAfterTable = false; % <- if true, only generates table, does not process data
isPreprocess = false; % <- if true full processing, otherwise just cropping
isReportDims = false; %if true, report dimensions of raw data
reprocessRest = true;
reprocessfMRI = false;
reprocessASL = false;
reprocessDTI = false;
reprocessVBM = false;

%outDir = '/media/UBU/Master_In/';
%baseDir = '/media/UBU/Master_DB/'; %'/Root'

if ~exist('baseDir','var') || isempty(baseDir)
    %baseDir = pwd; %
    baseDir = uigetdir('','Pick folder that contains all subjects');
    
end

%***Ignores directories containing '_' symbol
subjDirs = subFolderSub(baseDir);
subjDirs = sort(subjDirs);

%subjDirs = subjDirs(86:98);  % 1-50 of Polar, rest on other box - RN
%subjDirs = subjDirs(1); % temporary, for testing only!!! -- GY
%subjDirs = { 'M10413'; 'M10463'}; 
%subjDirs = {'M10432'};

modalityKeysVerbose = {'Lesion', 'T1', 'T2', 'DTI_',  'DTIrev', 'ASL', 'Rest_', 'fMRI'}; %DTIREV before DTI!!! both "DTIREV.nii" and "DTI.nii" have prefix "DTI"
modalityDependency = [0, 1, 1,  0, 4, 0, 0, 0]; %T1 and T2 must be from same study as lesion

modalityKeys = strrep(modalityKeysVerbose,'_',''); 
xperimentKeys = {'POLAR','SE', 'LIME', 'CT', 'R01', 'CAT'}; %order specifies priority: 1st item checked first!
%create empty structure
blank = [];
blank.subjName = [];
for i = 1: numel(modalityKeys)
    blank.nii.(modalityKeys{i}) =[];
end;
%1st: acquire data
nSubj = 0;
for s = 1: size(subjDirs,1)%1:nSubjDir2 %(nSubjDir2+1):nSubjDir
    subjName = deblank(subjDirs{s});
    if subjName(1) == '.', continue; end;
    %if (numel(subjName) > 1) && (subjName(2) == '4'), fprintf('SKIPPING %s\n', subjName); continue; end; %ignore folders with underscore, "M2015_needsmatfile"
    if isStringInKeySub (subjName,'_'), continue; end; %ignore folders with underscore, "M2015_needsmatfile"
    subjDir = [baseDir,filesep, subjName]; %no filesep
    %fprintf('%s\n', subjDir);
    nSubj = nSubj + 1;
    imgs(nSubj) = blank;
    imgs(nSubj).subjName = subjName;
    for m = 1:numel(modalityKeysVerbose)
        modality = modalityKeysVerbose{m};
        for x = 1: numel(xperimentKeys)
            xLabel = deblank(xperimentKeys{x}); %e.g. "R01"
            xDir = [subjDir,filesep, xLabel]; %no filesep
            if ~exist(xDir, 'file'), continue; end;
            %fprintf('%s\n', xDir);
            imgs(nSubj) = findImgsSub(imgs(nSubj), xDir, xLabel, modality, m, modalityDependency(m));
            %imgs(nSubj) = findImgsSub(imgs(nSubj), xDir, xLabel, modalityKeysVerbose, modalityDependency);
            %imgs(nSubj) = findImgsSub(imgs(nSubj), xDir, xLabel);
        end
    end
end
fprintf('Found %d subjects in %s\n', nSubj, baseDir);
if nSubj < 1, return; end;
if isReportDims
    reportDimsSub(imgs, nSubj); 
end;
%report results
startTime = tic;
% 1st row: describe values
f = fieldnames(imgs(1).nii);
str = 'n,subj';
for i = 1: numel(f)
   str = sprintf('%s\t%s',str, f{i} );
end
fprintf('%s\n', str);
% subsequent rows: source of images
for s = 1: nSubj
    subj = deblank(imgs(s).subjName);
    subjDir = fullfile(outDir, subj);
    matName = fullfile(subjDir, [subj, '_limegui.mat']);
    imgs(s) = findNovelImgs(subjDir, imgs(s), modalityKeysVerbose);  
    str = [int2str(s), ',', imgs(s).subjName];
    for i = 1: numel(f)
        x = '-';
        if ~isempty(imgs(s).nii.(f{i})) && isfield(imgs(s).nii.(f{i}), 'x')
           x = imgs(s).nii.(f{i}).x;
           if ~imgs(s).nii.(f{i}).newImg, x = ['~', x]; end;
        end
        str = sprintf('%s\t%s',str, x );
    end
    fprintf('%s\n', str);
end
fprintf('Table required %g seconds\n', toc(startTime));
%copy core files to new folder
if isExitAfterTable 
    fprintf('Disable isExitAfterTable for full analyses\n', str); 
    return; %return when we are done   
end
if exist(outDir, 'file') ~= 7, error('Unable to find folder %s', outDir); end;
%find images we have already processed
if isempty(spm_figure('FindWin','Graphics')), spm fmri; end; %launch SPM if it is not running
process1st = true;

for s =  1: nSubj 
    anyNewImg = false;
    subj = deblank(imgs(s).subjName);
    subjDir = fullfile(outDir, subj);
    if ~isfield(imgs(s).nii.T1,'img')
        fprintf('Skipping %s: no T1!\n', subj);
        continue;
    end
    %imgs(s) = findNovelImgs(subjDir, imgs(s), modalityKeysVerbose);
    global ForcefMRI;
    global ForceRest;
    global ForceASL;
    global ForceDTI;
    global ForceVBM;
    ForcefMRI=[];
    ForceRest=[];
    ForceASL=[];
    ForceDTI =[];
    ForceVBM = [];
    %666x - 
    %imgs(s).nii.fMRI.newImg = false;
    %imgs(s).nii.Rest.newImg = false;
    %666x <-
    if imgs(s).nii.fMRI.newImg, ForcefMRI = true; end;
    if imgs(s).nii.Rest.newImg, ForceRest = true; end;
    if imgs(s).nii.ASL.newImg, ForceASL = true; end;
    %666 if imgs(s).nii.DTI.newImg, ForceDTI = true; end;
      %to reprocess one modality for EVERYONE....
    if reprocessDTI && isfield(imgs(s).nii.DTI,'img')
        ForceDTI = true;
        anyNewImg = true;
    end
    if reprocessRest && isfield(imgs(s).nii.Rest,'img')
        ForceRest = true;
        anyNewImg = true;
    end
    if reprocessfMRI && isfield(imgs(s).nii.fMRI,'img')
        ForcefMRI = true;
        anyNewImg = true;
    end
    if reprocessASL && isfield(imgs(s).nii.ASL,'img')
        ForceASL = true;
        anyNewImg = true;
    end
    
     if reprocessVBM && isfield(imgs(s).nii.T1,'img')
        ForceVBM = true;
        anyNewImg = true;
    end
    
    %if imgs(s).nii.DTI.newImg, ForceDTI = true; end;
    
    %fprintf('%s %d\n', f{i}, imgs.nii.(f{i}).newImg);
    if exist(subjDir,'file') == 0, mkdir(subjDir); end;
    matName = fullfile(subjDir, [subj, '_limegui.mat']);
    mat = [];
    
    if exist(matName,'file'), mat = load(matName); end;
    
    for i =  1:numel(f)
        if ~isempty(imgs(s).nii.(f{i})) && (imgs(s).nii.(f{i}).newImg)
            m = f{i}; % modality: T1, T2..
            if ~isfield(imgs(s).nii.(f{i}),'img')
                mat.(m) = ''; %e.g. we used to have DTI+DTIrev, now we have DTI
                %warning('Expected field "img" for participant %s modality %s', subj, m);
                %error('123'); %fix this if it ever happens again
            elseif isempty(imgs(s).nii.(f{i}).img)
                mat.(m) = ''; %e.g. we used to have DTI+DTIrev, now we have DTI
            else
                anyNewImg = true;
                m = f{i}; % modality: T1, T2..
                x = imgs(s).nii.(f{i}).x; %e.g. experiment name "LIME", "CT"
                imgin = imgs(s).nii.(f{i}).img; %e.g. '~/dir/m2000/CT/T1.nii'
                imgout = fullfile(subjDir, sprintf('%s_%s_%s.nii',m, subj, x));
                fprintf('%s -> %s\n',imgin, imgout);
                moveImgUnGz(imgin, imgout);
                mat.(m) = imgout;
            end
        end
        if ~isempty(imgs(s).nii.(f{i})) && isfield(imgs(s).nii.(f{i}), 'x') && (~imgs(s).nii.(f{i}).newImg) %CR 2/2017: in case folder names have changed
                m = f{i}; % modality: T1, T2..
                x = imgs(s).nii.(f{i}).x; %e.g. experiment name "LIME", "CT"
                imgout = fullfile(subjDir, sprintf('%s_%s_%s.nii',m, subj, x));
                mat.(m) = imgout;
        end
    end
    
    if anyNewImg
        matNameGUI = fullfile(subjDir, [subj, '_limegui.mat']);
        fprintf('Creating %s\n',matNameGUI);
        save(matNameGUI,'-struct', 'mat');
        %determine T1 name even if output folder renamed...
        if imgs(s).nii.T1.newImg || imgs(s).nii.Lesion.newImg, setAcpcSubT1(matNameGUI); end;
        %if imgs(s).nii.T1.newImg, setAcpcSubT1(matNameGUI); end;
        %if imgs(s).nii.Lesion.newImg, setAcpcSubT1 (matNameGUI); end;
        
        if imgs(s).nii.DTI.newImg, setAcpcSubDTI (matNameGUI); end;

        %process the data
        nii_preprocess(mat,[],process1st);
        process1st = false; %only check for updates for first person
        %matName = fullfile(subjDir, sprintf('T1_%s_%s_lime.mat', subj, imgs(s).nii.T1.x));
        if isPreprocess
            nii_preprocess(mat,matName)
        else
            fprintf('Cropped but did not preprocess %s\n',matName);
        end
        
    end
end
fprintf('All done\n');
%end nii_harvest

function reportDimsSub(imgs,nSubj)
for s = 1: nSubj
    subj = deblank(imgs(s).subjName);
    f = fieldnames(imgs(1).nii);
    fprintf('\nID\tMRI\tstudy\tx\ty\tz\tvols\tTR\n');
    for i = 1: numel(f)
         x = '-';
        if ~isempty(imgs(s).nii.(f{i})) && isfield(imgs(s).nii.(f{i}), 'x') && exist(imgs(s).nii.(f{i}).img,'file')
            fnm = imgs(s).nii.(f{i}).img;
            hdr = readNiftiHdrSub(fnm);
            tr = 0; 
            if isfield (hdr(1).private, 'timing')
               tr =hdr(1).private.timing.tspace 
            end
            fprintf('%s\t%s\t%s\t%d\t%d\t%d\t%d\t%g\n', subj, f{i}, imgs(s).nii.(f{i}).x, hdr(1).dim(1),hdr(1).dim(2),hdr(1).dim(3),numel(hdr),tr  );
        end
    end
    
    
%     subjDir = fullfile(outDir, subj);
%     matName = fullfile(subjDir, [subj, '_limegui.mat']);
%     imgs(s) = findNovelImgs(subjDir, imgs(s), modalityKeysVerbose);  
%     str = [int2str(s), ',', imgs(s).subjName];
%     for i = 1: numel(f)
%         x = '-';
%         if ~isempty(imgs(s).nii.(f{i})) && isfield(imgs(s).nii.(f{i}), 'x')
%            x = imgs(s).nii.(f{i}).x;
%            if ~imgs(s).nii.(f{i}).newImg, x = ['~', x]; end;
%         end
%         str = sprintf('%s\t%s',str, x );
%     end
%     fprintf('%s\n', str);
end
%end reportDimsSub

function tf = endsWithSub(str,pattern) %endsWith only in Matlab 2016 and later
tf = false;
if numel(str) < numel(pattern), return; end;
strEnd = str(end-numel(pattern)+1:end);
tf = strncmpi(strEnd,pattern, numel(pattern));
%endsWithSub


function imgs = findNovelImgs(subjDir, imgs, modalityKeysVerbose)
f = fieldnames(imgs.nii);
for i = 1: numel(f)
        imgs.nii.(f{i}).newImg = true;%??
end
%'fMRI'
if ~isfield(imgs.nii,'T1') || isempty(imgs.nii.T1), return; end;
matname = dir(fullfile(subjDir,'T1_*_lime.mat'));
if isempty(matname), return; end;
matname = fullfile(subjDir,matname(1).name);
m = load(matname);
if isfield(m,'T1') 
    imgs.nii.T1.newImg = ~endsWithSub(m.T1.hdr.fname, ['_',imgs.nii.T1.x,'.nii']);
    imgs.nii.T2.newImg = imgs.nii.T1.newImg;
    imgs.nii.Lesion.newImg = imgs.nii.T1.newImg;
end;
if isfield(m,'cbf') && isfield(imgs.nii.ASL, 'x') 
   % m.cbf.hdr.fname
 
    imgs.nii.ASL.newImg = ~endsWithSub(m.cbf.hdr.fname, ['_',imgs.nii.ASL.x,'M0CSF.nii']);
end;
if isfield(m,'fa') && isfield(imgs.nii.DTI, 'x') 
    imgs.nii.DTI.newImg = ~endsWithSub(m.fa.hdr.fname, ['_',imgs.nii.DTI.x,'d_FA.nii']);
    imgs.nii.DTIrev.newImg = imgs.nii.DTI.newImg;
end;
if isfield(m,'RestAve') && isfield(imgs.nii.Rest, 'x')  
    imgs.nii.Rest.newImg = ~endsWithSub(m.RestAve.hdr.fname, ['_',imgs.nii.Rest.x,'.nii']);
end;
if isfield(m,'fMRIave') && isfield(imgs.nii.fMRI, 'x')
    imgs.nii.fMRI.newImg = ~endsWithSub(m.fMRIave.hdr.fname, ['_',imgs.nii.fMRI.x,'.nii']);
end;

%for i = 1: numel(f)
%    fprintf('%s %d\n', f{i}, imgs.nii.(f{i}).newImg);
%end
%end findNovelImgs()

function setAcpcSubT1 (matname)
m = load(matname);
if isfield(m,'T2') && isfield(m,'Lesion') && ~isempty(m.T2) && ~isempty(m.Lesion)
    nii_setOrigin12({m.T2,m.Lesion}, 2, true); %T2
    if isfield(m,'T1') && isfield(m,'Lesion') 
        nii_setOrigin12({m.T1}, 1, true); %T1 - crop with lesion
    end
    return;
end
if isfield(m,'T1') && isfield(m,'Lesion') && ~isempty(m.T1) && ~isempty(m.Lesion)
    nii_setOrigin12({m.T1,m.Lesion}, 1, true); %T1 - crop with lesion
    return;
end
if isfield(m,'T1') && ~isempty(m.T1)
    nii_setOrigin12(m.T1, 1, true); %T1 - crop
end

function setAcpcSubDTI (matname)
m = load(matname);
if isfield(m,'DTI') && isfield(m,'DTIrev') && ~isempty(m.DTI) && ~isempty(m.DTIrev)
    nii_setOrigin12({m.DTI,m.DTIrev}, 3, false); %DTI
elseif isfield(m,'DTI') && ~isempty(m.DTI)
    nii_setOrigin12(m.DTI, 3, false); %DTI
end
%end setAcpcSub();

function moveImgUnGz(inname, outname)
[ipth, inam,iext] = fileparts(inname);
[opth, onam,oext] = fileparts(outname);
%load data
if strcmpi(iext,'.gz') %unzip compressed data
	inname = gunzip(inname);
    inname = deblank(char(inname));
    [ipth, inam,iext] = fileparts(inname);
end;
copyfile(inname, outname);
if strcmpi(iext,'.gz') %fsl can not abide with coexisting img.nii and img.nii.gz
	delete(filename);
end;
%copy bvec
ibvec = fullfile(ipth, [inam, '.bvec']);
if exist(ibvec, 'file'),
    obvec = fullfile(opth, [onam, '.bvec']);
    copyfile(ibvec, obvec);
end;
%copy bval
ibval = fullfile(ipth, [inam, '.bval']);
if exist(ibval, 'file'),
    obval = fullfile(opth, [onam, '.bval']);
    copyfile(ibval, obval);
end;
%end moveImgUnGz()


function imgs = findImgsSub(imgs, xDir, xLabel, modalityKey, modalityNum, modalityDependency)
nameFiles = subImgSub(xDir);
if isempty(nameFiles), return; end;
%nameFiles = sort(nameFiles); %take first file for multiimage sequences, e.g. ASL
%nameFiles
%nameFiles(1)
f = fieldnames(imgs.nii);
if ~isempty(imgs.nii.(f{modalityNum})), return; end;
if exist('modalityDependency','var') && (modalityDependency ~= 0)
    dep = f{modalityDependency}; %e.g. T1 scan may depend on Lesion
    if ~isempty(imgs.nii.(dep))
        x = imgs.nii.(dep).x;
        if ~strncmpi(xLabel,x, numel(xLabel))
            %fprintf('"%s" must be from same experiment as "%s" (%s not %s)\n', modalityKey, dep, x, xLabel); 
            return;
        end;
        %fprintf('"%s" must be from same experiment as "%s" (%s)\n', modalityKey, dep, x);
    end;
end
for j = 1: numel(nameFiles)
    if strncmpi(modalityKey,nameFiles(j), numel(modalityKey))
       fname = fullfile(xDir, char(nameFiles(j)) );
       %fprintf('%d %s %s %s\n', i, char(f{i}), char(modalityKey), char(nameFiles(j)) );
       imgs.nii.(f{modalityNum}).x = xLabel;
       imgs.nii.(f{modalityNum}).img = fname;
       break;
    end
end;

function nVol = nVolSub(filename)
%report number of volumes
nVol = 0;
[hdr, img] = readNiftiSub(filename);
if isempty(hdr), return; end;
nVol = numel(hdr);
%end nVolSub

function hdr = readNiftiHdrSub(filename)
[p, n,x] = fileparts(filename);
if strcmpi(x,'.gz') %unzip compressed data
	filename = gunzip(filename);
    filename = deblank(char(filename));
    hdr = spm_vol(filename);
    error(fnm);
    delete(fnm);
    return;
end;
hdr = spm_vol(filename);

function [hdr, img] = readNiftiSub(filename)
%load NIfTI (.nii, .nii.gz, .hdr/.img) image and header
% filename: image to open
% open4d: if true all volumes are loaded
%To do:
%  endian: rare, currently detected and reported but not handled
%Examples
% hdr = nii_loadhdrimg('myimg.nii');
% [hdr, img] = nii_loadhdrimg('myimg.nii');
% [hdr, img] = nii_loadhdrimg('img4d.nii');
if ~exist('filename','var')  %fnmFA not specified
   [A,Apth] = uigetfile({'*.nii;*.gz;*.hdr;';'*.*'},'Select image');
   filename = [Apth, A];
end
[fpth, fnam,fext] = fileparts(filename);
if strcmpi(fext,'.img') %hdr/img pair
    filename = fullfile(fpth, [fnam, '.hdr']);
end
if ~exist(filename, 'file')
    error('Unable to find file %s', filename);
end
%load data
if strcmpi(fext,'.gz') %unzip compressed data
	filename = gunzip(filename);
    filename = deblank(char(filename));
end;
hdr = spm_vol(filename);
if hdr(1).dt(1) == 128
   fprintf('Skipping RGB image %s\n', filename);
   hdr = [];
   img = [];
   return;
end
img = spm_read_vols(hdr);
if strcmpi(fext,'.gz') %fsl can not abide with coexisting img.nii and img.nii.gz
	delete(filename);
end;
%end nii_loadhdrimg()

%end findImgsSub()

function nameFiles=subImgSub(pathFolder)
nameFiles=subFileSub(pathFolder);
if isempty(nameFiles), return; end;
n = nameFiles; nameFiles = [];
for i = 1: numel(n)
    [~,~,x] = fileparts(char(deblank(n(i))));
    if ~strncmpi('.gz',x, 3) && ~strncmpi('.nii',x, 4), continue; end;
    nameFiles = [nameFiles; n(i)]; %#ok<AGROW>
end
%end subFileSub()


function nameFiles=subFileSub(pathFolder)
d = dir(pathFolder);
isub = ~[d(:).isdir];
nameFiles = {d(isub).name}';
%end subFileSub()

function isKey = isStringInKeySub (str, imgKey)
isKey = true;
for k = 1 : size(imgKey,1)
    key = deblank(imgKey(k,:));
    pos = strfind(lower(char(str)),lower(key));
    if ~isempty(pos), isKey = pos(1); return; end;
end
isKey = false;
%isStringInKey()

function nameFolds=subFolderSub(pathFolder)
d = dir(pathFolder);
isub = [d(:).isdir];
nameFolds = {d(isub).name}';
nameFolds = nameFolds(cellfun(@(s)isempty(regexp(s,'_')),nameFolds)); %remove folders with underscores
nameFolds = nameFolds(cellfun(@(s)isempty(regexp(s,'\.')),nameFolds)); %remove folders with periods
nameFolds = nameFolds(cellfun(@(s)isempty(regexp(s,' ')),nameFolds)); %remove folders with spaces
nameFolds(ismember(nameFolds,{'.','..'})) = [];
%end subFolderSub()