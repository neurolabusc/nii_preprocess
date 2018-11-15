function nii_harvest_fast

baseDir = '/media/research/FAT1000CR/Master_In';
isExitAfterTable = true; % <- if true, only generates table, does not process data

global ForcefMRI;
global ForceRest;
global ForceASL;
global ForceDTI;
ForcefMRI=[];
ForceRest=[];
ForceASL=[];
ForceDTI =true;
    
if ~exist('baseDir','var') || isempty(baseDir)
    baseDir = uigetdir('','Pick folder that contains all subjects');
end
subjDirs = subFolderSub(baseDir);
subjDirs = sort(subjDirs);
%subjDirs = {'M4119'};
process1st = true;
%{'M4215','M4220','M4223','M4224'};
subjDirs = {'M4217';'M4218';'M2001'};
try
    for s = 1: size(subjDirs,1)%1:nSubjDir2 %(nSubjDir2+1):nSubjDir
        subjName = deblank(subjDirs{s});
        if  subjName(1) == '.', continue; end;
        subjDir = fullfile(baseDir, subjName);
        matName = fullfile(subjDir, [subjName, '_limegui.mat']);
        if ~exist(matName,'file'), continue; end;
        fprintf('>>>>>>>\t%d\t%s\n', s, matName);
        if isExitAfterTable, continue; end;
        mat = load(matName);
        %nii_preprocess(mat,[],process1st);
        nii_preprocess_gui(matName);
        process1st = false; %only check for updates for first person

    end
    doneSmsSub('OK', mfilename);
catch catcherror 
    doneSmsSub('Error', [catcherror.stack(1).file '->' catcherror.message]);
end
%end main

function doneSmsSub(title,msg)
p = Pushbullet();%api-key set in ~/Documents/Matlab/Pushbullet.m
p.pushNote([],title,msg);
%end doneSmsSub()

function nameFolds=subFolderSub(pathFolder)
d = dir(pathFolder);
isub = [d(:).isdir];
nameFolds = {d(isub).name}';
nameFolds = nameFolds(cellfun(@(s)isempty(regexp(s,'_')),nameFolds)); %remove folders with underscores
nameFolds = nameFolds(cellfun(@(s)isempty(regexp(s,'\.')),nameFolds)); %remove folders with periods
nameFolds = nameFolds(cellfun(@(s)isempty(regexp(s,' ')),nameFolds)); %remove folders with spaces
nameFolds(ismember(nameFolds,{'.','..'})) = [];
%end subFolderSub()