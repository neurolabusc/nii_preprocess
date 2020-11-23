function nii_preproc_scan(baseDir,subj)
% David Reddy 20201110
% front end for nii_preprocess()
% especially on Hyperion

    if ~exist('baseDir','var') || isempty(baseDir)
        error('specify basedir variable');
    end
    if ~exist(baseDir,'dir')
        error('baseDir %s is not a valid directory',baseDir);
    end

    if ~exist('subj','var')
        %***Ignores directories containing '_' symbol
        subjDirs = subFolderSub(baseDir);
        subjDirs = sort(subjDirs);
    else
        subjDirs = {subj};
    end
    
    nSubj = 0;
    for n = 1: size(subjDirs,1)
        fprintf('processing %s\n',subjDirs{n});
        nii_files = dir(fullfile(baseDir,subjDirs{n},'*.nii*'));
        
        imgs = struct();
        
        for m = 1: size(nii_files,1)
%preprocess data from multiple modalities3
% imgs.T1: filename of T1 scan - ONLY REQUIRED INPUT: ALL OTHERS OPTIONAL
% imgs.T2: filename used to draw lesion, if absent lesion drawn on T1
% imgs.lesion : lesion map
% imgs.ASL : pCASL or PASL sequence
% imgs.Rest : Resting state sequence
% DTI : Diffusion scan
% DTIrev : Diffusion scan with reverse phase encoding of 'DTIdum'
% fMRI : fMRI scan
%Examples
% imgs.T1 = 'T1.nii'; imgs.ASL = 'ASL.nii';
% nii_preprocess(imgs);            
            fname = nii_files(m).name;
            pname = fullfile(baseDir,subjDirs{n},nii_files(m).name);
            if (regexp(fname,'^T1_') == 1)
                imgs.T1 = pname;
            elseif (regexp(fname,'^T2_') == 1)
                imgs.T2 = pname;
            elseif(regexp(fname,'^Lesion_') == 1)
                imgs.lesion = pname;
            elseif (regexp(fname,'^ASL_') == 1)
                imgs.ASL = pname;
            elseif (regexp(fname,'^ASLrev_') == 1)
                imgs.ASLrev = pname;
            elseif (regexp(fname,'^Rest_') == 1)
                imgs.Rest = pname;
            elseif (regexp(fname,'^DTI_') == 1)
                imgs.DTI = pname;
            elseif (regexp(fname,'^DTIrev_') == 1)
                imgs.DTIrev = pname;
            elseif (regexp(fname,'^fMRI_') == 1)
                imgs.fMRI = pname;
            elseif (regexp(fname,'^FLAIR_') == 1)
                imgs.FLAIR = pname;
            end;
            %imgs
        end;
        nii_preprocess(imgs);
    end;
end
%nii_preproc_scan()

function nameFolds=subFolderSub(pathFolder)
    dir(pathFolder);
    d = dir(pathFolder);
    isub = [d(:).isdir];
    nameFolds = {d(isub).name}';
    nameFolds = nameFolds(cellfun(@(s)isempty(regexp(s,'^_')),nameFolds)); %remove folders starting with underscores
    nameFolds = nameFolds(cellfun(@(s)isempty(regexp(s,'\.')),nameFolds)); %remove folders with periods
    nameFolds = nameFolds(cellfun(@(s)isempty(regexp(s,' ')),nameFolds)); %remove folders with spaces
    nameFolds(ismember(nameFolds,{'.','..'})) = [];
end
%end subFolderSub()