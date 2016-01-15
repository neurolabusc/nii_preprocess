function nii_findmat (baseDir)
%find all the mat files
if ~exist('baseDir','var')
    baseDir = pwd;
end
files=dir(fullfile(baseDir,'*.mat'));
fprintf('Found %d mat files (subjects) in %s\n',numel(files), baseDir);
for f = 1: numel(files)
    fnm = fullfile(baseDir, files(f).name);
    m = load(fnm);
    if ~isfield(m,'dtifc_jhu') && isfield(m,'md')
       fprintf('%s missing dti\n', fnm); 
    end
end
%end nii_findmat()

