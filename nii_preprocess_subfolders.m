function nii_preprocess_subfolders(pth)
%find all limegui.mat files in subfolders and re-process them
% pth: parent folder
%Examples
% nii_preprocess_subfolders('~/a'); %would process ~/a/b/1_limegui.mat.mat, ~/a/c/2_limegui.mat
% nii_preprocess_subfolders; % search from current working directory
%Notes: will skip folders with "_" in name, e.g. will skip M2002:
% pth/M2000 pth/M2001 pth/M2002_dementia pth/M2003
if ~exist('pth','var'), pth = pwd; end;
f = subFolderSub(pth);
if isempty(f), error('No folders in parent folder %s', pdth); end;
global ForcefMRI;  ForcefMRI = true; warning('FORCED fMRI REPROCESSING'); %comment line for auto-processing
global ForceRest;  ForceRest = true; warning('FORCED REST REPROCESSING'); %comment line for auto-processing

t = tic;
n = 0;
%f = {'M2001'}; %for a single folder
for i = 1: numel(f)
   cpth = char(f(i)); %local child path
   if ~isempty(strfind(cpth,'_'))
      fprintf('Warning: "_" in folder name: skipping %s\n', char(cpth) );
      continue
   end
   cpth = char(fullfile(pth,f(i))); %full child path
   nii_preprocess_gui(cpth);
   n = n + 1;
end
fprintf('Processed %d *limegui.mat file in %gs\n', n, toc(t))

function nameFolds=subFolderSub(pathFolder)
d = dir(pathFolder);
isub = [d(:).isdir];
nameFolds = {d(isub).name}';
nameFolds(ismember(nameFolds,{'.','..'})) = [];
%end subFolderSub()