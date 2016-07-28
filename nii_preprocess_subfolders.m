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
global ForceRest;  ForceRest =  warning('FORCED REST REPROCESSING'); %comment line for auto-processing
%global ForceASL;  ForceASL = true; warning('FORCED ASL REPROCESSING'); %comment line for auto-processing

t = tic;
n = 0;
%f = {'M2082'}; %for a single folder
for i = 1: numel(f) %change 1 to larger number to restart after failure
   cpth = deblank(f(i)); %local child path 
   fprintf('===\t%s participant %d : %s\t===\n', mfilename, i, cpth);
   if ~isempty(strfind(cpth,'M2082')), error('all done'); end; %to stop at specific point
   if ~isempty(strfind(cpth,'_'))
      fprintf('Warning: "_" in folder name: skipping %s\n', char(cpth) );
      continue
   end
   cpth = char(fullfile(pth,cpth); %full child path
   nii_preprocess_gui(cpth);
   n = n + 1;
end
fprintf('Processed %d *limegui.mat file in %gs\n', n, toc(t))
%end nii_preprocess_subfolders()




function nameFolds=subFolderSub(pathFolder)
d = dir(pathFolder);
isub = [d(:).isdir];
nameFolds = {d(isub).name}';
nameFolds(ismember(nameFolds,{'.','..'})) = [];
%end subFolderSub()