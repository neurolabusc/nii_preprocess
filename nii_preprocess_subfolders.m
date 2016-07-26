function nii_preprocess_subfolders(pth)
%find all limegui.mat files in subfolders and re-process them
% pth: parent folder
%Examples
% nii_preprocess_subfolders('~/a'); %would process ~/a/b/1_limegui.mat.mat, ~/a/c/2_limegui.mat
% nii_preprocess_subfolders; % search from current working directory
if ~exist('pth','var'), pth = pwd; end;
f = subFolderSub(pth);
if isempty(f), error('No folders in parent folder %s', pdth); end;

global ForcefMRI;  ForcefMRI = true;
global ForceRest;  ForceRest = true;

t = tic;
n = 0;
f = {'M2127'};
for i = 1: numel(f)
   cpth = fullfile(pth,f(i)); %childpath
   mfile = dir(char(fullfile(cpth,'*limegui.mat')));
   if isempty(mfile), continue; end;
   mfile = char(fullfile(cpth, mfile(1).name));
   
   nii_preprocess_gui(mfile);
   n = n + 1;
end
fprintf('Processed %n *limegui.mat file in %gs\n', n, toc(t))

function nameFolds=subFolderSub(pathFolder)
d = dir(pathFolder);
isub = [d(:).isdir];
nameFolds = {d(isub).name}';
nameFolds(ismember(nameFolds,{'.','..'})) = [];
%end subFolderSub()