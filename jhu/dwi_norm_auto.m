function dwi_norm_auto (pth)
%normalize all dwi.img/dwi.voi pairs in a folder
% pth : folders with DWI images
%Examples
% dwi_norm_auto(pwd)

if ~exist('pth','var')
   pth = uigetdir(pwd,'Select folder with *.nii/*.voi pairs'); 
end
if isempty(pth) || ~isdir(pth), return; end;

fnms = dir(char(fullfile(pth,'*.voi')));
if isempty(fnms), fprintf('No *.voi files in %s\n', pth); return; end;
for f = 1: numel(fnms)
    lesion = fullfile(pth, fnms(f).name);
    [p,n] = fileparts(lesion);
    dwi = fullfile(p, [n, '.nii']);
    if ~exist(dwi,'file'); fprintf('Warning: Unable to find %s\n', dwi); continue; end;
    dwi_norm(dwi, lesion);
end
%end dwi_norm_auto()
