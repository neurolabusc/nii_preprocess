function fnms = nii_dir(pth, isDir, is1stOnly)
%return sorted list of files or folders
% pth       : search term
% isDir     : if true, return folders, else files
% is1stOnly : return at most one filename
%Note:
% The extension '.n*' will return BOTH .nii and .nii.gz files!
%Examples
% fnms = nii_dir(pwd) %all folders in current working directory
% fnms = nii_dir(fullfile(pwd,'M*')) %all folders that start with 'M'
% fnms = nii_dir(fullfile(pth,'T1_*.nii'), false) %all NIfTI that start with 'T1_'
% fnms = nii_dir(fullfile(pth,'T1_*.n*'), false) %all .nii AND .nii.gz files that start with 'T1_'
% fnms = nii_dir(pwd, true, true) %1st folder in current working directory

if ~exist('isDir', 'var') || (isDir)
    fnms = dirX(pth, true);
else
    [p, n, x] = fileparts(pth);
    isNIfTI = startsWith(x, '.nii', 'IgnoreCase', true);
    if isNIfTI
        pth = fullfile(p,n);
    end
    fnms = dirX(pth, false, isNIfTI);
end
if ~isempty(fnms) && exist('is1stOnly', 'var') && (is1stOnly)
    fnms = fnms(1);
end
%end dirs()

function fnms = dirX(pth, isDir, isNIfTI)
%sorted list of filenames, hidden and empty files removed
%return folders {'/home/chris/M1', '/home/chris/M3'}
% fnms = dirX(/home/chris/M*, true)
%return NIfTIs {'/home/chris/M1.nii', '/home/chris/M3.nii.gz'}
% fnms = dirX(/home/chris/M*, false, true)
%return JSONs
% fnms = dirX(/home/chris/*.json, false, true)
d = dir(pth);
d = d(~startsWith({d.name}, '.'));
if isDir
    isub = [d(:).isdir];
else
    
    isub = ~[d(:).isdir] & ([d(:).bytes] > 0);
    if exist('isNIfTI','var') && isNIfTI
       nii = endsWith({d.name}, '.nii') | endsWith({d.name}, '.nii.gz');
       isub = isub & nii;
    end
end
fnms = {d(isub).name}';
[~,idx] = sort(upper(fnms));
fnms = fnms(idx);
if ~exist(pth, 'dir')
   pth = fileparts(pth); 
end
fnms = strcat(pth, filesep, fnms);
%end dirX()