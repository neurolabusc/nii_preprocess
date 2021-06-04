function pth = nii_mrtrix_pth
%return path for MRtrix tools: Matlab does not source the bash PATH
%returns empty if path not found

pth = '/home/shared/miniconda/bin'; % <- customize: manually set
if isTrix(pth), return; end
pth = '/usr/local/bin';
if isTrix(pth), return; end
pth = char(java.lang.System.getProperty('user.home'));
pth = fullfile(pth, 'mrtrix3','bin');
if isTrix(pth), return; end
warning('Unable to find MRtrix3, please customize "%s"\n', mfilename);
pth = []; %failed to find mrtrix   
%end nii_mrtrix_pth()

function ret = isTrix(pth)
ret = exist(fullfile(pth,'mrdegibbs'), 'file');
%end isTrix()
