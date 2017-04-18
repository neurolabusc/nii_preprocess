function nii_preprocess_ct(ct, lesion, matName)
%Normalize CAT scan to standard space, use parameters to warp lesion
% ct : name of CT scan to normalize
% lesion : lesion image in same space as ct
% matName: name of mat file that with processed data.
%Example
% nii_preprocess_ct; %use GUI
% nii_preprocess_ct('CT_M2133.nii','LESION_M2133.nii');
if ~exist('ct','var') ||  ~exist(ct,'file')
 ct = spm_select(1,'image','Select CT to normalize');
end;
if ~exist('lesion','var') ||  ~exist(lesion,'file')
 lesion = spm_select(1,'image','Select matching lesion map');
end;
if ~exist('matName','var')
    [p,n] = spm_fileparts(ct);
    matName = fullfile(p,[n, '.mat']);
end;
if exist('clinical_ctnorm') ~= 2, error('Please install clinical toolbox for spm'); end;
vox = [1 1 1];
bb = [-78 -112 -70; 78 76 85];
vols = {ct, lesion};
nii_setOrigin12(vols, 4, true);
clinical_ctnorm(ct,lesion,vox,bb,true);
wct = prefixSub('w',ct);
wlesion = prefixSub('bws',lesion);
if ~exist(wct,'file') || ~exist(wlesion,'file'), error('Unable to find normalized images "%s" "%s"',wct, wlesion); end;
nii_nii2mat(wlesion, 'lesion', matName); %1



function nam = prefixSub (pre, nam)
if isempty(nam), return; end;
[p, n, x] = spm_fileparts(nam);
nam = fullfile(p, [pre, n, x]);
%end prefixSub()