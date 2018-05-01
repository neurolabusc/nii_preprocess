function isOK = nii_check_dims(vols)
%Returns error if all filenames exist but have different dimensions
% vols: list of filenames
%Examples
% nii_check_dims; %use GUI
% nii_check_dims({imgs.DTI; imgs.DTIrev})
% nii_check_dims(strvcat('DTI_CT184.nii','DTIrev_CT184.nii'));

isOK = true;
if ~exist('vols', 'var') || isempty(vols) 
    vols = spm_select(inf,'image','Selected images to compare');
end;
if ischar(vols), vols = cellstr(vols); end
if numel(vols) < 2, return; end;
h1 = [];
for i = 1 : numel(vols)
    if isempty(vols{i}), continue; end;
    [p, n, x] = spm_fileparts(deblank(vols{i}));
    if strcmpi(x,'.gz') %.nii.gz
        warning('%s Unable to check if DTI dimensions match (.gz format): "%s"\n', mfilename, vols{i});
        return;
    end;
    if ~exist(vols{i},'file')
        warning('%s Unable to find image: "%s"\n', mfilename, vols{i});
        isOK = false;
        return;   
    end
    h = spm_vol(fullfile(p, [n, x, ',1']));
    if isempty(h1), h1 = h; end;
    if ~all(h.dim == h1.dim)
        fprintf('Dimensions do not match (%dx%dx%d vs %dx%dx%d) "%s" "%s"\n', ...
            h1.dim(1), h1.dim(2), h1.dim(3), ...
            h.dim(1), h.dim(2), h.dim(3), ...
            h1.fname, h.fname); 
        isOK = false;
        return;
    end
end

