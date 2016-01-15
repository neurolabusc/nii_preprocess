function [prefix] = nii_detrend (old_filename, mask_filename, tissue_filenames, mp_filename, order, tissue_threshold,smoothFWHMmm)
% Removal of trends from a 4D NFTI file with Legendre polynomials.
% Trends to be removed:
% -  Linear (and, optionally, quadratic and cubic)
% -  Mean timecourses of white matter and of CSF
% -  Motion parameters
%
%
% Inputs:
%  old_filename      : name of the input 4D NIFTI file
%  mask_filename     : (optional) name of 3D NIFTI file with the 3D mask
%                             (if empty, the mask is computed using
%                             nonzero voxels in the first volume of the
%                             input)
%  tissue_filenames  : name(s) of 3D NIFTI mask of specific tissue (gray matter, white matter)
%  mp_filename       : name of the text file with motion parameters
%  order             : (optional) detrending order:
%                      1: removal of linear trends
%                      2: removal of linear and quadratic trends
%                      3: (default) removal of linear, quadratic and cubic trends
% tissue_threshold   : (optional) include voxels in tissue masks that
%                      exceed this value (default 0.5)
% smoothFWHMmm       : (optional) if provided and non-zero, image will be
%                      blurred with a Gaussian of given full-width half-maximum
% Outputs: prefix for detrended data, e.g. 'd' for 'fmri.nii' -> 'dfmri.nii'
% Examples
%   nii_detrend('swars.nii',[],'','rp_ars.txt') %detrend motion parameters
%   nii_detrend('swars.nii',[],'wc1t1.nii','') %detrend for Gray Matter
%   nii_detrend('swars.nii',[],strvcat('wc1t1.nii','wc2t1.nii' ),'') %detrend for Gray and White Matter
%   nii_detrend('swars.nii',[],'wc1t1.nii','rp_ars.txt') %detrend for Motion and Gray Matter
%   nii_detrend('swars.nii',[],strvcat('wc1t1.nii','wc2t1.nii' ),'rp_ars.txt') %detrend for Motion, Gray and White Matter
%   nii_detrend('swars.nii',[],strvcat('wc1t1.nii','wc2t1.nii'),'rp_ars.txt',3,0.9) % only gray and white matter voxels > 90% probability
% Chris Rorden and Grigori Yourganov <grigori.yourganov@gmail.com> 11/2013
% GNU General Public Licence (version 2).
if ~exist('old_filename') %no image
    old_filename = spm_select(1,'image','Select 4D data to detrend');
    [pth,nam,ext,vol] = spm_fileparts( old_filename);
    old_filename = fullfile(pth,[ nam, ext]); %'img.nii,1' -> 'img.nii'
end
if ~exist('tissue_filenames') %no image
    tissue_filenames = spm_select(inf,'image','Select tissue image(s)');
end
if ~exist('mp_filename') %no motion parameters specified
    [mp_filename,PathName] = uigetfile('*.txt','Select realignment parameters (rp_*.txt)');
    mp_filename = [PathName mp_filename];
    if (PathName == 0), mp_filename = ''; end; %user pressed cancel
end
if ~exist('order')
    order = 3;
end
if ~exist('tissue_threshold')
    tissue_threshold = 0.5;
end;
if exist(old_filename, 'file') ~= 2
    fprintf('%s error: unable to find NIfTI format image named %s)\n',mfilename,old_filename);
    return;    
end;
if ~exist('smoothFWHMmm')
    smoothFWHMmm = 0;
end
%%%if (length(mp_filename) < 1) && (length(tissue_filenames) < 1) 
%%%    fprintf('%s error: no predictors provided (specify motion parameters and/or tissue maps)\n',mfilename);
%%%    return;    
%%%end;
prefix = 'd'; %'d'etrended data
[pth,nam,ext,vol] = spm_fileparts( old_filename);
new_filename = fullfile(pth,[prefix nam, ext]); %'img.nii' -> 'dimg.nii'
% load the 4D image
nii_hdr = spm_vol (old_filename);
[nii_img] = spm_read_vols (nii_hdr);
n_vols = size (nii_img, 4);
if (n_vols < 2)
    fprintf('%s error: detrending requires a 4D image\n',mfilename);
    return;
end
% look at the first volume
firstvol = nii_img(:, :, :, 1);
volume_size = size (firstvol);
if ~exist('mask_filename') || isempty (mask_filename)
    % if mask is not specified, use the first volume instead
    brain_idx = find (firstvol ~= 0);
else
    % otherwise, use the mask
    mask_hdr = spm_vol (mask_filename);
    mask_img = spm_read_vols (mask_hdr);
    brain_idx = find (mask_img ~= 0);
end
% convert the image into a 2D data matrix; compute mean timecourses of CSF and WM 
img = reshape (nii_img, [], n_vols);
% next create predictor list
nonresidualized_pred = [];
%%%nonresidualized_pred = zeros(1,n_vols);
% add motion parameters as 6 predictors (3 dimensions of translations and rotations) 
if length(mp_filename) > 0
    fprintf('detrending with motion parameters frome %s\n', mp_filename);
    fp = fopen (mp_filename);
    mp = fscanf (fp, '%f ');
    fclose (fp);
    if size(mp,1) ~= (n_vols*6)
        fprintf('%s warning: parameter file %s reports %d values but expected %d (six per volume [for rotation and translation in 3 dimensions])\n',mfilename,mp_filename,size(mp,1), (n_vols*6));
    end
    mp = reshape (mp, [], n_vols);
    %%%nonresidualized_pred = [mp];
    nonresidualized_pred = [mp;  nonresidualized_pred];
end
% add mean signal from each tissue map as predictor
if length(tissue_filenames) > 0
    ntissue = length(tissue_filenames(:,1));
    for tissue = 1 : ntissue
        tissue_hdr = spm_vol (deblank (tissue_filenames(tissue,:)));
        tissue_img = spm_read_vols (tissue_hdr);
        tissue_idx = find (tissue_img > tissue_threshold);
        fprintf('detrending with %d voxels from image %s\n',sum(tissue_idx(:)), deblank (tissue_filenames(tissue,:)));
        tissue_mean_timecourse = mean (img (tissue_idx, :), 1);
        nonresidualized_pred = [tissue_mean_timecourse;  nonresidualized_pred];
    end
end
for k=size(nonresidualized_pred,1):-1:1
	if var(nonresidualized_pred(k,:)) == 0
		nonresidualized_pred(k,:) = [];
        fprintf('%s: parameter %d excluded (no variance)\n',mfilename, k);

	end
end
if (size(nonresidualized_pred,1) < 1)
    fprintf('%s error: no parameters left for detrending.\n',mfilename);
    return;
end
fprintf('%s: detrending for %d parameters\n',mfilename, size(nonresidualized_pred,1));
% create predictors; start with Legendre polynomials
Legendre_poly = zeros (4, n_vols);
x = (-1:2/(n_vols-1):1)';
Legendre_poly (1, :) = 1; 
Legendre_poly (2, :) = x; % linear trend
Legendre_poly (3, :) = 0.5 * (3*x.^2 - 1); % quadratic trend
Legendre_poly (4, :) = 0.5 * (5*x.^3 - 3*x); % cubic trend
n_pred = order + 1;
G = Legendre_poly (1:n_pred, :);
% Legendre polynomials are, by definition, mutually orthogonal.
% The other predictors (WM/CSF mean, 6 motion parameters) need to be
% residualized, that is, made orthogonal to each other and to Legendre
% polynomials. 
for i = 1:size (nonresidualized_pred, 1)
    [e, Beta] = residualize (G, nonresidualized_pred (i, :));
    G = [G; e];
    n_pred = n_pred + 1;
end % now all predictors are mutually orthogonal.
% extract the brain voxels; subtract the mean volume
masked_img = img (brain_idx, :);
brain_mean = mean (masked_img, 2);
masked_img = masked_img - repmat (brain_mean, [1 n_vols]);
% regress out the predictor matrix; add back the mean volume
[e, Beta] = residualize (G, masked_img);
img_detr = masked_img - Beta(:, 2:n_pred) * G(2:n_pred, :) + repmat (brain_mean, [1 n_vols]);
% put the detrended data into 4D image and save to disk
new_hdr = nii_hdr(1);%spm_vol([fnm ',1' ]); %load 4D dataset only once!
new_hdr.fname   = fullfile(pth,[prefix nam ext]);
if smoothFWHMmm > 0
	fprintf('Smoothing %d images with a %.1fmm FWHM Gaussian kernel\n',n_vols,smoothFWHMmm);
    smoothFWHMmm = [smoothFWHMmm smoothFWHMmm smoothFWHMmm];
    VOX = sqrt(sum(new_hdr.mat(1:3,1:3).^2));
    smoothFWHMvox = smoothFWHMmm/VOX; %for 3D arrays the FWHM is specified in voxels not mm     
    new_hdr.descrip = sprintf('%s - conv(%g,%g,%g)',new_hdr.descrip, smoothFWHMmm);
end;

for vol = 1:n_vols
    new_img = nii_img (:, :, :, vol);
    new_img (brain_idx) = img_detr (:, vol);
    new_hdr.n(1)=vol;
    %blur
    if smoothFWHMmm > 0
        presmooth = new_img+0; %+0 forces new matrix
        spm_smooth(presmooth,new_img,smoothFWHMvox,0);
    end %if blur
    spm_write_vol (new_hdr, squeeze (new_img ));
end %for each vol
% end nii_detrend()

function [e, Beta] = residualize (G, y)
% regress out the predictors G from the row vector y.
% model: y = Beta*G + e, where e is orthogonal to all rows of G.
G_pseudoinv = G' * inv (G * G');
Beta = y * G_pseudoinv;
e = y - Beta*G;
% end residualize