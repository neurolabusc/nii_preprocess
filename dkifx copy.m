function dkifx (DWI_nam, bval_nam, Mask_nam, Smooth, Erode)
%Comput mean kurtosis
% DWI_nam  : name of DWI file, either (img.nii or img.nii.gz)
% bval_nam : b-value name (img.bval)
% Mask_nam : name of masking image (mask.nii or mask.nii.gz)
% Smooth : optional, if true or absent data smoothed 1.5voxel FWHM
% Erode    : optional, if true then results are eroded
%instead of computing the diffusion and kurtosis tensor this type of fit
%estimates directly the values of the DKI standards Metrics mean
%diffusivity (MD) and mean kurtosis (MK). This procedure was shown to be
%the most robust procedure to estimate MD/MK (Neto Henriques et al., 2012a; 2012b),
% Neto Henriques, R., Correia, M., Cam-CAN, 2012a. Towards optimization of diffusion Kurtosis imaging to study brain changes with age. Poster presentation at the 29th annual meeting of the European Society for Magnetic Resonance in Medicine and Biology, Lisbon.
% Neto Henriques, R., Ferreira, H., Correia, M., 2012b. Diffusion Kurtosis Imaging of the Healthy Human Brain. Master Dissertation Bachelor and Master program in Biomedical Engineering and Biophysics, Faculty of Sciences, University of Lisbon.
%License: This code is was inspired by DiPy and inherits its license
% Copyright (c) 2016, Chris Rorden, Rafael Neto Henriques
% All rights reserved.
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
% * Redistributions of source code must retain the above copyright
% notice, this list of conditions and the following disclaimer.
% * Redistributions in binary form must reproduce the above
% copyright notice, this list of conditions and the following
% disclaimer in the documentation and/or other materials provided
% with the distribution.
% * The name of the developers nor the names of any
% contributors may be used to endorse or promote products derived
% from this software without specific prior written permission.
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
% "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
% LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
% A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
% OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
% SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
% LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
% DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
% THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
% (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
% OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
%Example
% dkifx; %use GUI
% dkifx('DTI_1k2ku.nii', 'DTI_1k2ku.bval', 'DTI_1k2ku_mask.nii')
% dkifx('sDTI_1k2ku.nii', 'DTI_1k2ku.bval', 'DTI_1k2ku_mask.nii', false);%skip smooting

%check requirements
if isempty(which('spm')) || ~strcmp(spm('Ver'),'SPM12'), error('SPM12 required'); end;
%check inputs - provide GUI if no files selected
if ~exist('DWI_nam','var')  %fnmFA not specified
   [A,Apth] = uigetfile({'*.nii;*.gz;*.hdr;';'*.*'},'Select DTI image');
   DWI_nam = [Apth, A];
end
if ~exist('bval_nam','var') %file not specified
   [A,Apth] = uigetfile({'*.bval';'*.*'},'Select b-value file');
   bval_nam = strcat(Apth,char(A));
end;
if ~exist('Mask_nam','var')  %fnmFA not specified
   [A,Apth] = uigetfile({'*.nii;*.gz;*.hdr;';'*.*'},'Select DTI mask');
   Mask_nam = [Apth, A];
end
%load data
fprintf('%s version 3/2016 loading data\n', mfilename);
[hdr_data, DWI_data] = read_volsSub (DWI_nam); %#ok<ASGLU>
[hdr_mask, DWI_mask] = read_volsSub (Mask_nam);
bval=load(bval_nam);
if numel(bval) ~= size(DWI_data,4)
    error(' Number of bvalues (%d) and volumes (%d) does not match', numel(bval), size(DWI_data,4));
end
%smooth data
if ~exist('Smooth','var') || Smooth  %default: smooth
    fprintf(' Smoothing %d volumes\n', size(DWI_data,4));
    for v = 1 : size(DWI_data,4) %smooth each volume
        %gaussian kernel with 1.5 voxel FWHM (Neto Henriques et al., 2012a; 2012b)
        src = DWI_data(:,:,:,v) + 0; %+0 forces new matrix
        dest = src + 0;
        spm_smooth(src,dest,1.5,0);
        DWI_data(:,:,:,v) = dest;
    end
end
%compute DKI
fprintf(' Computing mean kurtosis\n');
start = tic;
[MK, MD, S0]=DKI_DLS_Sub(DWI_data,DWI_mask,bval);
if exist('Erode','var') && Erode %remove noisy values near border
    [MK, MD, S0]=erodeMaskSub(MK, MD, S0,DWI_mask);
end
% https://github.com/nipy/dipy/blob/master/dipy/reconst/dki.py
% Because high amplitude negative values of kurtosis are not physicaly
%  and biologicaly pluasible, and these causes huge artefacts in kurtosis
%  based measures, directional kurtosis values than `min_kurtosis` are
%  replaced with `min_kurtosis`. defaut = -1
min_kurtosis = -1;
max_kurtosis = 128;
if sum(isfinite(MK(:))) > 0 
    fprintf(' Warning: clipping non-finite MK values (inf, -inf, NaN)\n'); 
    MK(MK == inf) = max_kurtosis;
    MK(MK == -inf) = min_kurtosis;
    MK(~isfinite(MK)) = 0;
end;
fprintf(' %g seconds, Mean kurtosis values range from %g..%g (clipped to %g..%g)\n', toc(start), min(MK(:)), max(MK(:)), min_kurtosis, max_kurtosis );
MK(MK < min_kurtosis) = min_kurtosis;
MK(MK > max_kurtosis) = max_kurtosis;
%save results
fprintf('Saving results\n');
[pth, nam] = filepartsSub(DWI_nam);
save_volSub(fullfile(pth, [nam, '_ldfDKI_MK.nii']), hdr_mask, MK)
save_volSub(fullfile(pth, [nam, '_ldfDKI_MD.nii']), hdr_mask, MD)
save_volSub(fullfile(pth, [nam, '_ldfDKI_S0.nii']), hdr_mask, S0)
%create smoothed, clipped MK-map without negative values: for warping to the individuals normalized, brain extracted T1
MKmask = MK;
%MKmask(~isfinite(MKmask)) = 0; %commented, since we do this with the raw MK
MKmask(MKmask < 0) = 0;
MKmask(MKmask > 4) = 4;
src = MKmask + 0;
spm_smooth(src,MKmask,1.5,0);
save_volSub(fullfile(pth, [nam, '_ldfDKI_MASK.nii']), hdr_mask, MKmask)
%end dkifx()

% --- local sub-functions follow ---

function [MK, MD, S0]=erodeMaskSub(MK, MD, S0,Mask)
%remove values near edge of mask
dims = size(Mask);
img = Mask(:); %convert 3D to 1D vector
img = (img > 0); %binarize 0,1
if (sum(img == 1) < 1) || (sum(img == 0) < 1), error('No variability in image'); end;
imgD = img; %dilated image
imgD = imgD + [img(2:end); 0]; %dilate left
imgD = imgD + [0; img(1:end-1)]; %dilate right
imgD = imgD + [img(dims(1)+1:end); zeros(dims(1),1)]; %dilate anterior
imgD = imgD + [zeros(dims(1),1); img(1:end-dims(1))]; %dilate posterior
sliceVox = dims(1) * dims(2); %voxels per slice
imgD = imgD + [img(sliceVox+1:end); zeros(sliceVox,1)]; %dilate inferior
imgD = imgD + [zeros(sliceVox,1); img(1:end-sliceVox)]; %dilate superior
imgD = (imgD > 6); %binarize ERODED
MK(imgD == 0) = 0;
MD(imgD == 0) = 0;
S0(imgD == 0) = 0;
%end erodeMaskSub()

function [pth,nam,ext,num] = filepartsSub(fname)
% extends John Ashburner's spm_fileparts.m to include '.nii.gz' as ext
num = '';
if ~ispc, fname = strrep(fname,'\',filesep); end
[pth,nam,ext] = fileparts(deblank(fname));
ind = find(ext==',');
if ~isempty(ind)
    num = ext(ind(1):end);
    ext = ext(1:(ind(1)-1));
end
if strcmpi(ext,'.gz')
   [pth nam ext] = fileparts(fullfile(pth, nam));
   ext = [ext, '.gz'];
end
%end filepartsSub()

function save_volSub(fnm, hdr, img)
h = hdr;
h.fname = fnm;
h.dim = hdr.dim(1:3); %only single voluue
h.pinfo = [1;0;0];
h.dt    =[16,0]; %32-bit real datatype
h.private.cal = [0 1];
spm_write_vol(h,img);
%end save_volSub()

function [hdr, img] = read_volsSub (fnm)
[fnm, isGz] = unGzSub (fnm); %convert FSL .nii.gz to .nii
hdr = spm_vol(fnm); %load header data
img = spm_read_vols(hdr); %load image data
if (isGz), delete(fnm); end; %remove .nii if we have .nii.gz
%end read_volsSub()

function [fnm, isGz] = unGzSub (fnm)
[pth,nam,ext] = spm_fileparts(fnm); %#ok<ASGLU>
isGz = false;
if strcmpi(ext,'.gz') %.nii.gz
    fnm = char(gunzip(fnm));
    isGz = true;
end;
%end unGzSub()

function [MK, MD,  S0]=DKI_DLS_Sub(data_in,data_mask,bval)
%Chris Rorden 3/2016
% Directly estimate mean diffusivity (MD) and mean kurtosis (MK).
%  Neto Henriques, R., Correia, M., Cam-CAN, 2012a. Towards optimization of diffusion Kurtosis imaging to study brain changes with age. Poster presentation at the 29th annual meeting of the European Society for Magnetic Resonance in Medicine and Biology, Lisbon.
%  Neto Henriques, R., Ferreira, H., Correia, M., 2012b. Diffusion Kurtosis Imaging of the Healthy Human Brain. Master Dissertation Bachelor and Master program in Biomedical Engineering and Biophysics, Faculty of Sciences, University of Lisbon.
% http://www.cam-can.com/publications/posters/ESMRMB_poster_Henriques.pdf
[Nx, Ny, Nz, Nvol]=size(data_in);
Nspace = Nx * Ny * Nz;
data_ok = reshape(data_in, Nspace, Nvol); %collapse space: 4D (XYZV) -> 2D (SV)
data_ok = data_ok(data_mask(:) > 0, :); %apply mask
A=[(-bval)' ((bval').^2)/6 ones(Nvol, 1)];
B=log(double(data_ok))';           
piA=pinv(A);
X=piA*B;
MD=real(X(1,:));  % Mean Diffusivity            
MK=real(X(2,:)./(MD.^2)); % Mean Kurtosis
S0=real(exp(X(3,:)));    % t2 imaging
%restore spatial shape
MK = restore1Dto3D(MK(:), data_mask);
MD = restore1Dto3D(MD(:), data_mask);
S0 = restore1Dto3D(S0(:), data_mask);
%end DKI_DLS_Sub()

function img3D = restore1Dto3D(img1D, mask) %convert 1D data to 3D, zeroing voxels in mask
img3D = zeros(size(mask));
img3D(mask > 0) = img1D;
%end restore1Dto3D
