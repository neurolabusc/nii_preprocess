function nii_makemask(template, thresh, outname)
%Create a binary mask of voxels in template that exceed threshold
% template : 4D NIfTI image (volume 1=graymatter, volume 2 = white matter)
% thresh   : threshold applied to template
% outname  : 
if ~exist('template','var'), template = fullfile(spm('Dir'),'tpm','TPM.nii'); end;
if ~exist('thresh','var'), thresh = 0.5; end;
if ~exist('outname','var'), 
    p = fileparts(which(mfilename));
    outname = fullfile (p,'mask.nii');
end;
if ~exist(template,'file'), error('Is SPM12 installed? Unable to find template %s',template); end;
hdr = spm_vol(template);
img = spm_read_vols(hdr);
if size(img,4) < 2, error('Expected 4D template'); end;
mhdr = hdr(1);
mhdr.fname = outname;
mimg = ((img(:,:,:,1) + img(:,:,:,2)) > thresh);
spm_write_vol(mhdr,mimg);