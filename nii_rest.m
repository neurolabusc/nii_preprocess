function nii_rest (imgs, TRsec, SliceOrder)
%preprocesses resting state data using SPM
%  imgs : structure of image names (imgs.T1, imgs.Rest)
%  TRsec      : Time between volumes in seconds
%  SliceOrder : see nii_sliceTime -1=no, 0=auto-detect,1=ascend,2=descend
%  doSliceTime : compute slice timing (true/false)
%Examples
% %linear Rest->T1, nonlinear T1->MNI: use when good linear fit between Rest and T1
% imgs.T1 = 'T1_M2129_LIME.nii';
% imgs.Rest = 'Rest_M2129_LIME copy.nii';
% nii_rest(imgs);
% %nonlinear Rest->MNI: use when Rest and T1
% imgs.T1 = [];
% imgs.Rest = 'Rest_M2129_LIME copy.nii';
% nii_rest(imgs);


if isempty(which('spm')) || ~strcmp(spm('Ver'),'SPM12'), error('SPM12 required'); end;
if ~exist('nii_batch12','file'), error('Make sure nii_batch12 is in path'); end;
if ~exist('imgs','var') %no input: select imgs[s]
 imgs.Rest = spm_select(inf,'image','Select 4D resting state volumes[s]');
 imgs.T1 = spm_select(1,'image','(Optional) Select T1 anatomical scan');
 imgs.SE = spm_select(1,'image','(Optional) select spin-echo scan for undistortion');
 imgs.SE = spm_select(1,'image','(Optional) select reversed-phase spin-echo scan for undistortion');
end
if ~exist('TRsec','var') %no input: select imgs[s]
    TRsec = 0;
    fprintf('Will attempt to detect TR from image file.\n');
end
if ~exist('SliceOrder','var') %no input: select imgs[s]
    SliceOrder = 0;
    fprintf('Will attempt to detect Slice Order from image file.\n');
end
clear matlabbatch
tic; %start timer
p.setOrigin = true;
p.t1name = imgs.T1; %'t1.nii'
%p.t1name = -1; %'t1.nii'

p.TRsec = TRsec; %repeat time off 10 seconds
p.slice_order = SliceOrder;
p.slice_order = -1;
for ses = 1 : length(imgs.Rest(:,1));
    p.fmriname = deblank(imgs.Rest(ses,:));
    %[pth,nam,ext] = spm_fileparts( deblank (fMRIname));
    [prefix, TRsec, slice_order] = nii_batch12(p);
    prefix
    TRsec
    slice_order
    error('123')
end; %for each session
fprintf('Done processing sessions in %0.3fsec\n', toc);
%end nii_rest()

