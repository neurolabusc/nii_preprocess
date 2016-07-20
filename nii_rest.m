function nii_rest (imgs, TRsec, SliceOrder, doSliceTime)
%preprocesses resting state data using SPM
%  imgs : (optional) structure of image names (imgs.T1, imgs.Rest)
%  TRsec      : Time between volumes in seconds
%  SliceOrder : see nii_sliceTime 0=auto-detect,1=ascend,2=descend
%Examples
% nii_rest; %use GUI
% imgs.T1 = 'T1.nii'; imgs.Rest = 'Rest.nii'; imgs.SE = 'RestSE.nii'; imgs.SErev = 'RestSErev.nii';
% nii_rest(imgs);
%
% History:
%
% TH 03/16: now pass in imgs structure rather than separate images for T1, rest, spin echo, and spin echo rev. 

if isempty(which('spm')) || ~strcmp(spm('Ver'),'SPM12'), error('SPM12 required'); end;
if ~exist('nii_batch12','file'), error('Make sure nii_batch12 is in path'); end;
if ~exist('imgs','var') %no input: select imgs[s]
 imgs.Rest = spm_select(inf,'image','Select 4D resting state volumes[s]');
 imgs.T1 = spm_select(1,'image','Select T1 anatomical scan');
 imgs.SE = spm_select(1,'image','(Optional) select spin-echo scan for undistortion');
 imgs.SE = spm_select(1,'image','(Optional) select reversed-phase spin-echo scan for undistortion');
end
clear matlabbatch
tic; %start timer
%p.t1name = imgs.T1; %'t1.nii'
p.t1name = -1; %'t1.nii'

p.TRsec = 0; %repeat time off 10 seconds
p.slice_order = 0;

for ses = 1 : length(imgs.Rest(:,1));
    p.fmriname = deblank(imgs.Rest(ses,:));
    %[pth,nam,ext] = spm_fileparts( deblank (fMRIname));
    [prefix, TRsec, slice_order] = nii_batch12(p);
    prefix
    TRsec
    slice_order
end; %for each session
fprintf('Done processing sessions in %0.3fsec\n', toc);
%end nii_rest()

