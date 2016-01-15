function nii_dti_prep (dtia,dtib,t1,lesion,newPath)
%Process DTI data
% dtia   : 4D DTI dataset with standard polarity (A->P direction)
% dtib   : 4D DTI dataset with opposite polatiry to dtia (P->A direction)
% t1     : name of T1 scan (previously normalized/segmented)
% lesion : (optional) lesion mask used for T1 norm/seg
%Output: DTI dataset
%
% steps
%  1: topup - undistort, create FA, MD, v1, v2, v3
%  2: create brain extracted T1
%  3: warp ROIs, Gray Matter and White Matter Maps to DTI space
%      -Reverse normalize ROI to native space
%      -Coregister extracted T1 to DTI
%      -Reslice ROI using nearest neighbor interpolation
%      -Reslice GM/WM tissue using trilinear interpolation 
%
% Example
%   nii_dti_prep('DTIA_LM1001.nii','DTIP_LM1001.nii','T1_LM1001.nii','LS_LM1001.nii','');
%   nii_dti_prep('DTIA_LM1001.nii','','T1_LM1001.nii','LS_LM1001.nii','');

if ~exist('lesion')
    lesion = '';
end
if isempty(newPath) 
    newPath = fullfile(pwd, 'temp');
    mkdir(newPath);
end
if (exist('newPath')) && (length(newPath) > 0)
    dtia = cpImgSub(newPath,dtia);
    dtib = cpImgSub(newPath,dtib);
    t1 = cpImgSub(newPath,t1);
    lesion = cpImgSub(newPath,lesion);
end
%roi = 'jhu1mm'; %name for region of interest
[pthm,namm,extm] = spm_fileparts( deblank (which(mfilename)));
[pth,nam,ext] = spm_fileparts(dtia);
[pthb,namb,extb] = spm_fileparts(dtib);
md = fullfile(pth,['v' nam '_MD.nii']); %mean diffusion map
if exist(md,'file')
    fprintf('Skipping topup: file exists named %s\n',md);
else
    %ensure b-vector and b-value files are in the correct folder
    bvec = fullfile(pth,[nam '.bvec']); 
    if exist(bvec) ~= 2 
        src = fullfile(pthm,['DTI.bvec']); 
        copyfile(src,bvec);
        src = fullfile(pthm,['DTI.bval']); 
        bval = fullfile(pth,[nam '.bval']); 
        copyfile(src,bval);
    end 
    bvec = fullfile(pthb,[namb '.bvec']); %mean diffusion map
    if exist(bvec) ~= 2 
        src = fullfile(pthm,['DTI.bvec']); %mean diffusion map
        copyfile(src,bvec);
        src = fullfile(pthm,['DTI.bval']); %mean diffusion map
        bval = fullfile(pthb,[namb '.bval']); %mean diffusion map
        copyfile(src,bval);
    end 
    %run topup to undistort images, compute MD and FA maps
    nii_topup(dtia,dtib,0.03465,2);
    md = fullfile(pth,['vtp' nam '_MD.nii.gz']); %mean diffusion map
    gunzip(md);
    md = fullfile(pth,['vtp' nam '_MD.nii']); %mean diffusion map
end;
%warp regions of interest to DTI data
%resliceROI(t1, md, lesion);
%refROIwarp =  nii_invflirtSub (anat, lesion, ref, refROI);
%end nii_dti_prep()

function newName = cpImgSub(newPath,oldName)
if length(oldName) < 1
   newName = '';
   return;
end
if exist(newPath) == 0
    mkdir(newPath);
end
[oldPath,nam,ext] = spm_fileparts( deblank (oldName));
newName = fullfile(newPath,[nam ext]);
doCpSub(oldPath,newPath,[nam ext]);
doCpSub(oldPath,newPath,['m' nam ext]);
doCpSub(oldPath,newPath,['c1' nam ext]);
doCpSub(oldPath,newPath,['c2' nam ext]);
doCpSub(oldPath,newPath,[nam '_seg_inv_sn.mat']);
doCpSub(oldPath,newPath,[nam '_seg_sn.mat']);
doCpSub(oldPath,newPath,[nam '.bvec']);
doCpSub(oldPath,newPath,[nam '.bval']);
%end cpImgSub()

function doCpSub(oldPath,newPath,namext);
oldName = fullfile(oldPath,namext);
if exist(oldName) ~= 2 
    return;
end
newName = fullfile(newPath,namext);
copyfile(oldName,newName);
%end doCpSub

function img = checkFilenameSub (pth, prefix, nam);
img = fullfile(pth,[prefix nam '.nii']);
if exist(img) ~= 2 
    fprintf('%s warning: unable to find image named %s\n',mfilename,img);
    img = '';
end  
%end checkFilenameSub()

function resliceROI(t1, dti, lesion) 
%coregister t1 to dti, use parameters to reslice c1, c2, tpm
roi = {'jhu1mm'; 'bro1mm'; 'catani1mm'}; %regions of interest
c1c2 = {['c1' t1]; ['c2' t1]};%tissue maps
betT1 = ['render' t1];
template2nativeSub(t1, roi, c1c2, dti, betT1, lesion); 
template2nativeAltSub(dti, roi); 

function template2nativeAltSub(dti, roi) 
[pthm,~,~] = spm_fileparts( deblank (which(mfilename)));
pthRoi = cellAddPth(pthm,'', roi);
spm_jobman('initcfg');
for i = 1: numel(pthRoi)
    %create binarized smoothed version of template
    hdr = spm_vol(char(pthRoi(i)));
    img = spm_read_vols(hdr);
    [spth,snam,sext] = spm_fileparts(char(pthRoi(i)));
    raw_imgName =  fullfile(pwd,[ snam sext]);
    hdr.fname   = raw_imgName;
    spm_write_vol (hdr, img); % write raw image in local directory    
    binImg = zeros(size(img));
    mx = max(img(:));
    binImg(img > 0) = mx;
    spm_smooth(binImg,img,[3 3 3],0); %blur 3-voxel FWHM
    smooth_imgName =  fullfile(pwd,['s' snam sext]);
    hdr.fname   = smooth_imgName;
    spm_write_vol (hdr, img);
    %1 align binarized smoothed template to DTI
    matlabbatch{1}.spm.spatial.coreg.estimate.ref = {[dti ',1']};
    matlabbatch{1}.spm.spatial.coreg.estimate.source = {[smooth_imgName ',1']};
    matlabbatch{1}.spm.spatial.coreg.estimate.other = {[raw_imgName ',1']};
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];
    %2 reslice template to DTI space using nearest neighbor interpolation
    matlabbatch{1}.spm.spatial.coreg.write.ref = {[dti ',1']};
    matlabbatch{1}.spm.spatial.coreg.write.source = {[raw_imgName ',1']};
    matlabbatch{1}.spm.spatial.coreg.write.roptions.interp = 0; %nearest neighbor
    matlabbatch{1}.spm.spatial.coreg.write.roptions.wrap = [0 0 0];
    matlabbatch{1}.spm.spatial.coreg.write.roptions.mask = 0;
    matlabbatch{1}.spm.spatial.coreg.write.roptions.prefix = 'rs';
    spm_jobman('run',matlabbatch);
end

function template2nativeSub(t1, roi, c1c2, dti, betT1, lesion) 
spm_jobman('initcfg');
[pth,nam,ext] = spm_fileparts(t1);
%1 transform regions of interest from standard to native space
[pthm,namm,extm] = spm_fileparts( deblank (which(mfilename)));
pthRoi = cellAddPth(pthm,'', roi);
matlabbatch{1}.spm.util.defs.comp{1}.sn2def.matname ={fullfile(pth, [ nam,'_seg_inv_sn.mat'])};
matlabbatch{1}.spm.util.defs.comp{1}.sn2def.vox = [NaN NaN NaN];
matlabbatch{1}.spm.util.defs.comp{1}.sn2def.bb = [NaN NaN NaN; NaN NaN NaN];
matlabbatch{1}.spm.util.defs.ofname = '';
matlabbatch{1}.spm.util.defs.fnames = pthRoi;
matlabbatch{1}.spm.util.defs.savedir.savepwd = 1;
matlabbatch{1}.spm.util.defs.interp = 0; %nearest neighbor
spm_jobman('run',matlabbatch);
pthRoi = cellAddPth(pwd,'w', roi);%regions of interest now in cwd
%2 use extracted T1 to align c1,c2, and TPM with DTI
matlabbatch{1}.spm.spatial.coreg.estimate.ref = {[dti ',1']};
matlabbatch{1}.spm.spatial.coreg.estimate.source = {[betT1 ',1']};
matlabbatch{1}.spm.spatial.coreg.estimate.other = {c1c2; pthRoi; lesion};
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7]; 
%3 reslice c1, c2 from native to DTI using trilinear interpolation
matlabbatch{2}.spm.spatial.coreg.write.ref = {[dti ',1']};
matlabbatch{2}.spm.spatial.coreg.write.source = [c1c2; {betT1}]; %c1c2;
matlabbatch{2}.spm.spatial.coreg.write.roptions.interp = 1; %linear
matlabbatch{2}.spm.spatial.coreg.write.roptions.wrap = [0 0 0];
matlabbatch{2}.spm.spatial.coreg.write.roptions.mask = 0;
matlabbatch{2}.spm.spatial.coreg.write.roptions.prefix = 'r';
%4 reslice tpm from native to DTI using nearest neighbor interpolation
matlabbatch{3}.spm.spatial.coreg.write.ref = {[dti ',1']};
matlabbatch{3}.spm.spatial.coreg.write.source = pthRoi;
matlabbatch{3}.spm.spatial.coreg.write.roptions.interp = 0; %nearest neighbor
matlabbatch{3}.spm.spatial.coreg.write.roptions.wrap = [0 0 0];
matlabbatch{3}.spm.spatial.coreg.write.roptions.mask = 0;
matlabbatch{3}.spm.spatial.coreg.write.roptions.prefix = 'r';
spm_jobman('run',matlabbatch);

function cellsOut = cellAddPth(pth,prefix,cells)
% cells = {'jhu1mm'} -> cellsOut = {'/home/d/jhu1mm.nii'
for i = 1: size(cells,1)
    cellsOut{i} = fullfile(pth, [prefix char(cells{i}), '.nii']);
    if exist(cellsOut{i},'file') ~= 2 
        error('%s warning: unable to find image named %s\n',mfilename,cellsOut{i});
    end;  
end
%cellAddPth()

function betT1 = extractSub(t1, c1, c2, c3, thresh, PreserveMask)   
%subroutine to extract brain from surrounding scalp
% t1: anatomical scan to be extracted
% c1: gray matter map
% c2: white matter map
% c3: [optional] spinal fluid map
% PreserveMask: [optional] any voxels with values >0 in this image will be spared
%Example
% extractBrain('mT1_LM1000.nii','c1T1_LM1000.nii','c2T1_LM1000.nii','c3T1_LM1000.nii');
if ~exist('t1') %no files
 t1 = spm_select(1,'image','Select T1 image');
end;
if ~exist('c1') %no files
    c1 = spm_select(1,'image','Select gray matter image');
end;
if ~exist('c2') %no files
    c2 = spm_select(1,'image','Select white matter image');
end;
if ~exist('c3') %no files
    c3 = '';%spm_select(1,'image','Select csf image');
end;
if ~exist('PreserveMask') 
    PreserveMask = '';
end
if ~exist('thresh')
    thresh = 0.05;
end
[pth,nam,ext] = spm_fileparts(t1);
%load headers
mi = spm_vol(t1);%bias corrected T1
gi = spm_vol(c1);%Gray Matter map
wi = spm_vol(c2);%White Matter map
%load images
m = spm_read_vols(mi);
g = spm_read_vols(gi);
w = spm_read_vols(wi);
if length(c3) > 0
   ci = spm_vol(c3);%CSF map
   c = spm_read_vols(ci);
   w = c+w; 
end;
w = g+w;
if  (length(PreserveMask) >0)
    mski = spm_vol(PreserveMask);%bias corrected T1
    msk = spm_read_vols(mski);
    w(msk > 0) = 1;
end;
if thresh <= 0
    m=m.*w;
else
    mask= zeros(size(m));
    for px=1:length(w(:)),
      if w(px) >= thresh
        mask(px) = 255;
      end;
    end;
    spm_smooth(mask,mask,1); %feather the edges
    mask = mask / 255;
    m=m.*mask;
end;
fprintf('creating render image based on %s %s %s %s %s\n',t1, c1, c2, c3, PreserveMask);
betT1 = fullfile(pth,['render',  nam, ext]);
mi.fname = betT1;
mi.dt(1) = 4; %16-bit precision more than sufficient uint8=2; int16=4; int32=8; float32=16; float64=64
spm_write_vol(mi,m);
%end extractSub()

function refROIwarp =  nii_invflirtSub (anat, lesion, ref, refROI)
%warp indexed image refROI to space of anat
% anat   : structural scan
% lesion : (optional) lesion mask in space of anat
% ref    : template structural scan
% refROI : indexed region of interest image in space of ref
%example
% nii_invflirt('wT1_LM1001.nii.gz','wLS_LM1001.nii.gz','catanianat.nii','catani1mm.nii');
fsldir= '/usr/local/fsl/';
if ~exist(fsldir,'dir'), error('%s: fsldir (%s) not found',mfilename,fsldir); end
if isempty(lesion)
    mask = '';
else
    %0 make inverted lesion mask
    [pth, nam, ext] = fileparts(lesion); [~, nam] = fileparts(nam); %file.nii.gz -> file.nii -> file
    mask = fullfile(pth,['i' nam '.nii.gz']);
    command=sprintf('sh -c ". ${FSLDIR}/etc/fslconf/fsl.sh; ${FSLDIR}/bin/fslmaths %s -thr 0.5 -binv %s"\n',lesion,mask);
    system(command);
    mask = [' -inweight ' mask];
end;
%1: compute matrix to go from anat -> ref
[pth, nam, ext] = fileparts(anat); [~, nam] = fileparts(nam); %file.nii.gz -> file.nii -> file
mat = fullfile(pth,[nam '.mat']);
setenv('FSLDIR', fsldir);
command=sprintf('sh -c ". ${FSLDIR}/etc/fslconf/fsl.sh; ${FSLDIR}/bin/flirt -in %s -ref %s -omat %s%s -bins 256 -cost corratio -searchrx -45 45 -searchry -45 45 -searchrz -45 45 -dof 12"\n',...
    anat,ref,mat,mask);
system(command);
if false %test of normalization
    [pth, nam, ext] = fileparts(ref); [~, nam] = fileparts(nam); %file.nii.gz -> file.nii -> file
    refwarp = fullfile(pth,['w' nam '.nii.gz']);
    command=sprintf('sh -c ". ${FSLDIR}/etc/fslconf/fsl.sh; ${FSLDIR}/bin/flirt -in %s -ref %s  -out %s -init %s -applyxfm"\n',anat,ref,refwarp, mat);
    system(command);
end
%2: compute inverse transform (ref -> anat)
invmat = fullfile(pth,['inv' nam '.mat']);
command=sprintf('sh -c ". ${FSLDIR}/etc/fslconf/fsl.sh; ${FSLDIR}/bin/convert_xfm -omat %s -inverse %s"\n',invmat,mat);
system(command);
%3 warp ROI to anat (use NEAREST NEIGHBOR INTERPOLATION)
[pth, nam, ext] = fileparts(refROI); [~, nam] = fileparts(nam); %file.nii.gz -> file.nii -> file
refROIwarp = fullfile(pth,['w' nam '.nii.gz']);
command=sprintf('sh -c ". ${FSLDIR}/etc/fslconf/fsl.sh; ${FSLDIR}/bin/flirt -in %s -ref %s  -out %s -init %s -applyxfm"\n',ref, anat,refROIwarp, invmat);
system(command);
% nii_invflirtSub()

