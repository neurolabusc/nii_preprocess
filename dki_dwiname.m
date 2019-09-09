function [dwi_name, bvec, bval] = dki_dwiname(dti)
%both.bval

dti = unGzNameSub(deblank(dti));
[p,n] = fileparts(dti);
dti = fullfile(p,n);
bvec = [dti 'du.eddy_rotated_bvecs'];
if ~exist(bvec,'file')
    bvec = [dti 'both.bvec'];
end
bval = [dti 'dboth.bval'];
if exist(bvec,'file') && exist(bval,'file') 
    dwi_name = getNiiSub(bvec); 
    return; 
end %combined AP/PA bvecs
%fall back to original bvecs and bvals
if ~exist(bvec,'file')
    bvec = [dti '.bvec'];
end
bval = [dti '.bval'];
if ~exist(bvec,'file') || ~exist(bval,'file')
    error('Can not find files %s %s', bvec, bval); 
end
dwi_name = getNiiSub(bvec);
%end getBVec()

function dwi_name = getNiiSub(bvec)
[p,n] = fileparts(bvec);
dwi_name = fullfile(p, [n,'.nii']);
if exist(dwi_name,'file')
   return; 
end
dwi_name = fullfile(p, [n,'.nii.gz']);
if ~exist(dwi_name,'file')
    error('Can not find files %s', dwi_name); 
end
%end getNii()


 function imgname = unGzNameSub(imgname)
[p, n, x] = fileparts(imgname);
if strcmpi(deblank(x),'.gz') %.nii.gz
    imgname = fullfile(p,n);
end
%end unGzNameSub()