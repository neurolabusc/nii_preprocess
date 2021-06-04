function MK_nam = dkifx2 (DWI_nam, bval_nam, Mask_nam, Smooth)
% DWI_nam  : name of DWI file, either (img.nii or img.nii.gz)
% bval_nam : b-value name (img.bval)
% Mask_nam : (optional) name of masking image (mask.nii or mask.nii.gz)
%Example
% dkifx2; %use GUI
% dkifx2('ddki_ecc.nii.gz', 'dki.bval', 'ddki_mask.nii.gz', false, true)
%check requirements
isDeleteTempImages = false; %if true, we will delete intermediate images
if isempty(which('spm')) || ~strcmp(spm('Ver'),'SPM12'), error('SPM12 required'); end;
%check inputs - provide GUI if no files selected
if ~exist('DWI_nam','var')  %fnmFA not specified
   fprintf('Select DTI image\n');
   [A,Apth] = uigetfile({'*.nii;*.gz;*.hdr;';'*.*'},'Select DTI image');
   if isnumeric(A), return; end;
   DWI_nam = [Apth, A];
end
if ~exist('bval_nam','var') %file not specified
   fprintf('Select b-value file\n');
   [A,Apth] = uigetfile({'*.bval';'*.*'},'Select b-value file');
   if isnumeric(A), return; end;
   bval_nam = strcat(Apth,char(A));
end;
if (nargin < 1) && ~exist('Mask_nam','var')
    fprintf('Select mask image (optional)\n');
   [A,Apth] = uigetfile({'*.nii;*.gz;*.hdr;';'*.*'},'Select mask image (optional)');
   if ~isnumeric(A), Mask_nam = [Apth, A]; end;    
end
if ~exist('Mask_nam','var')  %Mask_nam not specified
   Mask_nam = [];
else
    if ~exist(Mask_nam, 'file'), error('Unable to find %s', Mask_nam); end;
end
[p,n] = fileparts(bval_nam);
bvec_nam = fullfile(p, [n '.bvec']);
if ~exist(DWI_nam,'file'), error('Unable to find %s',DWI_nam); end;
if ~exist(bvec_nam,'file'), error('Unable to find %s',bvec_nam); end;
if ~exist(bval_nam,'file'), error('Unable to find %s',bval_nam); end;
%load data
fprintf('%s version 7/2017 loading data\n', mfilename);
%smooth data (optional)
sDWI_nam = DWI_nam;
isSmooth = false;
%if exist('Smooth','var') && Smooth  %default: no smooth - recommended by mrtrix
if ~exist('Smooth','var') || Smooth  %default: smooth - recommended by MUSC
    isSmooth = true;
    [hdr, img] = read_volsSub (DWI_nam);
    fprintf(' Smoothing %d volumes\n', size(img,4));
    [pth,nam] = filepartsSub(DWI_nam);
    sDWI_nam = fullfile(pth, ['s', nam, '.nii']);
    hdr = hdr(1);
    hdr.fname = sDWI_nam;
    for v = 1 : size(img,4) %smooth each volume
        %gaussian kernel with 1.25 voxel FWHM (Neto Henriques et al., 2012a; 2012b)
        %http://academicdepartments.musc.edu/cbi/dki/User%20Manual%202.6.pdf
        % Gaussian smoothing filter. An isotropic FWHM of roughly 1.25 times the voxel size is recommended.
        src = img(:,:,:,v) + 0; %+0 forces new matrix
        dest = src + 0;
        spm_smooth(src,dest,1.25,0);
        %img(:,:,:,v) = dest;
        hdr.n(1)=v;
        spm_write_vol(hdr,dest);
    end
end
fprintf(' Computing mean kurtosis\n');
start = tic;
[dt, dkt] = dwi2tensorSub(sDWI_nam, bvec_nam, bval_nam, Mask_nam);
tensor2metricSub(sDWI_nam, dt);
%fprintf(' %g seconds, Mean kurtosis values range from %g..%g (clipped to %g..%g)\n', toc(start), min(MK(:)), max(MK(:)), min_kurtosis, max_kurtosis );
fprintf('dwi2tensor %g seconds\n', toc(start) );
% function [X,param] = estimate_parameters(DT,KT,'mk')
[hdrT, imgT] = read_volsSub (dt);
[hdrK, imgK] = read_volsSub (dkt);
hdr = hdrT(1);
imgM = any(imgT,4)+any(imgK,4);
imgM(imgM > 0) = 1;
MK_nam = fullfile(p, [n, '_MKmask.nii']);
save_volSub(MK_nam, hdr, imgM);
imgT = reshape(imgT,[],6)';
imgK = reshape(imgK,[],15)';
%(6 x nvox) and KT is kurtosis tensor (15 x nvox)
[imgMK,param] = estimate_parameters(imgT,imgK,{'mk'});
%save MK map
imgMK = reshape(imgMK, hdr.dim(1), hdr.dim(2), hdr.dim(3));

imgMK(((imgM(:) > 0) .* (imgMK(:) == 0)) > 0) = nan; %make all MK = zero values = NaN
imgMK(isnan(imgMK(:))) = nan; %make all air NaN [optional]

MK_nam = fullfile(p, [n, '_MK.nii']);
save_volSub(MK_nam, hdr, imgMK)
if ~isDeleteTempImages, return; end; %next lines delete temp files
delete(dt);
delete(dkt);
if isSmooth, delete(sDWI_nam); end;
%end dkifx2()

% --- local sub-functions follow ---

function tensor2metricSub(sDWI_nam, dt)
if ~exist(dt, 'file'), warning('Unable to find %s', dt); return; end;
exenam = '/usr/local/bin/tensor2metric';
if ~exist(exenam, 'file')
	[~,exenam] = system('which tensor2metric');
	exenam = strip(exenam);
end
if ~exist(exenam, 'file')
   userDir = char(java.lang.System.getProperty('user.home'));
   exenam2 = fullfile(userDir, 'mrtrix3','bin','tensor2metric');
   if exist(exenam2, 'file')
    exenam = exenam2;
   else
    error('unable to find %s', exenam); 
   end
end
[pth,nam,ext] = filepartsSub(sDWI_nam);
fa = fullfile(pth, [nam, '_FAx.nii']);
v1 = fullfile(pth, [nam, '_V1x.nii']);
cmd = sprintf('%s -force -quiet -fa "%s" -vector "%s" "%s"', exenam, fa, v1,  dt);
status = systemSub(cmd);
if status ~= 0, error('unable to run command: %s', cmd); end;
%end dwidenoise()


function [dt, dkt] = dwi2tensorSub(dwi, bvec_nam, bval_nam,  mask)
if ~exist(dwi, 'file'), warning('Unable to find %s', dwi); return; end;
exenam = '/usr/local/bin/dwi2tensor';
if ~exist(exenam, 'file')
	[~,exenam] = system('which dwi2tensor');
	exenam = strip(exenam);
end
if ~exist(exenam, 'file')
   userDir = char(java.lang.System.getProperty('user.home'));
   exenam2 = fullfile(userDir, 'mrtrix3','bin','dwi2tensor');
   if exist(exenam2, 'file')
    exenam = exenam2;
   else
    error('unable to find %s', exenam); 
   end
end
[pth,nam,ext] = filepartsSub(dwi);
dt = fullfile(pth, [nam, '_DT.nii']);
dkt = fullfile(pth, [nam, '_DK.nii']);
cmd = '';
if ~isempty(mask)
    cmd = [' -mask "', mask, '" '];
end
%cmd = sprintf('%s %s -force -fslgrad "%s" "%s" -dkt "%s" "%s" "%s"', exenam, cmd, bvec_nam, bval_nam,  dkt, dwi,  dt);
cmd = sprintf('%s %s -quiet -force -fslgrad "%s" "%s" -dkt "%s" "%s" "%s"', exenam, cmd, bvec_nam, bval_nam,  dkt, dwi,  dt);
%cmd = fslCmdSub (cmd);
status = systemSub(cmd);
if status ~= 0, error('unable to run command: %s', cmd); end;
%end dwi2tensorSub()

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

function status  = systemSub (cmd)
% Save library paths
MatlabPath = getenv('LD_LIBRARY_PATH');
% Make Matlab use system libraries
setenv('LD_LIBRARY_PATH',getenv('PATH'));
fprintf('%s\n',cmd);
status = system(cmd);
% Reassign old library paths
setenv('LD_LIBRARY_PATH',MatlabPath);
%systemSub
