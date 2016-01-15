function [hdr, zimg] = nii_i3m (Imgname, Maskname,Thresh,SmoothFWHM, maskOutput,maskType, resliceMM)
%Create image where output has Z-score of brightness observed in mask
% Imgname      : Image to Z-transform
% Maskname     : Masking image (brain regions used for mean/stdev)
% Thresh       : Threshold for (normalized) masking image from 0..1
% SmoothFWHM   : Gaussian blur with Full-Width Half-Maximum width (mm)
% maskOutput   : (0..1) Threshold for masking output (0 for none)
% maskType     : 0 = use both hemispheres of mask, 
%                1=brighter hemisphere (e.g. on T1, injury TENDS dark)
%                2=darker hemisphere (on T2, injury is brighter)
%resliceMM     : (optional) rescale image to desired resolution
%Examples
% nii_i3m('T2_LM1001.nii','scmask.nii',0.5,10,0.25,2); %T2 scan
% nii_i3m('T1_LM1001.nii','scmask.nii',0.5,10,0.25,1); %T1 scan
% nii_i3m('wT1_P192.nii','',0.5,10,0.25,1); %i3m T1 image
% [hdr,img] = nii_i3m('T1_LM1001.nii','scmask.nii',0.5,10,0.25,1,3);

if ~exist('Imgname', 'var') %no files
 Imgname = spm_select(1,'image','Select image to transform');
end;
if ~exist('Maskname', 'var') || isempty(Maskname) %no files
    %Maskname = spm_select(1,'image','brain mask');
    p = fileparts(which(mfilename));
    Maskname = fullfile (p,'mask.nii');
    if ~exist(Maskname, 'file')
        nii_makemask; %generate a mask image
    end
end;
if ~exist('Thresh','var') %threshold not specified
 Thresh = 0.5;
end;
if ~exist('SmoothFWHM', 'var') %no files
 SmoothFWHM = 10;
end;
if ~exist('maskOutput', 'var') %no files
 maskOutput = 0;%none
end;
if ~exist('maskType', 'var')
    maskType = 0;
end
if ~exist('resliceMM', 'var')
    resliceMM = 0;
end
if ~exist(Maskname, 'file')
    fprintf('No images created (unable to find image named %s\n',Maskname);
    return;
end
Maskname = which(Maskname); %incase image not in CWD,SPM wants a full path
%load image
hdr = spm_vol(Imgname); %load img
img = spm_read_vols(hdr);
img(isnan(img)) = 0; % use ~isfinite instead of isnan to replace +/-inf with zero
[pth,nam,ext]=fileparts(hdr.fname);
hdr.fname = fullfile(pth,['z',  nam, ext]);
hdr.private.dat.scl_slope = 1;
hdr.private.dat.scl_inter = 0;
hdr.dt(1)=16; %save as 32-bit float uint8=2; int16=4; int32=8; float32=16; float64=64
mhdr = spm_vol(Maskname); %load mask
mask = spm_read_vols(mhdr);
if ~isequal(size(img), size(mask))
    %fprintf('Error: dimension mismatch %s %s\n',Imgname,Maskname);
    %return;
    %[mhdr, mask] = nii_reslice_target(mhdr, mask, hdr);
    fprintf('Note: reslicing image to match mask\n');
    [hdr, img] = nii_reslice_target(hdr, img, mhdr);
    hdr.fname = fullfile(pth,['z',  nam, ext]);
end
mask(isnan(mask)) = 0; % use ~isfinite instead of isnan to replace +/-inf with zero
mask = (mask - min(mask(:))) / ( max(mask(:)) - min(mask(:)) ); %normalize
maskbin = mask(:) > Thresh; %1D
maskbin = reshape(maskbin,size(mask)); %convert back to 3D
if (maskType > 0)
    for maskSide = 1:2 %check intensity of each hemisphere
        maskHemi = maskbin;    
        if (maskSide ~= 1) %only use right, make left = 0 
            Xlo = 1;
            Xhi = round(mhdr.dim(1)/2);
        else %only use left, make right = 0
            Xlo = round(mhdr.dim(1)/2);
            Xhi = mhdr.dim(1);
        end
        for Zi=1:mhdr.dim(3),
            for Yi=1:mhdr.dim(2),
                for Xi=Xlo:Xhi,
                    maskHemi(Xi,Yi,Zi) = 0;
                end; %for dim(1) = X
            end; %for dim(2) = Y
        end; %for dim(3) = Z
        imgmasked = img(maskHemi);
        mnSide(maskSide) = mean(imgmasked); %#ok<AGROW>
        stSide(maskSide) = std(imgmasked); %#ok<AGROW>
    end;%for left and right side
    selectSide = 1;
    if (maskType == 1) && (mnSide(2) > mnSide(1)) %maskType1=select brighter hemisphere
    	selectSide = 2;
    elseif  (maskType == 2) && (mnSide(2) < mnSide(1))  %maskType1=select brighter hemisphere
        selectSide = 2;      
    end
    mn=mnSide(selectSide);
    st=stSide(selectSide);
    fprintf('side 1 mean %f, side 2 mean %f, using side %d for mask\n',mnSide(1),mnSide(2),selectSide); 
 else %mask both hemispheres
    imgmasked = img(maskbin);
    mn=mean(imgmasked);
    st= std(imgmasked);
end %if mask one hemisphere
fprintf('%s has %d voxels, the %d defined by the mask %s have a mean intensity of %f\n',Imgname,length(img(:)), length(imgmasked(:)), Maskname,mn); 
%compute mean and stdev
if (st == 0) 
	fprintf('clinical_zintensity error: can not compute Z-scores when standard error is zero!');
	return;
end;
% next part saves transformed data
zimg=(img-mn)./ st;
if SmoothFWHM > 0 %
	fprintf('Smoothing by %.1fmm\n',SmoothFWHM);
    VOX = sqrt(sum(hdr.mat(1:3,1:3).^2)); %since we are passing a 3D array we need to define smooth in voxels not mm
    s  = SmoothFWHM./VOX;
    origz = zimg;
    spm_smooth(origz,zimg,s,16);
end
if resliceMM > 0 %reslice prior to binary mask
    if resliceMM == 3
        bb = [-78 -112 -50; 78 74 85]; %match bounding box of SPM
    else
        bb = [-78 -112 -50; 78 76 84]; %typical bounding box
    end   
    resliceMM = [resliceMM resliceMM resliceMM];
 
    [hdr,zimg] = nii_resliceSub(hdr,zimg,resliceMM,bb); 
    [mhdr,mask] = nii_resliceSub(mhdr,mask,resliceMM,bb); %#ok<ASGLU> %save to disk 
end
if maskOutput > 0 %length(PreserveName) > 0
	maskbin = mask(:) > maskOutput; %1D
    maskbin = reshape(maskbin,size(mask)); %convert back to 3D
    zimg((maskbin==0)) = NaN;%0; %<- only masked regions....
end;
%save data
if nargout > 1
    return; %save to memory, not disk
end
spm_write_vol(hdr,zimg);
%end nii_i3m()

function [outhdr, outimg] = nii_resliceSub(inhdr,inimg, Voxdim, BB)
% resample images to have specified voxel dims and bounding box
%   works in memory, no data saved to disk
%
%John Ashburner's reorient.m adapted by Ged Ridgway and Chris Rorden
% This version doesn't check spm_flip_analyze_images -- the handedness of
% the output image and matrix should match those of the input.
%Examples
%  nii_reslice('zT1_LM1001.nii','',[1 1 1],[-78 -112 -50; 78 76 84]); %save to disk
%  nii_reslice('zT1_LM1001.nii','',[2 2 2],[-78 -112 -50; 78 76 84]); %save to disk
%  nii_reslice('zT1_LM1001.nii','',[3 3 3],[-78 -112 -50; 78 74 85]); %save to disk note 74mm!!!!
% Next example shows how to use this without writing to disk
%  inhdr = spm_vol(zT1_LM1001.nii); %load header
%  inimg = spm_read_vols(inhdr); %load volume
%  [outhdr,outimg] = nii_reslice(inhdr,inimg, Voxdim, BB); %resize to memory
%  spm_write_vol(outhdr,outimg); %save resized image

if ~isstruct(inhdr)
    inhdr = spm_vol(inhdr); %load img
	inimg = spm_read_vols(inhdr);
end
voxdim = Voxdim;
bb = BB;
% default voxdim to current volume's voxdim, (from mat parameters)
if any(isnan(voxdim))
    vprm = spm_imatrix(inhdr.mat);
    vvoxdim = vprm(7:9);
    voxdim(isnan(voxdim)) = vvoxdim(isnan(voxdim));
end
voxdim = voxdim(:)';
mn = bb(1,:);
mx = bb(2,:);
% default BB to current volume's
if any(isnan(bb(:)))
    vbb = world_bb(inhdr);
    vmn = vbb(1,:);
    vmx = vbb(2,:);
    mn(isnan(mn)) = vmn(isnan(mn));
    mx(isnan(mx)) = vmx(isnan(mx));
end
% voxel [1 1 1] of output should map to BB mn
% (the combination of matrices below first maps [1 1 1] to [0 0 0])
mat = spm_matrix([mn 0 0 0 voxdim])*spm_matrix([-1 -1 -1]);
% voxel-coords of BB mx gives number of voxels required
% (round up if more than a tenth of a voxel over)
imgdim = ceil(mat \ [mx 1]' - 0.1)';
% output image
[pth,nam,ext] = fileparts(inhdr.fname);
outhdr            = inhdr;
outhdr.fname      = fullfile(pth,['r' nam ext]);
outhdr.dim(1:3)   = imgdim(1:3);
outhdr.mat        = mat;
%outhdr = spm_create_vol(outhdr);
outimg = zeros(outhdr.dim(1:3));
for i = 1:imgdim(3)
    M = inv(spm_matrix([0 0 -i])*inv(outhdr.mat)*inhdr.mat); %#ok<MINV>
    outimg(:,:,i) = spm_slice_vol(inimg, M, imgdim(1:2), 1); % (linear interp)
end
if nargout < 2
    fprintf('saving resliced image %s\n',outhdr.fname);
    spm_write_vol(outhdr,outimg);
end
%end nii_reslice()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function bb = world_bb(V)
%  world-bb -- get bounding box in world (mm) coordinates
d = V.dim(1:3);
% corners in voxel-space
c = [ 1    1    1    1
    1    1    d(3) 1
    1    d(2) 1    1
    1    d(2) d(3) 1
    d(1) 1    1    1
    d(1) 1    d(3) 1
    d(1) d(2) 1    1
    d(1) d(2) d(3) 1 ]';
% corners in world-space
tc = V.mat(1:3,1:4)*c;
% bounding box (world) min and max
mn = min(tc,[],2)';
mx = max(tc,[],2)';
bb = [mn; mx];
%end world_bb()