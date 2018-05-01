function nii_famask(fnm, isSaveCopy)
%remove bright voxels on the rim of the cortex (e.g. clean up FA maps)
%  fnm : filename of image to correct
%  isSaveCopy : (optional) if true then original is saved
%Examples
% nii_famask; %use gui
% nii_famask(DTI_M2001_R01_FA.nii)

if ~exist('fnm','var')
	fnm = spm_select(1,'image','Select image[s] for NaN removal'); 
end
%load image
hdr = spm_vol(fnm);
maskKey = 'FAmask'; %we will save this in the header so we do not run twice
if ~isempty(strfind(hdr.descrip, maskKey))
    fprintf('%s skipping image (already masked): %s\n', mfilename, fnm);
    return;
end
img = spm_read_vols(hdr);
if exist('isSaveCopy','var') && isSaveCopy
    [pth, nm, ext] = spm_fileparts(fnm);
    hdro = hdr;
    hdro.fname = fullfile(pth, [ nm, '_premask', ext]);  
    spm_write_vol(hdro,img);
end
img(isnan(img)) = 0;%max(img(:)); % use ~isfinite instead of isnan to replace +/-inf with zero
imgErode = dilateErodeSub(img, false); %1st erosion
imgErode = dilateErodeSub(imgErode, false); %2nd erosion
imgDiff = diffSub(img); %compute edges
imgDiff(img == 0) = 0; %only voxels in brain
thresh = threshSub(img, 3); %compute 3 standard-deviation
imgDiff(imgErode > 0) = 0; %only erode surface voxels
imgDiff(imgDiff < thresh) = 0; %only erode voxels brighter than threshold
img(imgDiff > 0) = 0;
fprintf('%s masked %d voxels from %s\n', mfilename, sum(imgDiff(:) > 0), fnm);
%save image
hdr.descrip = [maskKey,' ', hdr.descrip]; %mark this image as masked
%uncomment next lines if you do not want to overwrite image
% [pth, nm, ext] = spm_fileparts(fnm);
% hdr.fname = fullfile(pth, ['q', nm, ext]);  
spm_write_vol(hdr,img);
%end nii_fingemask()


function thresh = threshSub(img, stdev)
%e.g. stdev 1.96 will report 1.96 standard deviations
im = img(img ~= 0);
im = im(:);
thresh = std(im) * stdev;
%end threshSub()


function sumimg = diffSub(img)
% how different is voxel brightness from neighbors
dims = size(img);
if size(img,3) < 3, error('Not a 3D volume'); end;
img = img(:); %convert 3D to 1D vector
sumimg = zeros(size(img));
sumimg = abs(img - [img(2:end); 0]); %dilate left
sumimg = sumimg + abs(img - [0; img(1:end-1)]); %dilate right
sumimg = sumimg + abs(img - [img(dims(1)+1:end); zeros(dims(1),1)]); %dilate anterior
sumimg = sumimg + abs(img - [zeros(dims(1),1); img(1:end-dims(1))]); %dilate posterior
sliceVox = dims(1) * dims(2); %voxels per slice
sumimg = sumimg + abs(img - [img(sliceVox+1:end); zeros(sliceVox,1)]); %dilate inferior
sumimg = sumimg + abs(img - [zeros(sliceVox,1); img(1:end-sliceVox)]); %dilate superior
sumimg = reshape(sumimg, dims); %convert from 1D to 3D
%end dilateErodeSub();

function imgD = dilateErodeSub(img, doDilate)
dims = size(img);
if size(img,3) < 3, error('Not a 3D volume'); end;
img = img(:); %convert 3D to 1D vector
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
if (~doDilate) %no files
    imgD = (imgD > 6); %binarize ERODED
else
    imgD = (imgD > 0); %binarize DILATED
end
imgD = reshape(imgD, dims); %convert from 1D to 3D
%end dilateErodeSub();