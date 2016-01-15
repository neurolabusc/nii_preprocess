function [ret] = nii_temporal_filter(fnm, tr, saveToDisk,doDetrend, saveAs32BitFloat)
%Bandpass temporal filtering
% fnm        : name of 4D image for temporal filtering
%              (alternatively: if fnm is an numeric image, return filtered volume but do not save to disk)
% tr         : time per volume (Repeat Time, TR) in seconds
% saveToDisk : filtered image saved to disk
% doDetrend  : apply LINEAR detrending (use nii_detrend for nonlinear)
%
%Output : depends on saveToDisk
%   SaveToDisk == True : returns prefix of output volume
%   SaveToDisk == False : returns filtered dataset 
% Chris Rorden (2013), inspired by http://restfmri.net/
%Example: 
%  nii_temporalfilter('4d.nii'); 
%  nii_temporalfilter('4d.nii', 2.2); %TR=2.2sec
%  nii_temporalfilter('4d.nii', 2.2, false); %TR=2.2sec, don't save to disk
%  nii_temporalfilter(imgXYZT, 2); %process 4D data 'XYZT'
    
if ~exist('fnm','var'), fnm = spm_select(1,'image','Select 4D image for filtering'); end
if (~exist('tr','var')) || (tr==0), tr = 1; fprintf('%s assuming TR (sample period) is %f seconds.\n',mfilename,tr); end;
if ~exist('saveToDisk','var'), saveToDisk= true; end;
if ~exist('doDetrend','var'), doDetrend= true; end;
if ~exist('saveAs32BitFloat','var'), saveAs32BitFloat = false; end;
HighPass_LowCutoff=0.01;%Hz; 
LowPass_HighCutoff =0.1;%Hz 
fprintf('%s TR=%.3f seconds, band pass filtering higher than %.4f Hz and lower than %.4f Hz.\n',mfilename,tr,HighPass_LowCutoff,LowPass_HighCutoff);
if isnumeric(fnm)
    filtImg = fnm;
    saveToDisk = false;
else
	[pth nam ext vol] = spm_fileparts(fnm);
	fnm = fullfile(pth, [nam ext]); %remove volume index 'vol' 
	if (exist(fnm) == 0); fprintf('%s unable to find image %s\n',which(mfilename),fnm);  return; end;
	hdr = spm_vol(fnm); %load 4D dataset only once!
	filtImg = spm_read_vols(hdr);
end
[nX nY nZ nV] = size(filtImg);
if (nV < 2); fprintf('%s requires 4D volumes, %s is a 3D volume.\n',which(mfilename),fnm);  return; end;
%subtract mean to baseline correct to zero
meanImg = mean(filtImg,4);
for v = 1: nV
    filtImg(:,:,:,v) = filtImg(:,:,:,v) - meanImg;
end
%apply simple linear detrending
if (doDetrend)
    fprintf('Removing linear trends (use nii_detrend for nonlinear effects)\n');
    filtImg = reshape(filtImg,nX*nY*nZ,nV);
    filtImg = detrend(filtImg');
    filtImg = reshape(filtImg', nX, nY, nZ, nV);
else
    fprintf('%s assumes data already detrended (nii_detrend)\n',mfilename);
end;
%next apply temporal filter
filtImg =permute(filtImg,[4 1 2 3]); %make time first dimension
rt = tr;
%filter=[0.01, 0.1]; %HighPass_LowCutoff=0.01Hz; LowPass_HighCutoff =0.1Hz 
filtImg=fft(filtImg,[],1);
f=(0:size(filtImg,1)-1);
f=min(f,size(filtImg,1)-f);
idx=find(f<HighPass_LowCutoff*(rt*size(filtImg,1))|f>LowPass_HighCutoff*(rt*size(filtImg,1)));
idx=idx(idx>1);
filtImg(idx,:)=0;
filtImg=real(ifft(filtImg,[],1));
filtImg=permute(filtImg,[2 3 4 1]); %make time last dimension   
%next: add mean image
for v = 1: nV
    filtImg(:,:,:,v) = filtImg(:,:,:,v) + meanImg;
end
if (~saveToDisk)
    ret = filtImg;
    return;
end;
%save to disk
ret = 'f'; %prefix for new image
hdr = hdr(1);%spm_vol([fnm ',1' ]); %load 4D dataset only once!
hdr.fname   = fullfile(pth,[ret nam ext]);;
hdr.private.timing.toffset= 0;
hdr.private.timing.tspace= tr;
%next: save data as 32-bit real
if saveAs32BitFloat 
    hdr.dt(1) = 16; %make 32 bit real
    hdr.private.dat.dtype = 'FLOAT32-LE';
    hdr.private.dat.scl_slope = 1;
    hdr.private.dat.scl_inter = 0;
    hdr.pinfo = [1;0;352]; %slope=1, intercept = 0
elseif (spm_type(hdr.dt(1),'intt')) %integer input
    dmx  = spm_type(hdr.dt,'maxval');
    dmn  = spm_type(hdr.dt,'minval');
    mx = max(filtImg(:));
    mn = min(filtImg(:));
    if dmn < 0
    	sf = max(mx/dmx,-mn/dmn);
    else
    	sf = mx/dmx;
    end
    hdr.pinfo = [sf;0;0];
    %filtImg = round(filtImg);
end
for vol=1:nV
    hdr.n(1)=vol;
    spm_write_vol(hdr,filtImg(:, :, :, vol));
end;
%end nii_temporalfilter()
