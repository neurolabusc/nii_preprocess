%GY, Oct 5, 2017
function nii_fmri60(fmriname, t1name, t2name, lesion_name)
%analyze data in a naming task
% t1name   : filename for T1 scan
% fmriname : filename for fMRI scan


%if ~exist('t1name','var')
%	t1name = spm_select(1,'image','Select T1 scan'); 
%end
if ~exist('fmriname','var')
	fmriname = spm_select(1,'image','Select fMRI scan'); 
end
if ~exist('t1name','var') || isempty(t1name)
	p.t1name = -1;% t1name; 
else
    p.t1name = t1name;
end
if exist('t2name','var') && ~isempty(t2name)
	p.t2name = t2name; 
    %GY, Oct 5, 2017
    if exist('lesion_name','var') && ~isempty(lesion_name)
        p.lesion_name = lesion_name; 
    end
    
end
p.fmriname = fmriname;
p.setOrigin = true; 
p.TRsec = 10; %repeat time is 1.92 seconds per volume
p.slice_order = -1; %SKIP
p.phase =  ''; %phase image from fieldmap
p.magn = ''; %magnitude image from fieldmap
%statistical information (optional: if not provided not statitics)
p.names{1} = 'naming';
p.names{2} = 'abstract';
%onsets for 1st session, onsets{session, cond}
p.duration{1} = 3; %duration 1 for events, longer for blocks
%condition 1: named pictures
p.onsets{1,1} = [12.728 22.185 36.612 47.654 66.134 84.632 105.48 124.528 134.769 165.825 182.822 197.149 204.871 216.03 246.686 255.876 262.364 274.957 303.245 324.594 333.234 345.543 355.083 363.873 373.581 384.239 397.282 406.905 417.013 436.427 447.486 464.582 473.639 502.594 516.771 537.836 545.776 577.149 583.421 593.595 602.535 623.701 646.217 673.756 697.172 714.769 722.241 735.517 743.874 754.498 787.156 792.994 804.636 817.345 827.353 835.11 844.816 856.224 882.994 902.776 922.273 943.456 952.996 966.607 972.778 985.921 996.913 1017.511 1026.351 1046.016 1053.338 1064.696 1083.11 1096.822 1102.691 1113.833 1124.575 1132.18 1154.097 1166.806];
%condition 2: silent abstract images
p.onsets{1,2} = [7.541 53.809 76.792 92.255 115.021 143.892 154.017 176.667 224.02 233.376 283.163 297.824 316.889 422.767 452.773 487.716 492.688 526.645 552.364 565.857 612.075 637.978 657.893 665.265 685.58 707.18 763.989 776.565 866.182 873.053 897.572 913.1 937.902 1004.836 1037.376 1074.07 1144.122 1174.011 1183.47 1196.762];

p.onsets{1,1} = p.onsets{1,1} - 1.0;
p.onsets{1,2} = p.onsets{1,2} - 1.0;
%p.onsets{1,1} = p.onsets{1,1} - (p.TRsec/2); %we will time from start not middle of acq
%p.onsets{1,2} = p.onsets{1,2} - (p.TRsec/2); %we will time from start not middle of acq
%CR 10/2017 SPARSE DESIGN subtract one second, TA=2.0, no slice timing correction, T0=1
p.onsets{1,1} = p.onsets{1,1} - 1.0;
p.onsets{1,2} = p.onsets{1,2} - 1.0;
p.statAddSimpleEffects = true;
cropVolumesSub(fmriname, 0, 60); %Crop to precisely 60 volumes
p.mocoRegress = false; %should motion parameters be included in statistics?
nii_batch12(p);
%end nii_fmri60()

function cropVolumesSub(fnm, skipVol, nVol)
%ensure that image 'fnm' has no more than nVol volumes
[p,n,x] = spm_fileparts(fnm);
matName = fullfile(p,[n,'.mat']);
if exist(matName,'file'), fprintf('Deleting mat file (motion correction will replace this) : %\n', matName); delete(matName); end;
rpName = fullfile(p,['rp_', n,'.txt']);
if exist(rpName,'file'), fprintf('Deleting text file (motion correction will replace this) : %\n', rpName); delete(rpName); end;
hdr = spm_vol(fnm); 
img = spm_read_vols(hdr);
[nX nY nZ nV] = size(img);
if (skipVol < 0), skipVol = 0; end;
if ( (skipVol+nVol) > nV), nVol = nV - skipVol; end;
if ((skipVol == 0) && (nVol == nV)), return; end;
if ~strcmpi(x,'.nii') %unzip compressed data
    error('Only able to process ".nii" images %s', mfilename);
end;
fnmOrig = fullfile(p, [n, '_v', int2str(nV), x]);
movefile(fnm, fnmOrig);
fprintf('%s has volumes %d..%d volumes from %s\n',fnm,skipVol+1,skipVol+nVol,fnmOrig); 
hdr = hdr(1);
hdr.fname   = fnm;
for vol=1:nVol
    hdr.n(1)=vol;
    imgMod = img(:, :, :, skipVol+vol);
    spm_write_vol(hdr,imgMod(:, :, :, 1));
end;
%end cropVolumesSub();