function prefix = nii_sliceTime(fMRIname, slice_order, TRsec)
%Adjust for slice time order 
% fMRIname    : NIFTI format 4D fMRI image
% slice_order : Slice order [0..6] use 0 for auto-detect if images
%               from Siemens and converted with dcm2nii from Nov 2013 or later
% TRsec       : Repeat time (seconds), 0 for autodetect
%Examples
% nii_sliceTime('img4d.nii');  %slice timing with automatic slice order and TR detection
% nii_sliceTime('img4d.nii',1); %ascending slice timing with auto TR detection
% nii_sliceTime('img4d.nii',1,2.2); %ascending slice timing with 2200ms TR 
%Siemens have unusual interleaving
% http://cbs.fas.harvard.edu/node/559#slice_order
% https://wiki.cimec.unitn.it/tiki-index.php?page=MRIBOLDfMRI
%these are the possible slice_orders http://nifti.nimh.nih.gov/pub/dist/src/niftilib/nifti1.h
kNIFTI_SLICE_UNKNOWN =  0; %AUTO DETECT
kNIFTI_SLICE_SEQ_INC = 1; %1,2,3,4
kNIFTI_SLICE_SEQ_DEC = 2; %4,3,2,1
%kNIFTI_SLICE_ALT_INC = 3; %1,3,2,4 ascending interleaved
%kNIFTI_SLICE_ALT_DEC = 4; %4,2,3,1 descending interleaved
kNIFTI_SLICE_ALT_INC2 = 5; %2,4,1,3 Siemens interleaved with even number of slices 
kNIFTI_SLICE_ALT_DEC2 = 6; %3,1,4,2 Siemens interleaved descending with even number of slices
if ~exist('fMRIname','var')
    fMRIname = spm_select(inf,'image','Select 4D data for slice timing');
end
if ~exist('slice_order','var')
    slice_order = 0; %<- set your desired slice order, 0 for auto-detect
end
if ~exist('TRsec','var')  
    TRsec = 0; %0 for autodetect
end;
[pth,nam,ext] = spm_fileparts( deblank(fMRIname(1,:)));
fMRIname1 = fullfile(pth,[ nam, ext]); %'img.nii,1' -> 'img.nii'
if slice_order == 0 %attempt to autodetect slice order
    fid = fopen(fMRIname1);
    fseek(fid,122,'bof');
    slice_order = fread(fid,1,'uint8');
    fclose(fid);
    if (slice_order > kNIFTI_SLICE_UNKNOWN) && (slice_order <= kNIFTI_SLICE_ALT_DEC2)
        fprintf('Auto-detected slice order as %d\n',slice_order);
    else
        fprintf('%s error: unable to auto-detect slice order. Please manually specify slice order or use recent versions of dcm2nii.\n',mfilename);
        return;
    end;
end

hdr = spm_vol([fMRIname1 ',1']);
if TRsec == 0
    TRsec = hdr.private.timing.tspace;
    if TRsec == 0
        fprintf('Error: unable to auto-detect repeat time. Please manually specify TR or use recent versions of dcm2nii.\n');
        return;
    end; 
end
nslices = hdr.dim(3);
if nslices <= 1 %automatically detect TR
    fprintf('Fatal Error: image %s does not have multiple slices per volume - slice time correction inappropriate. Please edit m-file named %s\n',fMRIname,which(mfilename));
    return;
end;
if (slice_order == kNIFTI_SLICE_ALT_INC2) || (slice_order == kNIFTI_SLICE_ALT_DEC2) %sequential
    isSiemens = true;
else
    isSiemens = false;
end;
if (slice_order == kNIFTI_SLICE_SEQ_INC) || (slice_order == kNIFTI_SLICE_SEQ_DEC) %sequential
	so = 1:1:nslices;
else % if sequential else Interleaved
	if (mod(nslices,2) == 0) && (isSiemens) %even number of slices, Siemens
		so =[2:2:nslices 1:2:nslices ];
	else
		so =[1:2:nslices 2:2:nslices];
	end
end
if (mod(slice_order,2) == 0) %isDescending
	so = (nslices+1)-so;
end; %isDescending
if TRsec <= 0 %automatically detect TR
    fprintf('%s: TR not specified.\n',mfilename);
    return;
end;
TA = (TRsec/nslices)*(nslices-1);
fprintf('Slice time correction assumptions (from %s)\n',which(mfilename));

fprintf('  Slice order=%d, slices=%d, TR= %0.3fsec, TA= %fsec, referenced to 1st slice.\n', slice_order,nslices, TRsec,TA);
if (TRsec < 1.0) || (TRsec > 5.0) 
    fprintf('  Aborting: strange Repeat Time (TR). Please edit the m-file.\n');
    if  (TRsec > 5.0) 
          fprintf('  Long TR often used with sparse imaging: if this is a sparse design please set the TA manually.\n');
    end; 
    if  (TRsec < 0.5) 
          fprintf('  Short TR may be due to DICOM-to-NIfTI conversion. Perhaps use dcm2nii.\n');
    end; 
    return;
end; %unusual TR

if size(fMRIname,1) == 1
    sesname = fullfile(pth,[nam, ext]);
	hdr = spm_vol(sesname);
	nvol = length(hdr);
	if (nvol < 2), fprintf('Slice time correction requires multiple volumes %s\n', sesname); return; end;
    sesvols = cellstr([sesname,',1']);
    for vol = 2 : nvol
        sesvols = [sesvols; [sesname,',',int2str(vol)]  ]; %#ok<AGROW>
    end;
    matlabbatch{1}.spm.temporal.st.scans = {sesvols}';
else
    matlabbatch{1}.spm.temporal.st.scans = {cellstr(fMRIname)};
end                                
matlabbatch{1}.spm.temporal.st.nslices = nslices;
matlabbatch{1}.spm.temporal.st.tr = TRsec;
matlabbatch{1}.spm.temporal.st.ta = TA;
matlabbatch{1}.spm.temporal.st.so = so;
matlabbatch{1}.spm.temporal.st.refslice = so(1); %set slice order to the first slice http://www.alivelearn.net/?p=1037
%so(1) is first acquired slice, set it as the refernece
%  to see how this works add the line "fprintf('spm_slice_timing slice %d shift %f\n',k, shift amount);' if the 'for k =' loop of spm_slice_timing.m
fprintf('Setting reference slice as %d\n',so(1));
prefix = 'a';
matlabbatch{1}.spm.temporal.st.prefix = prefix;
spm_jobman('run',matlabbatch);
%end subfunction slicetimeSub