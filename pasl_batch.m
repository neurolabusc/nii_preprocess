function pasl_batch (ASLnames, T1name, Delaytime,mm);
%preprocesses Siemens PASL data - requires SPM8 and ASLtbx
%  ASLnames  : 4D PASL image(s), should have odd number of volume (M0 followed by  labelled/unlabelled pairs)
%  T1name    : T1-weighted anatomical scan from participant
%  Delaytime : Labelling delay time in SECONDS "inversion time 2"-"inversion time 1"
%  mm        : Reslice to this resolution
% Reference: Ze Wang, Geoffrey K. Aguirre, Hengyi Rao, Jiongjiong Wang,
% Maria A. Fernandez-Seara, Anna R. Childress, and John A. Detre, Epirical
% optimization of ASL data analysis using an ASL data processing toolbox:
% ASLtbx, Magnetic Resonance Imaging, 2008, 26(2):261-9.
%
% for details https://www.cfn.upenn.edu/~zewang/ASLtbx_manual.pdf
%Siemens product PASL (PICORE): post labeling delay time is"inversion time 2"-""inversion time 1"
%Delaytime = 1.1 for LIME protocol "inversion time 1 = 700ms" "inversion time 2 = 1800ms"
%Examples
% pasl_batch('1crPASL1000ms.nii','crT1.nii',1); %single session ASL with T1
       
if exist('asl_perf_subtract')~=2; fprintf('%s requires ASLtbx\n',which(mfilename)); return; end;
if exist('spm_realign_asl')~=2; fprintf('%s requires ASLtbx\n',which(mfilename)); return; end;
if exist('spm')~=2; fprintf('%s requires SPM\n',which(mfilename)); return; end;
if ~exist('Delaytime'), 
    Delaytime=1.1; 
    fprintf('%s warning: assuming label delay time of %.3f seconds.\n',which(mfilename),Delaytime);  
else
    fprintf('Assuming label delay time of %.3f SECONDS.\n',Delaytime);
end;
if ~exist('mm')
	mm = 2;
end
spm('Defaults','fMRI');
spm_jobman('initcfg'); % useful in SPM8 only
clear matlabbatch

if nargin <1 %no input: select ASL file[s]
 ASLnames = spm_select(inf,'image','Select 4D ASL volumes[s]');
end
if nargin <2 %no input: select T1 file[s]
 T1name = spm_select(1,'image','Select T1 anatomical scan');
end
tic; %start timer...
% 0 - ONCE PER INDIVIDUAL normalize T1 scan using unified normalization-segmentation...
[pthT1,namT1,extT1,vol] = spm_fileparts( deblank (T1name));
mT1name = priorSegnormSub(T1name); %use precomputed segmentation if available
if length(mT1name) < 1 %no previous segmentation - this allows prior lesion masked segmentation
	mT1name = segnormSub(T1name);
	segnormwriteSub(T1name,{mT1name}, 2); %warp bias corrected image
	c1name = fullfile(pthT1,['c1', namT1, extT1]); %segmented gray matter
	c2name = fullfile(pthT1,['c2', namT1, extT1]); %segmented white matter
	c3name = fullfile(pthT1,['c3', namT1, extT1]); %segmented CSF
	segnormwriteSub(T1name,{c1name}, mm); %warp gm with normalization parameters
	segnormwriteSub(T1name,{c2name}, mm); %warp wm with normalization parameters
	segnormwriteSub(T1name,{c3name}, mm); %warp csf with normalization parameters
end
c1name = binarizeSub(fullfile(pthT1,['wc1', namT1, extT1]), 0.90,true);%resliced GM
c2name = binarizeSub(fullfile(pthT1,['wc2', namT1, extT1]), 0.95,true);%resliced WM
c3name = binarizeSub(fullfile(pthT1,['wc3', namT1, extT1]), 0.90,true);%resliced CSF
%REMAINING STEPS ONCE PER ASL SCAN
nses = length(ASLnames(:,1));
for ses = 1 : nses
    Filename = deblank (ASLnames(ses,:));
    nVol = numVolumesSub(Filename);
    if nVol < 3 %this script uses 4d data
        fprintf('ERROR: %s requires multiple volumes\n',which(mfilename));
        return;
    end;
    if mod(nVol,2) == 0
        fprintf('%s error: Siemens PASL files should have an odd number of volumes (M0 image followed by ASL pairs)\n',which(mfilename));
        return;
    end;
    M0name = cropVolumesSub (Filename, 0, 1,'m0');
    Filename = cropVolumesSub (Filename, 1, inf); %remove first volume
    Filename = vol4DSub(Filename); %'img.nii'-> 'img.nii,1','img.nii,2'...    
    fprintf('Pre-proocessing %d volumes of session %d\n',length(Filename(:,1)),ses);
    fprintf(' Assuming origin was manually set to Anterior Commissure (SPM''s ''Display'' function).\n');
    % 1 - motion correct ASL images
    spm_realign_asl(Filename); %estimate motion correction - n.b. ASL labelled and unlabelled are different 
    % 2 - reslice motion corrected images (and create mean)
    [meanname, Filename] = applyRealignSub(Filename);
    % 3a - coregister ASL to T1 image
    coregSub(meanname, T1name, Filename);
    % 3b1 - coregister 1st volume (M0 image) to mean ASL...
    coregSub(char(M0name),meanname, M0name);
    % 3b2 - warp segmented T1 images to  ASL space
    M0seg = fullfile(pthT1,['c2', namT1, extT1]); %segmented white matter
    coregWriteSub(mT1name, M0name,  {M0seg});
    M0seg = binarizeSub(fullfile(pthT1,['rc2', namT1, extT1]), 0.95,true);%resliced WM
    
    %c1name = fullfile(pthT1,['c1', namT1, extT1]); %segmented gray matter
    %c2name = fullfile(pthT1,['c2', namT1, extT1]); %segmented white matter
    %c3name = fullfile(pthT1,['c3', namT1, extT1]); %segmented CSF
    %coregWriteSub(mT1name, M0name,  cellstr(strvcat(c1name,c2name,c3name)));
    %c1name = binarizeSub(fullfile(pthT1,['rc1', namT1, extT1]), 0.90,true);%resliced GM
    %c2name = binarizeSub(fullfile(pthT1,['rc2', namT1, extT1]), 0.95,true);%resliced WM
    %c3name = binarizeSub(fullfile(pthT1,['rc3', namT1, extT1]), 0.90,true);%resliced CSF

    % 4 - smooth with 6mm FWHM Gaussian
    Filename = smoothSub(Filename, 6);
    % 5 - conduct CBF estimates...
    fprintf('\nWARNING: next values must be correct - please match Siemens protocol PDF with %s\n',which(mfilename));
    FirstimageType =1; %the first PASL volume after the M0 is labelled
    fprintf('FirstimageType = %d (0=label,1=control)\n',FirstimageType);%0:label; 1:control; for the sequence (PASL and CASL) distributed by CFN, the first image is set to be label.
    SubtractionType = 0;
    fprintf('SubtractionType = %d (0=simple,1=surround,2=sinc)\n',SubtractionType);
    SubtractionOrder = 0; %***** -> 0, 1 for pCASL
    fprintf('SubtractionOrder = %d (0=label-control[FAIR-PASL],1=control-label[CASL])\n',SubtractionType);
    fprintf(' NOTE: if resulting CBF maps have negative values in gray matter, you have set the subtraction order incorrectly.');
    Flag = [1 1 1 0 0 1 1 1 1];
    %Flag = [1MaskFlag,2MeanFlag,3CBFFlag,4BOLDFlag,5OutPerfFlag,6OutCBFFlag,7QuantFlag,8ImgFormatFlag,9D4Flag]
    fprintf('Mask perfusion images? = %d (0=no, 1=yes)\n',Flag(1));
    fprintf('Create mean image? = %d (0=no, 1=yes)\n',Flag(2));
    fprintf('CBF quantifications? = %d (0=no, 1=yes)\n',Flag(3));
    fprintf('Save BOLD images? = %d (0=no, 1=yes)\n',Flag(4));
    fprintf('Save Perf images? = %d (0=no, 1=yes)\n',Flag(5));
    if Flag(3) %CBFFlag
        fprintf('Save CBF images? = %d (0=no, 1=yes)\n',Flag(6));
        fprintf('Using a unique M0 value for all voxels? = %d (0=no, 1=yes)\n',Flag(7));
    end;
    fprintf('Save as NIfTI format (else Analyze)? = %d (0=no, 1=yes)\n',Flag(8));
    fprintf('Save as 4D? = %d (0=no, 1=yes)\n',Flag(9));
    Timeshift = 0.5; %NOT USED: only for sinc interpolation
    AslType = 0; %***** -> 0, 1 for pCASL
    fprintf('AslType = %d (0=PASL,1=CASL/pCASL)\n',AslType);
    labeff = 0.95; %  labeling efficiency, this should be measured for onsite scanner. %***** -> 0.95
    fprintf('labeling efficiency = %g (0.95 for PASL, 0.68 for CASL, 0.85 for pCASL)\n',labeff); 
    MagType = 1;% 0=1.5T, 1=3T
    fprintf('Field Strength = %dT\n',(MagType+1)*1.5);
    Labeltime =  80*0.0185; %The CFN pCASL RF block duration is ALWAYS = 0.0185s (20 RF pulses with gaps) - 18500us 
    fprintf('Num RF Blocks = %g\n',Labeltime/0.0185);  % *********
    
    fprintf('Post Label Delay = %dusec\n',Delaytime*1000000);
    Slicetime = 35.88235294;
    fprintf('Slicetime = %gms\n',Slicetime); % ****
    fprintf(' To compute slicetime get the minimal TR by clicking the "TR" window in the protocol panel, then slicetime=(minTR-labelingtime-delaytime)/#slices.');
    fprintf(' For example, the MCBI default sequence has a minimum TR= 2090ms plus the delay time (e.g. with a 1200ms delay, the minimum time is 3290ms)' ); 
    fprintf(' For example if MinTR=3090, labelingtime=1480, delaytime=1000, slices=17, then slicetime= 35.88235294');
      %CR note this is much longer than the EPI readout time  15.04ms = 0.47ms *32 lines [64 matrix, iPAT=2] 
    TE=11; %*** 12ms for pCASL
    fprintf('Echo Time (TE) = %gms\n',TE); 
    %M0seg = segmentedC2; %PASL only segmented white matter M0 image
    maskimg = ''; %predefined mask image 
    FullFilename = strvcat(FullpathSub('', Filename));
    % next line uses the parameters from the previous lines
    asl_perf_subtract(FullFilename,FirstimageType, SubtractionType,SubtractionOrder,Flag,Timeshift,AslType,labeff,MagType,Labeltime,Delaytime,Slicetime,TE,M0name,M0seg,maskimg);
    meanname = addprefix('meanCBF_0_', strvcat(deblank (Filename(1,:))  ) );
    % 6 - apply normalization parameters to mean cbf image
    segnormwriteSub(T1name,{meanname}, mm);
    meanname = addprefix('w', meanname);
    maskSub(meanname,0.5);
    descriptiveStatisticsSub (meanname, strvcat(c1name,c2name,c3name)); %report mean GM, WM, CSF
end;
fprintf('Done processing %d sessions in %0.3fsec\n',nses, toc);
%%END function pasl_process

%%%%%%% SUBFUNCTIONS FOLLOW

function [sesvols] = vol4DSub1st(sesvol1);
% input: filename from single image in 4D volume, output: first volume in 4D dataset
%  example 4D file with 3 volumes input= 'img.nii', output= {'img.nii,1'}
[pth,nam,ext,vol] = spm_fileparts( deblank (sesvol1));
sesname = fullfile(pth,[nam, ext]);
hdr = spm_vol(sesname);
nvol = length(hdr);
if (nvol < 1), fprintf('vol4DSub1st error: image data required %s\n', sesname); return; end;
sesvols = cellstr([sesname,',1']);
%%END vol4DSub1st()


function coregWriteSub(sourceName, refName, inname)
fprintf('Coregistering %s to %s, and reslicing %d images.\n',sourceName,refName,length(inname(:,1)) );
matlabbatch{1}.spm.spatial.coreg.estwrite.ref = {refName};
matlabbatch{1}.spm.spatial.coreg.estwrite.source ={sourceName};
matlabbatch{1}.spm.spatial.coreg.estwrite.other = inname;
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'nmi';
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.sep = [4 2];
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.fwhm = [7 7];
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.interp = 1;
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.mask = 0;
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.prefix = 'r';
spm_jobman('run',matlabbatch); %make mean motion corrected
%%END coregWriteSub()

function [sesvols] = vol4DSub(sesvol1);
% input: filename from single image in 4D volume, output: list of filenames for all volumes in 4D dataset
%  example 4D file with 3 volumes input= 'img.nii', output= {'img.nii,1';'img.nii,2';'img.nii,3'}
[pth,nam,ext,vol] = spm_fileparts( deblank (sesvol1));
sesname = fullfile(pth,[nam, ext]);
hdr = spm_vol(sesname);
nvol = length(hdr);
if (nvol < 2), fprintf('Error 4D fMRI data required %s\n', sesname); return; end;
sesvols = cellstr([sesname,',1']);
for vol = 2 : nvol
        sesvols = [sesvols; [sesname,',',int2str(vol)]  ];
end;
%END vol4DSub()

function [longname] = addprefix4DSub(prefix, shortname);
% input: filenames from single image in 4D volume, output: list of filenames for all volumes in 4D dataset
%  example 4D file with 3 volumes input= {'img.nii,1';'img.nii,2';'img.nii,3'}, output= {'simg.nii,1';'simg.nii,2';'simg.nii,3'}
nsessions = length(shortname(:,1));
longname = cell(nsessions,1);
for s = 1 : nsessions 
	[pth,nam,ext,vol] = spm_fileparts(strvcat( deblank (shortname(s,:))) );
	sesname = fullfile(pth,[prefix, nam, ext,vol]);
    longname(s,1) = {sesname};
end;
%END addprefix4DSub()

function [longname] = addprefix(prefix, shortname);
%adds path if not specified
[pth,nam,ext,vol] = spm_fileparts(shortname);
longname = fullfile(pth,[prefix nam, ext]);
if exist(longname)~=2; fprintf('Warning: unable to find image %s - cd to approrpiate working directory?\n',longname); end;
longname = [longname, ',1'];
%%END addprefix()

function [meanname, outname] =applyRealignSub(inname)
fprintf('Reslicing %d image with motion correction parameters.\n',length(inname(:,1)));
matlabbatch{1}.spm.spatial.realign.write.data =inname;
matlabbatch{1}.spm.spatial.realign.write.roptions.which = [2 1]; 
matlabbatch{1}.spm.spatial.realign.write.roptions.interp = 4;
matlabbatch{1}.spm.spatial.realign.write.roptions.wrap = [0 0 0];
matlabbatch{1}.spm.spatial.realign.write.roptions.mask = 1;
matlabbatch{1}.spm.spatial.realign.write.roptions.prefix = 'r';
spm_jobman('run',matlabbatch); %make mean motion corrected
meanname = addprefix('mean', strvcat(deblank (inname(1,:))  ) );
outname = addprefix4DSub('r', inname); %add 'r'ealigned prefix
%%END  applyRealignSub()

function [outName] = binarizeSub(imgName, thresh,zeroEdgeSlices)
%creates binary version of imgName where any voxel >= thresh is set to one
[pth,nam,ext,vol] = spm_fileparts( deblank (imgName));
hdr = spm_vol(imgName);
[img] = spm_read_vols(hdr);
imgBin= zeros(size(img));
imgBin((img >= thresh)) = 1; 
outName = fullfile(pth,['b', nam, ext]); %'b'inary
hdr.fname   = outName;
if zeroEdgeSlices %next zero top and bottom 3rd
    z = size(imgBin,3);
    z3rd = floor(z /3);
    if z3rd > 0
        fprintf('Zeroing top and bottom slices of %s\n',outName);
        for sl = 1:z3rd
            imgBin(:,:,sl) = 0; %bottom slices
            imgBin(:,:,z-sl+1) = 0; %top slices
        end;
    end;
end;
spm_write_vol(hdr,imgBin);
%end binarizeSub()

function coregSub(sourcename, refname, inname)
fprintf('Coregistering %s to %s, and applying transforms to %d images.\n',sourcename,refname,length(inname(:,1)) );
matlabbatch{1}.spm.spatial.coreg.estimate.ref = {refname};
matlabbatch{1}.spm.spatial.coreg.estimate.source ={sourcename};
matlabbatch{1}.spm.spatial.coreg.estimate.other = inname;
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];
spm_jobman('run',matlabbatch); %make mean motion corrected
%%END coregSub()

function outname = smoothSub(inname, FWHMmm);
fprintf('Smoothing %d image[s] with a %fmm FWHM Gaussian kernel\n',length(inname(:,1)),FWHMmm);
matlabbatch{1}.spm.spatial.smooth.data = inname;
matlabbatch{1}.spm.spatial.smooth.fwhm = [FWHMmm FWHMmm FWHMmm];
matlabbatch{1}.spm.spatial.smooth.dtype = 0;
matlabbatch{1}.spm.spatial.smooth.im = 0;
matlabbatch{1}.spm.spatial.smooth.prefix = 's';
spm_jobman('run',matlabbatch);
outname = addprefix4DSub('s', inname); %add 's'moothed prefix
%%END smoothSub()

function [bias_t1] = priorSegnormSub(t1);
%returns name of bias file if normalization was already computed
[pth,nam,ext,vol] = spm_fileparts( deblank (t1));
bias_t1 = fullfile(pth,['m', nam, ext]);
if  ~exist(bias_t1, 'file')
    bias_t1 = deblank (t1); 
end
if  ~exist(bias_t1, 'file')
    error('Unable to find T1 image %s', bias_T1);
end
segmat =fullfile(pth,[ nam, '_seg_sn.mat']);
if  (exist(segmat, 'file') == 2)
    fprintf('Skipping normalization (will use existing values from %s)\n',segmat);
    return;
end
bias_t1 = [];
%end priorSegnormSub()

function [bias_t1] = segnormSub(t1);
%estimate normalization based on unified segmentation normalization of T1
fprintf('Unified segmentation normalization of %s\n',t1);
matlabbatch{1}.spm.spatial.preproc.data = {t1};
matlabbatch{1}.spm.spatial.preproc.output.GM = [0 0 1];
matlabbatch{1}.spm.spatial.preproc.output.WM = [0 0 1];
matlabbatch{1}.spm.spatial.preproc.output.CSF = [0 0 1];
matlabbatch{1}.spm.spatial.preproc.output.biascor = 1;
matlabbatch{1}.spm.spatial.preproc.output.cleanup = 0;
matlabbatch{1}.spm.spatial.preproc.opts.tpm = {fullfile(spm('Dir'),'tpm','grey.nii');fullfile(spm('Dir'),'tpm','white.nii');fullfile(spm('Dir'),'tpm','csf.nii')};
matlabbatch{1}.spm.spatial.preproc.opts.ngaus = [2;2;2;4];
matlabbatch{1}.spm.spatial.preproc.opts.regtype = 'mni';
matlabbatch{1}.spm.spatial.preproc.opts.warpreg = 1;
matlabbatch{1}.spm.spatial.preproc.opts.warpco = 25;
matlabbatch{1}.spm.spatial.preproc.opts.biasreg = 0.0001;
matlabbatch{1}.spm.spatial.preproc.opts.biasfwhm = 60;
matlabbatch{1}.spm.spatial.preproc.opts.samp = 3;
matlabbatch{1}.spm.spatial.preproc.opts.msk = {''};
spm_jobman('run',matlabbatch);
%next - return name of bias-corrected image...
[pth,nam,ext,vol] = spm_fileparts( deblank (t1));
bias_t1 = fullfile(pth,['m', nam, ext]);
%%END segnormSub()

function segnormwriteSub(t1,mod, mm);
%reslice ASL/fMRI data based on previous segmentation-normalization
[pth,nam,ext,vol] = spm_fileparts( deblank (t1));
fprintf('Applying unified segmentation normalization parameters from %s to %d image[s], resliced to %fmm\n',t1,length(mod(:,1)),mm);
matlabbatch{1}.spm.spatial.normalise.write.subj.matname = {fullfile(pth,[ nam, '_seg_sn.mat'])};
matlabbatch{1}.spm.spatial.normalise.write.subj.resample = mod;
matlabbatch{1}.spm.spatial.normalise.write.roptions.preserve = 0;
matlabbatch{1}.spm.spatial.normalise.write.roptions.bb = [-78 -112 -50; 78 76 85];
matlabbatch{1}.spm.spatial.normalise.write.roptions.vox = [mm mm mm];
matlabbatch{1}.spm.spatial.normalise.write.roptions.interp = 1;
matlabbatch{1}.spm.spatial.normalise.write.roptions.wrap = [0 0 0];
matlabbatch{1}.spm.spatial.normalise.write.roptions.prefix = 'w';
spm_jobman('run',matlabbatch);
%%END segnormwriteSub()

function outname = maskSub(inname, Thresh);
mask = fullfile(spm('Dir'),'apriori','brainmask.nii');
[pth,nam,ext,vol] = spm_fileparts( inname(1,:));
innameX = fullfile(pth,[nam, ext]); %remove volume label
outname = fullfile(pth,['msk' nam, ext]);
tempname = fullfile(pth,['sk' nam, ext]); 
fprintf('Masking %s with %s at a threshold of %g, resulting in %s\n',innameX,mask,Thresh,outname);
%1 Reslice mask to match image
copyfile(mask,tempname); %move mask - user may not have write permission to SPM folder
matlabbatch{1}.spm.spatial.coreg.write.ref = {innameX};
matlabbatch{1}.spm.spatial.coreg.write.source = {tempname};
matlabbatch{1}.spm.spatial.coreg.write.roptions.interp = 1;
matlabbatch{1}.spm.spatial.coreg.write.roptions.wrap = [0 0 0];
matlabbatch{1}.spm.spatial.coreg.write.roptions.mask = 0;
matlabbatch{1}.spm.spatial.coreg.write.roptions.prefix = 'm';
spm_jobman('run',matlabbatch);
delete(tempname);
%2 Apply mask to image
VImg = spm_vol(innameX);
Img = spm_read_vols(VImg);
VMsk = spm_vol(outname); 
MskImg = spm_read_vols(VMsk);
delete(outname); %<- this was the mask, we will overwrite it with masked image
VImg.fname = outname; %we will overwrite the mask
for x = 1 : size(Img,1)
	for y = 1 : size(Img,2)
		for z = 1 : size(Img,3)
			if MskImg(x,y,z) < Thresh
				Img(x,y,z) = 0;
			end;
		end; %z 
	end; %y 
end; %x 
spm_write_vol(VImg,Img);
%end maskSub()

function [longname] = FullpathSub(prefix, shortname);
% input: appends path to all files (if required)
nsessions = length(shortname(:,1));
longname = cell(nsessions,1);
for s = 1 : nsessions 
	[pth,nam,ext,vol] = spm_fileparts(strvcat( deblank (shortname(s,:))) );
    if length(pth)==0; pth=pwd; end;
	sname = fullfile(pth,[prefix, nam, ext]);
    if exist(sname)~=2; fprintf('Warning: unable to find image %s - cd to approrpiate working directory.\n',sname); end;
	sname = fullfile(pth,[prefix, nam, ext,vol]);
    longname(s,1) = {sname};
end;
%end FullpathSub()

function [nVol] = numVolumesSub (imgName);
%returns number of volumes (timepoints, directions) in a NIfTI image
% imgName : name of source image 
%Examples
% nii_cropVolumes('img.nii'); %returns volumes in img.nii
if ~exist('imgName')
 imgName = spm_select(1,'image','Select 4D image to crop');
end
[pth,nam,ext,vol] = spm_fileparts( deblank (imgName));
imgName = fullfile(pth,[nam, ext]); %remove volume is 'img.nii,1' -> 'img.nii'
hdr = spm_vol([imgName]); 
nVol = length(hdr);
%end numVolumesSub()

function [cropName4D] = cropVolumesSub (imgName4D, skipVol, nVol,prefix);
%given 4D NIfTI image creates NIfTI image with only selected input volumes
% imgName4D : name of source image (4D)
% skipVol   : number of initial volumes deleted from output
% nVol      : number of input volumes
% prefix    : string appended to filename
%Examples
% nii_cropVolumes('img.nii',1,inf); %drop first volume
if ~exist('imgName4D')
 imgName4D = spm_select(1,'image','Select 4D image to crop');
end
if ~exist('prefix')
 prefix = 'c';
end
[pth,nam,ext,vol] = spm_fileparts( deblank (imgName4D));
imgName4D = fullfile(pth,[nam, ext]); %remove volume is 'img.nii,1' -> 'img.nii'
if ~exist(imgName4D), return; end;
if ~exist('skipVol') || ~exist('nVol')
    prompt = {'Skip first N volumes:','Retain N volumes'};
    dlg_title = 'Values for cropping';
    num_lines = 1;
    def = {'0','1'};
    answer = inputdlg(prompt,dlg_title,num_lines,def);
    skipVol = str2num(answer{1});
    nVol = str2num(answer{2});
end;
if (nVol < 1), fprintf('%s quitting: you must retain at least one volume\n',mfilename); return; end;
cropName4D = fullfile(pth,[prefix, nam, ext]); %remove volume is 'img.nii,1' -> 'img.nii'
hdr = spm_vol([imgName4D]); 
[img] = spm_read_vols(hdr);
[nX nY nZ nV] = size(img);
if (skipVol < 0), skipVol = 0; end;
if ( (skipVol+nVol) > nV), nVol = nV - skipVol; end;
if ((skipVol == 0) & (nVol == nV)), fprintf('%s quitting: image only has %d volumes\n',mfilename,nV); return; end;
fprintf('%s has volumes %d..%d volumes from %s\n',cropName4D,skipVol+1,skipVol+nVol,imgName4D); 
hdr = hdr(1);
hdr.fname   = cropName4D;
for vol=1:nVol
    hdr.n(1)=vol;
    imgMod = img(:, :, :, skipVol+vol);
    spm_write_vol(hdr,imgMod(:, :, :, 1));
end;
%edn cropVolumesSub()

function [stats] = descriptiveStatisticsSub (imgName, maskName);
%input: provided with continous image (imgName) and binary mask(s)
%output: mean and stand
% imgName  : continous image
% maskName : binary masking image(s)
%Examples
% nii_descriptiveStatistics('beta.nii','mask1.nii')
% nii_descriptiveStatistics('beta.nii',strvcat('mask1.nii','mask2.nii'))
if ~exist('imgName')
 imgName = spm_select(1,'image','Select continous image');
end
if ~exist('maskName')
 maskName = spm_select(inf,'image','Select binary mask(s)');
end
hdr = spm_vol([imgName]); 
[img] = spm_read_vols(hdr);
nMask = length(maskName(:,1));
stats = zeros(3,nMask);
fprintf('image:\t%s',imgName);
saveText = true;
if saveText
    myfile = fopen('results.txt' ,'at');  
    fprintf(myfile,'image:\t%s',imgName);
end;
for m = 1 : nMask
    msk = deblank (maskName(m,:));
    mhdr = spm_vol([msk]); 
    [mimg] = spm_read_vols(mhdr);
    if (size(img) ~= size(mimg)), fprintf('%s error: dimensions of %s do not match %s\n',mfilename,msk,imgName); end;
    vox = img(mimg ~= 0);
    sd = std(vox(:));
    mn = mean(vox(:));
    stats(1,m) = mn; %mean
    stats(2,m) = sd; %standard deviation
    stats(3,m) = length(vox); %number of voxels in mask
    fprintf('\tmaskName=\t%s\tmaskVoxels=\t%d\tmean=\t%f\tstdev=\t%f',msk,length(vox),mn,sd);
    if saveText
     fprintf(myfile,'\tmaskName=\t%s\tmaskVoxels=\t%d\tmean=\t%f\tstdev=\t%f',msk,length(vox),mn,sd);
    end;   
end;
fprintf('\n');
if saveText, fprintf(myfile,'\n'); fclose(myfile); end;
%end descriptiveStatisticsSub()