function [cbfName, c1L, c1R, c2L, c2R] = nii_pasl12 (ASLnames, T1name)
%preprocesses Siemens PASL data - requires SPM12 and ASLtbx
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
% nii_pasl12('/Users/rorden/Desktop/pre/LM1010/ASL_LM1010.nii', '/Users/rorden/Desktop/pre/LM1010/T1_LM1010.nii', 1.1, [2 2 2]);
% nii_pasl12('ASL_P022.nii','T1_P022.nii')
% [cbfName, c1L, c1R, c2L, c2R] = nii_pasl12('ASL_P022.nii','T1_P022.nii')
cbfName = ''; c1L =''; c1R = ''; c2L = ''; c2R = '';
maskLeft = true;
isSPM12orNewerSub;
addpath(fullfile(spm('Dir'),'compat')); %ASLtbx uses spm_chi2_plot 
if exist('asl_perf_subtract','file')~=2; fprintf('%s requires ASLtbx\n',which(mfilename)); return; end;
if exist('spm_realign_asl','file')~=2; fprintf('%s requires ASLtbx\n',which(mfilename)); return; end;
% if ~exist('Delaytime','var'), 
%     Delaytime=1.1; 
%     fprintf('%s warning: assuming label delay time of %.3f seconds.\n',which(mfilename),Delaytime);  
% else
%     fprintf('Assuming label delay time of %.3f SECONDS.\n',Delaytime);
% end;
%if ~exist('mm','var')
mm = [2 2 2];%[1.5 1.5 1.5];
%end
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
[T1nameBet, T1name] = priorSegnormSub(T1name); %bT1.nii, y_T1.nii 
if isempty(T1name), return; end;
%REMAINING STEPS ONCE PER ASL SESSION
nses = length(ASLnames(:,1));
for ses = 1 : nses
    Filename = deblank (ASLnames(ses,:));
    setOriginSub(Filename, 3);
    [nVol, nSlices] = numVolumesSub(Filename);
    if nVol < 5 %this script uses 4d data
        fprintf('ERROR: %s requires multiple volumes %s\n',which(mfilename), Filename);
        return;
    end;
    if mod(nVol,2) == 0
        %fprintf('%s error: Sicemens PASL files should have an odd number of volumes (M0 image followed by ASL pairs)\n',which(mfilename));
        %return;
        isPASL = false;
    else
        isPASL = true;
    end;
    M0name = cropVolumesSub (Filename, 0, 1,'m0');
    Filename = cropVolumesSub (Filename, 1, inf); %remove first volume
    Filename = vol4DSub(Filename); %'img.nii'-> 'img.nii,1','img.nii,2'...    
    fprintf('Pre-processing %d volumes of session %d\n',length(Filename(:,1)),ses);
    % 1 - motion correct ASL images
    spm_realign_asl(Filename); %estimate motion correction - n.b. ASL labelled and unlabelled are different 
    % 2 - reslice motion corrected images (and create mean)
    [meanname, Filename] = applyRealignSub(Filename);
    % 3b1 - coregister 1st volume (M0 image) to mean ASL...
    coregSub(char(M0name),T1nameBet, [{meanname}; Filename]);
    % 3b2 - warp segmented T1 images to  ASL space
    [pthT1,namT1, extT1] = spm_fileparts( deblank (T1name));
    %c1
    c1Lm = zeroHemisphereSub (addprefix ('c1', T1name), true, meanname);
    c1Rm = zeroHemisphereSub (addprefix ('c1', T1name), false, meanname);
    c2Lm = zeroHemisphereSub (addprefix ('c2', T1name), true, meanname);
    c2Rm = zeroHemisphereSub (addprefix ('c2', T1name), false, meanname);
    
    if maskLeft
       M0seg = c2Lm;
    else
        M0seg = coregWriteSub(addprefix ('c2', T1name), meanname);
    end
    M0seg = binarizeSub(M0seg, 0.95,true);%resliced WM
    % 4 - smooth with 6mm FWHM Gaussian
    Filename = smoothSub(Filename, 6);
    % 5 - conduct CBF estimates...
    fprintf('\nWARNING: next values must be correct - please match Siemens protocol PDF with %s\n',which(mfilename));
    fprintf('White matter mask is %s\n', M0seg);
    if ~isPASL %pCASL
        
        fprintf('USING pCASL DEFAULTS\n');
        
        labeff = 0.85; %  labeling efficiency, this should be measured for onsite scanner. %***** -> 0.95
        AslType = 1; %***** -> 0, 1 for pCASL
        fprintf(' To compute slicetime get the minimal TR by clicking the "TR" window in the protocol panel, then slicetime=(minTR-labelingtime-delaytime)/#slices.');
        fprintf(' For example, the MCBI default sequence has 17 slices, minimum TR= 2090ms, 80 rf blocks, plus the delay time (e.g. with a 1200ms delay, the minimum time is 3290ms)' ); 
        fprintf(' labelingtime=NumRFBlocks * 18.5');
        fprintf(' For example if MinTR=3090, labelingtime=1480, delaytime=1000, slices=17, then slicetime= 35.88235294');

        if nSlices == 17 %MCBI protocol
            fprintf('Assuming MCBI defaults\n');
            %default MCBI uses 1200ms or 1000ms delay... 17 slices
            % note you get the same slice time regardless of delay (adds to minTRms) 
            MinTRms = 3290; %minTR in milliseconds, includes PostLabelDelay
            Delaytime = 1.2; %in seconds
            RFBlocks = 80;
            TE=12; %Echo time, in ms
        elseif nSlices == 16
            fprintf('Assuming migraine study defaults');
            %Souviks study uses 82 RF blocks...typically 16 slices
            MinTRms = 2940; %minTR in milliseconds, includes PostLabelDelay
            Delaytime = 1; %PostLabelDelay seconds (not microseconds)
            RFBlocks = 82;
            TE=6.7; %Echo time, in ms
        else
            error('Unknown pCASL sequence (%d slices): %s\n', nSlices, Filename{1});
        end
        Labeltime =  RFBlocks*0.0185; %in seconds, The CFN pCASL RF block duration is ALWAYS = 0.0185s (20 RF pulses with gaps) - 18500us 
        Slicetime = (MinTRms - (Delaytime *1000) - (Labeltime * 1000) ) / nSlices; %Slicetime in ms, not sec!!!
    else
        fprintf('Using PASL settings\n');
        %FirstimageType =1; %the first PASL volume after the M0 is labelled
        Delaytime = 1.8;
        AslType = 0; %***** -> 0, 1 for pCASL
        labeff = 0.95; %  labeling efficiency, this should be measured for onsite scanner. %***** -> 0.95
        Slicetime = 35.88235294;
        Labeltime =  80*0.0185; %Not used for PASL
        TE=11; %*** 12ms for pCASL
    end
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

    fprintf('AslType = %d (0=PASL,1=CASL/pCASL)\n',AslType);
    fprintf('labeling efficiency = %g (0.95 for PASL, 0.68 for CASL, 0.85 for pCASL)\n',labeff); 
    MagType = 1;% 0=1.5T, 1=3T
    fprintf('Field Strength = %dT\n',(MagType+1)*1.5);
    fprintf('Post Label Delay = %dusec\n',Delaytime*1000000);
    
    fprintf('Slicetime = %gms\n',Slicetime); % ****
    fprintf(' To compute slicetime get the minimal TR by clicking the "TR" window in the protocol panel, then slicetime=(minTR-labelingtime-delaytime)/#slices.\n');
    fprintf(' For example, the MCBI default sequence has a minimum TR= 2090ms plus the delay time (e.g. with a 1200ms delay, the minimum time is 3290ms)\n' ); 
    fprintf(' For example if MinTR=3090, labelingtime=1480, delaytime=1000, slices=17, then slicetime= 35.88235294\n');
      %CR note this is much longer than the EPI readout time  15.04ms = 0.47ms *32 lines [64 matrix, iPAT=2] 
    
    fprintf('Echo Time (TE) = %gms\n',TE); 
    %M0seg = segmentedC2; %PASL only segmented white matter M0 image
    maskimg = ''; %predefined mask image 
    FullFilename = strvcat(FullpathSub('', Filename));
    % next line uses the parameters from the previous lines
    %M0name
    %FullFilename
    for FirstimageType = 1:-1:0  %):1
        fprintf('FirstimageType = %d (0=label,1=control)\n',FirstimageType);%0:label; 1:control; for the sequence (PASL and CASL) distributed by CFN, the first image is set to be label.
        asl_perf_subtract(FullFilename,FirstimageType, SubtractionType,SubtractionOrder,Flag,Timeshift,AslType,labeff,MagType,Labeltime,Delaytime,Slicetime,TE,M0name,M0seg,maskimg);
        [pthasl,namasl, extasl] = spm_fileparts( deblank (FullFilename(1,:)));
        meanname=fullfile(pthasl, ['meanCBF_' num2str(SubtractionType) '_' namasl  extasl]);
        if ~exist(meanname, 'file') %new versions of ASL toolbox create a different name
            meanname=fullfile(pthasl, ['meanCBF_' num2str(SubtractionType) '_' namasl 'M0CSF' extasl]);
        end
        c1R = maskMeanSub(meanname, c1Lm);%gray matter right used masked left
        c1L = maskMeanSub(meanname, c1Rm);%gray matter left used masked right
        c2R = maskMeanSub(meanname, c2Lm);%white matter right used masked left
        c2L = maskMeanSub(meanname, c2Rm);%white matter left used masked right
        if ((c1R+c1L) < (c2R+c2L)) 
            fprintf('CBF values do not make sense - rerunning with reverse tag/label order %s\n', char(M0name) );
            %fid = fopen('errors.txt','a');
            %fprintf(fid, 'ASL CBF higher in white matter\t%s\n', matName);
            %fclose(fid);
        else
            break;
        end
    end
    
    % 6 - apply normalization parameters to mean cbf image
    meanname = newSegWriteSub(T1name, meanname, mm);
    cbfName = maskSub(meanname,0.5);
    %descriptiveStatisticsSub (meanname, strvcat(c1name,c2name,c3name)); %report mean GM, WM, CSF
end;
fprintf('Done processing %d sessions in %0.3fsec\n',nses, toc);
%end fnii_pasl12()

%hdr = spm_vol([fMRIname1 ',1']);
%if TRsec == 0
%    TRsec = hdr.private.timing.tspace;

function  targetname = newSegWriteSub(t1name, targetname, vox, bb)
%reslice img using pre-existing new-segmentation deformation field
if isempty(targetname) || isempty(t1name), return; end;
if ~exist('bb','var'), bb = [-78 -112 -50; 78 76 85]; end;
if ~exist('vox','var'), vox =[2 2 2]; end;
[pth,nam,ext, vol] = spm_fileparts(t1name); %#ok<NASGU>
defname = fullfile(pth,['y_' nam ext]);
if ~exist(defname,'file')
    error('Unable to find segment deformation image %s',defname);
end
fprintf('Warping %s based on segment of %s\n', targetname, t1name);
matlabbatch{1}.spm.spatial.normalise.write.subj.def = {defname};
matlabbatch{1}.spm.spatial.normalise.write.subj.resample = {targetname};
matlabbatch{1}.spm.spatial.normalise.write.woptions.bb = bb;
matlabbatch{1}.spm.spatial.normalise.write.woptions.vox = vox;
matlabbatch{1}.spm.spatial.normalise.write.woptions.interp = 4;
spm_jobman('run',matlabbatch);
targetname = addprefix('w', targetname); 
%end newSegWriteSub()

function [T1nameBet, T1name] = priorSegnormSub(T1name) %use precomputed segmentation 
[pth,nam,ext,vol] = spm_fileparts( deblank (T1name));
T1nameBet = fullfile(pth,['b', nam, ext]);
T1nameWarp = fullfile(pth,['y_', nam, ext]);
if  ~exist(T1nameWarp, 'file')
    T1nameWarp = fullfile(pth,['y_e', nam, ext]);
    T1name = fullfile(pth,['e', nam, ext]);
end
if  ~exist(T1name, 'file') ||~exist(T1nameWarp, 'file') || ~exist(T1nameBet, 'file')
    fprintf('Please normalize T1 first (could not find %s or %s)\n', T1nameWarp, T1nameBet);
    T1name = '';
end
%end priorSegnormSub()

function sourceName = coregWriteSub(sourceName, refName)
fprintf('Reslicing %s to match %s.\n',sourceName,refName);
matlabbatch{1}.spm.spatial.coreg.write.ref = {refName};
matlabbatch{1}.spm.spatial.coreg.write.source = {sourceName};
matlabbatch{1}.spm.spatial.coreg.write.roptions.interp = 4;
matlabbatch{1}.spm.spatial.coreg.write.roptions.wrap = [0 0 0];
matlabbatch{1}.spm.spatial.coreg.write.roptions.mask = 0;
matlabbatch{1}.spm.spatial.coreg.write.roptions.prefix = 'r';
spm_jobman('run',matlabbatch); %make mean motion corrected
sourceName = addprefix('r',sourceName);%coregistered
%end coregWriteSub()

function coregSub(sourcename, refname, inname)
%core('m0ASL_LM1010.nii','bT1_LM1010.nii', strvcat('cASL_LM1010.nii','meancASL_LM1010.nii'))
% core('m0ASL_LM1010.nii','bT1_LM1010.nii', {'cASL_LM1010.nii','meancASL_LM1010.nii'})
matlabbatch{1}.spm.spatial.coreg.estimate.ref = {refname};
matlabbatch{1}.spm.spatial.coreg.estimate.source ={sourcename};
if exist('inname','var') && ~isempty(inname)
    if ischar(inname), inname = cellstr(inname); end;
    fprintf('Coregistering %s to %s, and applying transforms to %d images.\n',sourcename,refname,numel(inname) );
    matlabbatch{1}.spm.spatial.coreg.estimate.other = inname;
else
     fprintf('Coregistering %s to %s\n',sourcename,refname);
    matlabbatch{1}.spm.spatial.coreg.estimate.other = {''};
end
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];
spm_jobman('run',matlabbatch); %make mean motion corrected
%end coregSub()

function [sesvols] = vol4DSub(sesvol1)
% input: filename from single image in 4D volume, output: list of filenames for all volumes in 4D dataset
%  example 4D file with 3 volumes input= 'img.nii', output= {'img.nii,1';'img.nii,2';'img.nii,3'}
[pth,nam,ext,vol] = spm_fileparts( deblank (sesvol1));
sesname = fullfile(pth,[nam, ext]);
hdr = spm_vol(sesname);
nvol = numel(hdr);
if (nvol < 2), fprintf('Error 4D fMRI data required %s\n', sesname); return; end;
sesvols = cellstr([sesname,',1']);
for vol = 2 : nvol
        sesvols = [sesvols; [sesname,',',int2str(vol)]  ];
end;
%end vol4DSub()

function [longname] = addprefix4DSub(prefix, shortname)
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

function [longname] = addprefix(prefix, shortname)
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

function outname = maskSub(inname, Thresh);
pth = fileparts(which(mfilename));
mask = fullfile(pth,'mask.nii');
if ~exist(mask,'file'), error('unable to find mask %s',mask); end;
%mask = fullfile(spm('Dir'),'apriori','brainmask.nii');
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
hdr = spm_vol(innameX);
img = spm_read_vols(hdr);
VMsk = spm_vol(outname); 
MskImg = spm_read_vols(VMsk);
delete(outname); %<- this was the mask, we will overwrite it with masked image
hdr.fname = outname; %we will overwrite the mask
img(MskImg < Thresh) = nan;
img(img == 0) = nan;
% for x = 1 : size(Img,1)
% 	for y = 1 : size(Img,2)
% 		for z = 1 : size(Img,3)
% 			if MskImg(x,y,z) < Thresh
% 				Img(x,y,z) = 0;
% 			end;
% 		end; %z 
% 	end; %y 
% end; %x 
spm_write_vol(hdr,img);
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

function [nVol, nSlices] = numVolumesSub (imgName);
%returns number of volumes (timepoints, directions) in a NIfTI image
% imgName : name of source image 
%Examples
% nii_cropVolumes('img.nii'); %returns volumes in img.nii
if ~exist('imgName')
 imgName = spm_select(1,'image','Select 4D image to crop');
end
[pth,nam,ext] = spm_fileparts( deblank (imgName)); %remove volume index, 'img.nii,1' ->'img.nii'
imgName = fullfile(pth,[nam, ext]); %remove volume is 'img.nii,1' -> 'img.nii'
hdr = spm_vol([imgName]); 
nVol = length(hdr);
nSlices = hdr(1).dim(3);
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

function isSPM12orNewerSub
%check that SPM is installed and is at least release 6225
if exist('spm','file') ~= 2, error('Please install SPM12 or later'); end;
[v,r] = spm('Ver','',1); r = str2double(r); %#ok<ASGLU>
if r < 6225, error('Please update your copy of SPM'); end;
%end isSPM12orNewer()

function coivox = setOriginSub(vols, modality)
%Align images so that origin and alignment roughly match MNI space
%  vols : cell string of image name(s) - first image used for estimate, others yoked
%  modality : modality of first image 1=T1, 2=T2, 3=EPI
%Example
%  setOrigin('T1.nii',1); %align T1 scan
%  setOrigin({'T1s005.nii', 'fmriblocks009.nii'},1); %use T1 to align T1 and fMRI data
%  setOrigin %use graphical interface
%Chris Rorden 12/2014 (now supports SPM12)
if ~exist('vols','var') %no files specified
 vols = spm_select(inf,'image','Reset origin for selected image(s) (estimated from 1st)');
end
if ischar(vols)
    vols = cellstr(vols);
end
if ~exist('modality','var') %no files specified
 modality = 1;
 fprintf('%s Modality not specified, assuming T1\n', mfilename);
end
coivox = ones(4,1);
%extract filename 
[pth,nam,ext, ~] = spm_fileparts(deblank(vols{1}));
fname = fullfile(pth,[nam ext]); %strip volume label
%report if filename does not exist...
if (exist(fname, 'file') ~= 2) 
 	fprintf('%s error: unable to find image %s.\n',mfilename,fname);
	return;  
end;
hdr = spm_vol([fname,',1']); %load header 
img = spm_read_vols(hdr); %load image data
img = img - min(img(:));
img(isnan(img)) = 0;
%find center of mass in each dimension (total mass divided by weighted location of mass
% img = [1 2 1; 3 4 3];
sumTotal = sum(img(:));
coivox(1) = sum(sum(sum(img,3),2)'.*(1:size(img,1)))/sumTotal; %dimension 1
coivox(2) = sum(sum(sum(img,3),1).*(1:size(img,2)))/sumTotal; %dimension 2
coivox(3) = sum(squeeze(sum(sum(img,2),1))'.*(1:size(img,3)))/sumTotal; %dimension 3
XYZ_mm = hdr.mat * coivox; %convert from voxels to millimeters
fprintf('%s center of brightness differs from current origin by %.0fx%.0fx%.0fmm in X Y Z dimensions\n',fname,XYZ_mm(1),XYZ_mm(2),XYZ_mm(3)); 
for v = 1:   numel(vols) 
    fname = deblank(vols{v});
    if ~isempty(fname)
        [pth,nam,ext, ~] = spm_fileparts(fname);
        fname = fullfile(pth,[nam ext]); 
        hdr = spm_vol([fname ',1']); %load header of first volume 
        fname = fullfile(pth,[nam '.mat']);
        if exist(fname,'file')
            destname = fullfile(pth,[nam '_old.mat']);
            copyfile(fname,destname);
            fprintf('%s is renaming %s to %s\n',mfilename,fname,destname);
        end
        hdr.mat(1,4) =  hdr.mat(1,4) - XYZ_mm(1);
        hdr.mat(2,4) =  hdr.mat(2,4) - XYZ_mm(2);
        hdr.mat(3,4) =  hdr.mat(3,4) - XYZ_mm(3);
        spm_create_vol(hdr);
        if exist(fname,'file')
            delete(fname);
        end
    end
end%for each volume
coregEstTemplateSub(vols, modality);
for v = 1:   numel(vols) 
    [pth, nam, ~, ~] = spm_fileparts(deblank(vols{v}));
    fname = fullfile(pth,[nam '.mat']);
    if exist(fname,'file')
        delete(fname);
    end
end %for each volume
%end setOriginSub()

function coregEstTemplateSub(vols, modality)
%vols: images to coregister - first used for estimate
if modality == 2
   template = fullfile(spm('Dir'),'canonical','avg152T2.nii');
elseif modality == 3
    template  = fullfile(spm('Dir'),'toolbox','OldNorm','EPI.nii');
else
    template = fullfile(spm('Dir'),'canonical','avg152T1.nii');
end
if ~exist(template,'file')
    error('Unable to find template named %s\n', template);
end
if ischar(vols)
    vols = cellstr(vols);
end
matlabbatch{1}.spm.spatial.coreg.estimate.ref = {template};
matlabbatch{1}.spm.spatial.coreg.estimate.source = {[deblank(vols{1}),',1']};%{'/Users/rorden/Desktop/3D.nii,1'};
if  numel(vols) > 1
    matlabbatch{1}.spm.spatial.coreg.estimate.other = vols(2:end);% {''};
else
    matlabbatch{1}.spm.spatial.coreg.estimate.other = {''};
end
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];
spm_jobman('run',matlabbatch);
%end coregEstTemplateSub()

function fnm = zeroHemisphereSub (fnm, zeroLeft, meanname)
%set all voxels in one hemisphere to zero
% fnm: name of 3D volume
% zeroLeft: if true left hemisphere masked, else right
%Example:
% zeroHemisphere('c2eT1_P022.nii', true)
% zeroHemisphere; %use GUI
if ~exist('fnm','var')
	fnm = spm_select(1,'image','Select image[s] for NaN removal'); 
end
hdr = spm_vol(fnm);
img = spm_read_vols(hdr);
[pth nm ext] = spm_fileparts(fnm);
if ~exist('zeroLeft','var') || zeroLeft
    fnm = fullfile(pth, ['r' nm ext]); %left masked, so 'r'ight survivies 
    img(1:floor(size(img,1)/2),:,:) = 0;
else
    fnm = fullfile(pth, ['l' nm ext]); %right masked, so 'l'eft survives
    img(ceil(size(img,1)/2):end,:,:) = 0;
end
hdr.fname = fnm;
img(isnan(img)) = 0;%max(img(:)); % use ~isfinite instead of isnan to replace +/-inf with zero
spm_write_vol(hdr,img);
if ~exist('meanname','var') || isempty(meanname), return; end;
fnm =coregWriteSub(fnm, meanname);
%end zeroHemisphere()            

function mean = maskMeanSub(fnm, roi)
%report mean intensity of image 'fnm' modulated by roi
%Example
% maskMean('meanCBF_0_srcASL_P022.nii','rzc1eT1_P022.nii')
% maskMean('meanCBF_0_srcASL_P022.nii','rzc2eT1_P022.nii')
hdr = spm_vol(fnm);
img = spm_read_vols(hdr);
rhdr = spm_vol(roi);
rimg = spm_read_vols(rhdr);
if (max(rimg(:)) > (1.01)) || (min(rimg(:)) < -0.01)
    fprintf('Regions of interest should have voxels in the range 0..1 %s\n', roi);
    mean = 0;
    return;
end
img = img .* rimg;
mean =  sum(img(:)) / sum(rimg(:));
%end maskMeanSub()