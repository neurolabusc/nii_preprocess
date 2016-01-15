function nii_lime
mm = 3; %reslicing 


baseDir = '/Users/rorden/desktop/LIME/Processing/';

preprocess = true;
preprocessDti = true;

% subj_names = strvcat('LM1001','LM1002','LM1003','LM1004','LM1005','LM1006','LM1007',...
% 'LM1008','LM1009','LM1010','LM1011','LM1014','LM1015','LM1016',...
% 'LM1017','LM1018','LM1019','LM1020','LM1021','LM1022','LM1023',...
% 'LM1024','LM1025','LM1026','LM1027','LM1028','LM1029','LM1030',...
% 'LM1031','LM1032','LM1033','LM1034','LM1035','LM1036','LM1037','LM1039','LM1040');

subj_names = strvcat('LM1001');
%[pth,~,~] = spm_fileparts( deblank (which(mfilename)));
%nameTabFile = fullfile(pth,['LIME_12_16_2013.xlsx']); %excel file
t1Prefix = 'T1_';
lesionPrefix = 'LS_';
t2Prefix = 'T2_'; 
restPrefix = 'REST_';
paslPrefix = 'PASL_'; 
dtiAPrefix = 'DTIA_';%phase in A->P
dtiPPrefix = 'DTIP_';%phase in P->A
restTRsec = 1.85;
restSliceOrder = 2; %2=descending
spm_jobman('initcfg');
spm fmri %start SPM so .ps files will be created for preprocessing


for s = 1:   size(subj_names,1) 
    subj = deblank(subj_names(s,:));
    pth = [baseDir,subj filesep];
    nameMat = fullfile(pth,[ subj, '.mat']);
    nameT1 = checkFilename (pth, t1Prefix, subj);
    nameT2 = checkFilename (pth, t2Prefix, subj);
    nameLesion = checkFilename (pth, lesionPrefix, subj);
    nameRest = checkFilename (pth, restPrefix, subj);
    namePasl = checkFilename (pth, paslPrefix, subj);
    nameDtiA = checkFilename (pth, dtiAPrefix, subj);
    nameDtiP = checkFilename (pth, dtiPPrefix, subj);
    if ~lesionMatchT2Sub (nameT2,nameLesion);
        return
    end
    if preprocess  
        
        %normalize T1 -   
        if ~skipIfFilename(pth, ['m' t1Prefix], subj)
            %next set origin
            originSub(nameT1,nameT2,nameLesion,nameRest,namePasl,nameDtiA,nameDtiP);
            clinical_mrnormseg (nameT1,nameLesion,nameT2, true, [mm mm mm; 1 1 1], [-78 -112 -50; 78 76 85], false, 0.005, 2);  
        end
        if ~skipIfFilename(pth, ['zw1m' t1Prefix], subj) %compute i3m maps
            nii_i3m(appendPrefixSub('w1m',nameT1),maskSub,0.5,10,0.25,1); %i3m T1 image
            nii_i3m(appendPrefixSub('w1',nameT2),maskSub,0.5,10,0.25,2); %i3m T2 image
        end
        if preprocessDti
        	if false %exist([pth 'dti']) ~= 0 %#ok<UNRCH>
        		fprintf('Skipping dti preprocessing for %s as dti folder exists\n',subj);
        	else
            	if (~isempty(nameDtiA)) && (~isempty(nameDtiP))
               		%process DTI 
                    ditpth = [pth 'dti']; %save files in new directory
                    nameLesionSR = checkFilename (pth,['sr' lesionPrefix], subj);
                    dtiSub(nameDtiA,nameDtiP,nameT1,nameLesionSR, ditpth);
               		
               		%
               		%nii_dti_prep (nameDtiA,nameDtiP,nameT1,nameLesionSR, ditpth);
            	end
        	end
        end %if preprocessDti
        %preprocess ASL data
        if (~isempty(namePasl)) &&  ~skipIfFilename(pth, ['rc' paslPrefix], subj)
            pasl_batch (namePasl, nameT1, 1.1,mm); %1.1s = "inversion time 2 = 1800ms"-"inversion time 1 = 700ms" 
        end
        if (~isempty(nameRest)) &&  ~skipIfFilename(pth, ['fdwa' restPrefix], subj)
            nii_rest_batch(nameRest,nameT1,restTRsec,restSliceOrder);
        end
    end; %if preprocess
    error('xx');
    extract2MatSub( nameT1,nameT2, nameLesion,nameRest,namePasl, nameMat);  
    moveProcessedFilesSub(baseDir,subj,nameMat,t1Prefix,lesionPrefix,t2Prefix,restPrefix,paslPrefix);
end %for each subject
%end nii_lime()

function dtiSub(dtia,dtib,t1,lesion,newPath)
if isempty(newPath) 
    newPath = fullfile(pwd, 'temp');
    mkdir(newPath);
end
if (exist('newPath','var')) && (~isempty(newPath))
    dtia = cpImgSub(newPath,dtia);
    dtib = cpImgSub(newPath,dtib);
    t1 = cpImgSub(newPath,t1);
    lesion = cpImgSub(newPath,lesion);
end
fsldir= '/usr/local/fsl/';
if ~exist(fsldir,'dir'), error('%s: fsldir (%s) not found',mfilename,fsldir); end
flirt = [fsldir 'bin/flirt'];
if ~exist(flirt,'file')
	error('%s: fsl not installed (%s)',mfilename,flirt);
end
setenv('FSLDIR', fsldir);
curpath = getenv('PATH');
setenv('PATH',sprintf('%s:%s',fullfile(fsldir,'bin'),curpath));
%./dti.sh "DTIA_LM1001" "" "LS_LM1001" "T1_LM1001"
command= [fileparts(which(mfilename)) filesep 'dti.sh'];
command=sprintf('sh -c ". %setc/fslconf/fsl.sh; %s "%s" "%s" "%s" "%s""\n',fsldir,command, dtia,dtib,lesion,t1);
command
system(command);
%end dtiSub

function newName = cpImgSub(newPath,oldName);
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


function originSub (nameT1,nameT2,nameLesion,nameRest,namePasl,nameDtiA,nameDtiP)
vols = nameT1;
if (~isempty(nameT2)), vols = strvcat(vols, nameT2); end; %#ok<*REMFF1>
if (~isempty(nameLesion)), vols = strvcat(vols, nameLesion); end;
if (~isempty(nameRest)), vols = strvcat(vols, nameRest); end;
if (~isempty(namePasl)), vols = strvcat(vols, namePasl); end;
if (~isempty(nameDtiA)), vols = strvcat(vols, nameDtiA); end;
if (~isempty(nameDtiP)), vols = strvcat(vols, nameDtiP); end;
setOriginSub(vols, 1);
%end originSub()            

function vols = vol1OnlySub(vols)
%only select first volume of multivolume images '/dir/img.nii' -> '/dir/img.nii,1', '/dir/img.nii,33' -> '/dir/img.nii,1'
oldvols = vols;
vols = [];
for v = 1:   size(oldvols,1) 
    [pth,nam,ext, ~] = spm_fileparts(deblank(oldvols(v,:)));
    vols = strvcat(vols, fullfile(pth, [ nam ext ',1']) ); %#ok<REMFF1>
end
%end vol1OnlySub()

function coivox = setOriginSub(vols, modality)
% see nii_detect_Origin
%Example
%  nii_centerinten('4D.nii');

coivox = ones(4,1);
if ~exist('vols','var') %no files specified
 vols = spm_select(inf,'image','Reset origin for selected image(s) (estimated from 1st)');
end
vols = vol1OnlySub(vols); %only process first volume of 4D datasets...
if ~exist('modality','var') %no files specified
 modality = 1;
 fprintf('%s Modality not specified, assuming T1\n', mfilename);
end
%extract filename 
[pth,nam,ext, ~] = spm_fileparts(deblank(vols(1,:)));
fname = fullfile(pth,[nam ext]); %strip volume label
%report if filename does not exist...
if (exist(fname, 'file') ~= 2) 
 	fprintf('%s set origin error: unable to find image %s.\n',mfilename,fname);
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
for v = 1:   size(vols,1) 
    fname = deblank(vols(v,:));
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
coregSub(vols, modality);
for v = 1:   size(vols,1) 
    [pth, nam, ~, ~] = spm_fileparts(deblank(vols(v,:)));
    fname = fullfile(pth,[nam '.mat']);
    if exist(fname,'file')
        delete(fname);
    end
end%for each volume
%end setOriginSub()

function coregSub(vols, modality)
if modality == 2
    template = fullfile(spm('Dir'),'templates','T2.nii');
elseif modality == 3
    template  = fullfile(spm('Dir'),'toolbox','Clinical','scct.nii');
else
    template = fullfile(spm('Dir'),'templates','T1.nii');
end
if ~exist(template,'file')
    error('Unable to find template named %s\n', template);
end
matlabbatch{1}.spm.spatial.coreg.estimate.ref = {template};
matlabbatch{1}.spm.spatial.coreg.estimate.source = {[deblank(vols(1,:)),',1']};%{'/Users/rorden/Desktop/3D.nii,1'};
matlabbatch{1}.spm.spatial.coreg.estimate.other = cellstr(vols(2:end,:));% {''};
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];
spm_jobman('run',matlabbatch);
%end coregSub()

function extract2MatSub( nameT1,nameT2,nameLesion,nameRest,namePasl, nameMat)
%behavioralData2MatSub(nameTabFile,subj,nameMat); %extract behavioral data

i3mVoxels2MatSub(appendPrefixSub ('w1m',nameT1),1,nameMat); %voxelswise T1
i3mVoxels2MatSub(appendPrefixSub ('w1',nameT2),2,nameMat);
lesionVoxels2MatSub (appendPrefixSub ('bwsr',nameLesion),nameMat);

kROIs = nii_roi_list ();
for roiIndex = 1: size(kROIs,1)
    template = deblank(kROIs(roiIndex,:));
    nii_roi2stats(template,appendPrefixSub ('zw1m',nameT1), '','i3mT1_',nameMat);
    nii_roi2stats(template,appendPrefixSub ('zw1',nameT2), '','i3mT2_',nameMat);
    nii_roi2stats(template,appendPrefixSub ('bw1sr',nameLesion),'','lesion_',nameMat);
    nii_roi2stats(template, appendPrefixSub ('alf_dwa',nameRest), '','rest_',nameMat);
    nii_roi2stats(template, appendPrefixSub ('fdwa',nameRest), '', 'rest_',nameMat);
    nii_roi2stats(template,appendPrefixSub ('mskwmeanCBF_0_src',namePasl),'','cbf_',nameMat);
end

%end extract2Mat()

function moveProcessedFilesSub(baseDir,subj,nameMat,t1Prefix,lesionPrefix,t2Prefix,restPrefix,paslPrefix);
%archive final images into a separate folder
pth = [baseDir,subj filesep];
fpth = makeDirSub([baseDir 'Standard']);
fpth = makeDirSub([fpth subj]);
moveFileSub(nameMat,fullfile(fpth,[ subj, '.mat']) );
moveFileSub(fullfile(pth,['renderw1m' t1Prefix subj '.nii']), fullfile(fpth,['b' t1Prefix subj '.nii']) );
moveFileSub(fullfile(pth,['w1m' t1Prefix subj '.nii']), fullfile(fpth,[ t1Prefix subj '.nii']) );
moveFileSub(fullfile(pth,['zw1m' t1Prefix subj '.nii']), fullfile(fpth,['z' t1Prefix subj '.nii']) );
moveFileSub(fullfile(pth,['w1sr' lesionPrefix subj '.nii']), fullfile(fpth,[ lesionPrefix subj '.nii']) );
moveFileSub(fullfile(pth,['w1sr' lesionPrefix subj '.nii']), fullfile(fpth,[ lesionPrefix subj '.nii']) );
moveFileSub(fullfile(pth,['w1' t2Prefix subj '.nii']), fullfile(fpth,[ t2Prefix subj '.nii']) );
moveFileSub(fullfile(pth,['zw1' t2Prefix subj '.nii']), fullfile(fpth,['z' t2Prefix subj '.nii']) );
moveFileSub(fullfile(pth,['fdwa' restPrefix subj '.nii']), fullfile(fpth,[ restPrefix subj '.nii']) );
moveFileSub(fullfile(pth,['alf_dwa' restPrefix subj '.nii']), fullfile(fpth,['alf' restPrefix subj '.nii']) );
moveFileSub(fullfile(pth,['mskwmeanCBF_0_src' paslPrefix subj '.nii']), fullfile(fpth,[paslPrefix subj '.nii']) );        
%end moveFilesSub

function outPath = makeDirSub(Path)
if exist(Path) ~= 7
    mkdir( Path );
end
outPath = [Path filesep];
%makeDir()

function moveFileSub(src, dest)
if exist(src) ~= 2
    %fprintf('Unable to find file %s\n',src);
    return;
end
copyfile(src,dest);
%end moveFileSub()

function i3mVoxels2MatSub (nameNormalizedAnat,modality,diskName)
%filename ('wT1_LM1001.nii,1,'LM1001.mat');  ('wT2_LM1001.nii,2,'LM1001.mat');
if exist(nameNormalizedAnat) ~= 2 
    fprintf('%s warning: unable to find image named %s\n',mfilename,nameNormalizedLesion);
    return;
end
[hdr,img] = nii_i3m(nameNormalizedAnat,maskSub,0.5,10,0.25,modality,2);
stat = [];
fieldname = ['i3mT' int2str(modality)];
stat.(fieldname).hdr = hdr;
stat.(fieldname).dat = img;
if exist(diskName)
    old = load(diskName);
    stat = nii_mergestruct(stat,old);
end
save(diskName,'-struct', 'stat');
%end i3mData2Mat()

function lesionVoxels2MatSub (nameNormalizedLesion,diskName)
if exist(nameNormalizedLesion) ~= 2 
    fprintf('%s warning: unable to find image named %s\n',mfilename,nameNormalizedLesion);
    return;
end
hdr = spm_vol(nameNormalizedLesion); %load dataset only once!
stat = [];
stat.lesion.hdr = hdr;
stat.lesion.dat = spm_read_vols(hdr);
if exist(diskName)
    old = load(diskName);
    stat = nii_mergestruct(stat,old);
end
save(diskName,'-struct', 'stat');
%end lesionData2Mat()

function behavioralData2MatSub (tabName,subjName,diskName);
[~,~,ext] = spm_fileparts( deblank (which(tabName)));
if strcmpi(ext,'.tab')
	error('Please check that nii_tab2mat still works... many changes');
    stat.behav = nii_tab2mat(tabName,subjName);
else
    if length(length(subjName) > 2) && (strcmpi(subjName(1),'L')) && (strcmpi(subjName(1),'L'))
        subjName = subjName(3:length(subjName));
    end
	stat.behav = nii_xls2mat(tabName , 'Data (2)', subjName);
end
if (length(stat) < 1) || (length(diskName) < 1), return; end
if exist(diskName)
    old = load(diskName);
    stat = nii_mergestruct(stat,old);
end
save(diskName, '-struct', 'stat');
%end behavioralData2Mat()

function dimsMatch = lesionMatchT2Sub (T2,lesion)
dimsMatch = true;
if (length(T2) < 1) || (length(lesion) < 1), return; end
lhdr = spm_vol(lesion); %lesion header
t2hdr = spm_vol(T2); %pathological scan header
if ~isequal(lhdr.dim,t2hdr.dim);
    dimsMatch = false;
    fprintf('%s ERROR: Dimension mismatch %s %s: %dx%dx%d %dx%dx%d\n',mfilename, T2,lesion, t2hdr.dim(1),t2hdr.dim(2),t2hdr.dim(3), lhdr.dim(1),lhdr.dim(2),lhdr.dim(3));  
end
%end dimsMatch()

function newName =appendPrefixSub (prefix,origName)
[pth,nam,ext] = spm_fileparts( deblank (origName));
newName = fullfile(pth,[prefix nam ext]);
%end appendPrefix

function isExisting = skipIfFilename (pth, prefix, nam);
imgName = fullfile(pth,[prefix nam '.nii']);
isExisting = (exist(imgName) == 2);
if isExisting
    fprintf('%s exists: will skip preprocessing step.\n',imgName);
end
%end skipIfFilename

function msk = maskSub()
msk= [fileparts(which(mfilename)) filesep 'scmask.nii'];
%end maskSub()

function img = checkFilename (pth, prefix, nam)
img = fullfile(pth,[prefix nam '.nii']);
if exist(img) ~= 2 
    fprintf('%s warning: unable to find image named %s\n',mfilename,img);
    img = '';
end  
%end checkFilename()