function nii_findimg_Mar2015 (baseDir)
%find relevant images for patients
% baseDir: root folder that holds images
%Example Structure (%s is basefolder, with two subjects)
% ~\p001\T2_P001.nii
% ~\p001\LS_P001.nii
% ~\p001\T1_P001.nii
% ~\p001\fMRI_P001.nii
% ~\p002\T2_P002.nii
% ~\p002\LS_P002.nii
% ~\p002\T1_P002.nii
%Example
% nii_findimg_q %use gui
% nii_findimg_q(pwd)
% nii_findimg_q('/Users/rorden/Desktop/pre')

if ~exist('baseDir','var') || isempty(baseDir)
    baseDir = pwd; %uigetdir('','Pick folder that contains all subjects');
end
subjDirs = subFolderSub(baseDir);
subjDirs = sort(subjDirs);
%subjDirs = {'LM1003'};
% 'LM1052', 'LM1054', 'LM1056RHD', 'P001', 'P011', 'P012', 'P021', 'P023', 'P024', 'P026', 'P027', 'P098', 'P106', 'P134', 'P171'

%subjDirs = {'LM1052', 'LM1054', 'LM1056RHD', 'P001', 'P011', 'P012', 'P021', 'P023', 'P024', 'P026', 'P027', 'P098', 'P106', 'P134', 'P171'
%subjDirs = {'LM1052', 'LM1054', 'LM1056RHD'}
nSubjDir = numel(subjDirs);
%nSubjDir2 = floor(numel(subjDirs) /2);
fslParallelSub;
fprintf('Processing %d folders (subjects) in %s\n', nSubjDir, baseDir);
%fprintf('Processing %d of %d folders (subjects) in %s\n',nSubjDir2,nSubjDir, baseDir);
for s = nSubjDir: -1 : 1 %1:nSubjDir2 %(nSubjDir2+1):nSubjDir    
    subjDir = [baseDir,filesep, deblank(subjDirs{s}) ]; %no filesep
    matName = [subjDirs{s} '.mat'];
    %reportMatDetails(matName); continue;
    
    imgs = subjStructSub(deblank(subjDirs{s}));
    imgs.Lesion = imgfindSub(imgs.Lesion,strvcat('LS_'), subjDir); %#ok<REMFF1>
    imgs.T1 = imgfindSub(imgs.T1,'T1_', subjDir); 
    imgs.T2 = imgfindSub(imgs.T2,'T2_', subjDir); 
    imgs.ASL = imgfindSub(imgs.ASL,'ASL_', subjDir);
    imgs.fMRI = imgfindSub(imgs.fMRI,'fMRI_', subjDir);
    imgs.Rest = imgfindSub(imgs.fMRI,'rest_', subjDir);
    imgs.fMRI = imgfindSub(imgs.fMRI,'fMRI_', subjDir);
    imgs.DTI = imgfindSub(imgs.DTI,'APDTI_', subjDir);
    if isempty(imgs.DTI) %cannot find 'APDTI_' look for standard DTI
        imgs.DTI = imgfindSub(imgs.DTI,'DTI_', subjDir);
    else
        imgs.DTIrev = imgfindSub(imgs.DTIrev,'PADTI_', subjDir);
    end
    imgs = unGzAllSub(imgs); %all except DTI - fsl is OK with nii.gz
    if ~isempty(imgs.T1) %&& ~isempty(imgs.Lesion) % ~isempty(imgs.DTI)
        
        %doT1Sub(imgs, matName); %normalize T1
        %doI3MSub(imgs, matName);
        %doDtiSub(imgs, matName);
        %
        %fprintf('completed\t%s\n', subjDirs{s});
        
        %if s > 4, error('123'); end;
        %error('123 %s', matName)
        doAslSub(imgs, matName);
        %doRestSub(imgs, matName, 1.85, 2); %TR= 1.850 sec, descending; 
        %doRestSub(imgs, matName, 2.05, 5); %Souvik study - 2.0sec, interleaved even number of slices
        %printSub(imgs,deblank(subjDirs{s}),true );
        
        %dofMRISub(imgs, matName);
        %doFaMdSub(imgs, matName);
        %reportMatDetails(matName);
        
    end
    %checkDims(imgs);
end
%end nii_findimg()

function reportMatDetails(matName)
file_list = nii_modality_list;
fldName = 'Participant';
[p,n] = fileparts(matName);
fldExist = n;
for i = 1 : size(file_list,1)
    fld = deblank(file_list(i,:));
    if strcmpi(fld,'rest') || strcmpi(fld,'dtifc') || strcmpi(fld,'dti') %connectivity maps always ROI specific
        fld = [fld, '_jhu'];
    end
    fldName = [fldName, sprintf('\t%s',fld)];
    x = isFieldSub(matName,fld);
    fldExist = [fldExist, sprintf('\t%d',x)];
    %fprintf('Field %s exists: %d\n', fld, x);
end
logname = 'modalities.txt';
logexists = exist(logname, 'file');
fid = fopen(logname,'a');
if ~logexists
    fprintf(fid, sprintf('%s\n',fldName));
end
fprintf(fid, sprintf('%s\n',fldExist));
fclose(fid);
%end reportMatDetails()

function dofMRISub(imgs, matName)
if isempty(imgs.T1) || isempty(imgs.fMRI), return; end; %required
[~,n] = fileparts(matName);
[p] = fileparts(imgs.fMRI);
cstat = fullfile(p,['fMRI_', n], 'con_0002.nii');
bstat = fullfile(p,['fMRI_', n], 'beta_0001.nii');
if ~exist(cstat, 'file') || ~exist(bstat,'file')
    nii_fmri60(imgs.T1, imgs.fMRI); %666
end
if ~exist(cstat, 'file') || ~exist(bstat,'file')
    error('fMRI analysis failed');
end
if isFieldSub(matName, 'fmri'), return; end; %stats already exist
nii_nii2mat(cstat, 'fmri' , matName); %12
nii_nii2mat(bstat, 'fmrib', matName); %13
%end doI3MSub()


function printSub(imgs, nam, doT1)
if isempty(imgs.T1) || isempty(imgs.Lesion), return; end; %required
T1 = prefixSub('wb',imgs.T1); %warped brain extracted image
LS = prefixSub('wsr', imgs.Lesion);
if ~exist(T1,'file') || ~exist(LS,'file') , return; end; %required
if ~exist('doT1', 'var') || doT1
    spm_clf;
    spm_figure('Clear', 'Graphics');
    spm_orthviews('Reset');
    fig = spm_figure('FindWin','Graphics');
    ax  = axes('Position',[0.1 0.5 0.35 0.3],'Visible','off','Parent',fig);
    spm_orthviews('Image',LS,[0.01 0.01 .48 .49]);
    spm_orthviews('Image',T1,[0.51 0.01 .48 .49]);
    fMRI = prefixSub('sw', imgs.fMRI);
    if exist(fMRI, 'file')
        spm_orthviews('Image',fMRI,[0.01 0.38 .48 .49]);
    end
    i3m = prefixSub('zw',imgs.T1); %warped brain extracted image
    if exist(i3m, 'file')
        spm_orthviews('Image',i3m,[0.51 0.38 .48 .49]);
    end
    [~, nam] = fileparts(T1);
    text(0,0.97, sprintf('normalized %s',nam),'Parent',ax, 'FontSize',24, 'fontn','Arial');
    XYZmm = getCenterOfIntensitySub(LS);
    spm_orthviews('setcoords',XYZmm)
    spm_print;
    return;
end
%next DTI
if  isempty(imgs.DTI) , return; end; %required
FA = prepostfixSub('n', '_FA', imgs.DTI);
ROI = prepostfixSub('', '_roi', imgs.DTI);
if ~exist(FA,'file') || ~exist(ROI,'file') , return; end; %required
FA = unGzSub (FA);
ROI = unGzSub (ROI);
spm_clf;
spm_figure('Clear', 'Graphics');
spm_orthviews('Reset');
fig = spm_figure('FindWin','Graphics');
ax  = axes('Position',[0.1 0.5 0.35 0.3],'Visible','off','Parent',fig);
spm_orthviews('Image',FA,[0.01 0.33 .48 .49]);
spm_orthviews('Image',ROI,[.51 0.33 .48 .49]);
[~, nam] = fileparts(FA);
text(0,0.95, sprintf('DTI %s',nam),'Parent',ax, 'FontSize',24, 'fontn','Arial');
spm_print;
%end printSub()
 
function XYZmm = getCenterOfIntensitySub(vols)
XYZmm = ones(3,1);
if ischar(vols), vols = cellstr(vols); end;
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
coivox = ones(4,1);
coivox(1) = sum(sum(sum(img,3),2)'.*(1:size(img,1)))/sumTotal; %dimension 1
coivox(2) = sum(sum(sum(img,3),1).*(1:size(img,2)))/sumTotal; %dimension 2
coivox(3) = sum(squeeze(sum(sum(img,2),1))'.*(1:size(img,3)))/sumTotal; %dimension 3
XYZmm = hdr.mat * coivox; %convert from voxels to millimeters
XYZmm = XYZmm(1:3);
%end setCenterOfIntensitySub()

function doFaMdSub(imgs, matName)
if isempty(imgs.T1) || isempty(imgs.DTI), return; end; %required
T1 = prefixSub('wb',imgs.T1); %warped brain extracted image
FA = prepostfixSub('', '_FA', imgs.DTI);
MD = prepostfixSub('', '_MD', imgs.DTI);
if ~exist(T1,'file') || ~exist(FA,'file') || ~exist(MD,'file') , return; end; %required
FA = unGzSub (FA);
MD = unGzSub (MD);
wFA = prepostfixSub('w', '', FA);
wMD = prepostfixSub('w', '', MD);
if ~exist(wFA,'file') || ~exist(wMD,'file')
    nFA = rescaleSub(FA);
    %nii_setOrigin12({nFA, FA, MD}, 1,false); %rFA looks like T1
    oldNormSub( {nFA, FA,MD}, T1, 8, 10 );
end

if isFieldSub(matName, 'fa'), return; end; %stats already exist
nii_nii2mat(wFA, 'fa', matName); %6
nii_nii2mat(wMD, 'md', matName); %8
%end doFaMdSub()

function nam = prepostfixSub (pre, post, nam)
[p, n, x] = filepartsSub(nam);
nam = fullfile(p, [pre, n, post, x]);
if ~exist(nam, 'file') && strcmpi(x,'.nii') 
   nam = fullfile(p, [pre, n, post, '.nii.gz']);
   if ~exist(nam, 'file') 
        nam = fullfile(p, [pre, n, post, x]);
   end
end
if ~exist(nam, 'file') && strcmpi(x,'.nii.gz') 
   nam = fullfile(p, [pre, n, post, '.nii']);
   if ~exist(nam, 'file') 
        nam = fullfile(p, [pre, n, post, x]);
   end
end
%end prefixSub()

function [p,n,x] = filepartsSub (fnm)
if isempty(fnm), return; end;
fnm = char(deblank(fnm));
[p,n,x] = fileparts(fnm);
if strcmpi(x,'.gz') %.nii.gz
    [p2,n2,x2] = spm_fileparts(n);
    n = n2;
    x = [x2, x];
end;
%end filepartsSub()

function fname = rescaleSub(fname)
hdr = spm_vol(fname);
img = spm_read_vols(hdr);
if min(img(:)) < 0, error('Image intensity must be positive'); end;
if min(img(:)) == max(img(:)), error('Image has no variability'); end;
img = img - min(img(:)); %scale from zero
img = img/max(img(:)); %scale from 0..1
[pth nm ext] = spm_fileparts(fname);
img = power(img, 0.5);
fname = fullfile(pth, ['n' nm ext]);
hdr.fname = fname;  
spm_write_vol(hdr,img);
%end rescaleSub()

function checkDims(imgs)
if isempty(imgs.T1) || isempty(imgs.Lesion), return; end; %we need these images
hdr = spm_vol(imgs.Lesion);
[pth nm] = spm_fileparts(imgs.Lesion);
fprintf('%s\t%d\t%d\t%d\n', nm, hdr.dim(1), hdr.dim(2), hdr.dim(3) );
%end checkDims()


function doDtiSub(imgs, matName)
if isempty(imgs.T1) || isempty(imgs.DTI), return; end; %required
betT1 = prefixSub('b',imgs.T1); %brain extracted image
eT1 = prefixSub('e',imgs.T1); %enantimorphic image
if ~exist(betT1,'file') || ~exist(eT1,'file'), return; end; %required
if isDtiDone(imgs), fprintf('Skipping DTI processing (probtrackx done)'); return; end;
n = bvalCountSub(imgs.DTI);
if (n < 1)
    fprintf('UNABLE TO FIND BVECS/BVALS FOR %s\n', imgs.DTI);
    return
end
if (n < 12)
    fprintf('INSUFFICIENT BVECS/BVALS FOR %s\n', imgs.DTI);
    return
end
dti_u=prepostfixSub('', 'u', imgs.DTI);
if exist(dti_u, 'file')
    fprintf('Skipping DTI preprocessing: found %s\n', dti_u);
else
    clipSub (imgs); %topup requires images with even dimensions
    %1 - eddy current correct
    command= [fileparts(which(mfilename)) filesep 'dti_1_eddy.sh'];
    if isempty(imgs.DTIrev)
        command=sprintf('%s "%s"',command, imgs.DTI);
    else
        nr = bvalCountSub(imgs.DTIrev);
        if (nr ~= n)
            fprintf('BVECS/BVALS DO NOT MATCH %s %s\n', imgs.DTI, imgs.DTIrev);
            return    
        end
        command=sprintf('%s "%s" "%s"',command, imgs.DTI, imgs.DTIrev);
    end
    cleanupDtiDir(imgs.DTI);
    doFslCmd (command);
end
doDtiBedpostSub(imgs.DTI);
doDtiWarpJhuSub(imgs); %warp atlas to DTI
doDtiTractSub(imgs.DTI); %tractography
%end doDtiSub()

function done = isDtiDone(imgs)
done = false;
if isempty(imgs.DTI), return; end;
p = fileparts( imgs.DTI );
pDir = fullfile(p,'probtrackx');
if exist(pDir, 'file')
    done = true; 
end;
%end isDtiDone()

function clipSub (imgs)
if isempty(imgs.DTIrev)
    nii_clipeven(imgs.DTI, true);
else
    nii_clipeven({imgs.DTI, imgs.DTIrev}, true);
end
%end clipSub()

% function doDtiBedpostSub(dti) %warp atlas to DTI
% t_start=tic;
% pth = fileparts(dti);
% bed_dir=fullfile(pth, 'bedpost');
% if ~exist(bed_dir, 'file'), mkdir(bed_dir); end;
% bed_dirX=fullfile(pth, 'bedpost.bedpostX');
% if ~exist(bed_dirX, 'file'), mkdir(bed_dirX); end;
% bed_done=fullfile(bed_dirX, 'xfms', 'eye.mat');
% if exist(bed_done,'file'), fprintf('Skipping bedpost (already done)\n'), return, end;
% dti_u=prepostfixSub('', 'u', dti);
% dti_x=fullfile(bed_dir, 'data.nii.gz');
% if ~exist(dti_u,'file'), error('Bedpost unable to find %s', dti_u); end;
% copyfile(dti_u, dti_x);
% [bvec, bval] = getBVec(dti);
% dti_x=fullfile(bed_dir, 'bvecs');
% copyfile(bvec, dti_x);
% dti_x=fullfile(bed_dir, 'bvals');
% copyfile(bval, dti_x);
% dti_faThr=prepostfixSub('', '_FA_thr', dti);
% dti_x=fullfile(bed_dir, 'nodif_brain_mask.nii.gz');
% copyfile(dti_faThr, dti_x);
% command = fullfile(fslDirSub, 'bin', 'bedpostx1');
% if ~exist(command,'file'), error('You need the special file %s', command); end;
% command=sprintf('bedpostx1 "%s" ', bed_dir);
% doFslCmd (command);
% bed_cmd=fullfile(bed_dirX, 'commands.txt');
% if ~exist(bed_cmd, 'file'), error('Unable to find %s', bed_cmd); end;
% maxThreads = feature('numCores');
% if maxThreads > 15, maxThreads = maxThreads - 1; end;
% nPerm = 5000;
% %cmds = textread(bed_cmd, '%s', 'delimiter', '\n'); %deprecated - use textscan
% cmds = dataread('file', bed_cmd, '%s', 'delimiter', '\n'); %deprecated - use textscan
% %fid = fopen(bed_cmd);
% %cmds = textscan(fid,'%s','delimiter','\n');
% %fclose(fid);
% threads = {};
% for i = 1: numel(cmds)
%     if mod(i,maxThreads) == 0, fprintf('bedpost %d/%d (%d threads)\n', i, numel(cmds), maxThreads); end;
%     doFslCmd ([cmds{i} ' &'], false);
%     isDone =fullfile(bed_dirX, 'diff_slices', ['data_slice_' sprintf('%4.4d',i-1)],'dyads1.nii.gz'); %file created when thread completes
%     threads = [threads; {isDone}]; %#ok<AGROW>
%     threads = threadSub(threads, maxThreads);
% end
% threadSub(threads, 0);
% pause(2);
% command = fullfile(fslDirSub, 'bin', 'bedpostx_postproc.sh');
% if ~exist(command,'file'), error('You need the special file %s', command); end;
% command=[command ' ', bed_dir];
% doFslCmd (command);
% fprintf ('Bedpost took %f seconds to run.\n', toc(t_start) );
%end doDtiBedpostSub()

function doDtiBedpostSub(dti) %warp atlas to DTI
t_start=tic;
pth = fileparts(dti);
bed_dirX=fullfile(pth, 'bedpost.bedpostX');
%if exist(bed_dirX, 'file'), rmdir(bed_dirX, 's'); end; %666
bed_done=fullfile(bed_dirX, 'xfms', 'eye.mat');
if exist(bed_done,'file'), fprintf('Skipping bedpost (already done)\n'), return, end;
if ~exist(bed_dirX, 'file'), mkdir(bed_dirX); end;
bed_dir=fullfile(pth, 'bedpost');
%if exist(bed_dir, 'file'), rmdir(bed_dir, 's'); end; %666
if ~exist(bed_dir, 'file'), mkdir(bed_dir); end;
dti_u=prepostfixSub('', 'u', dti);
dti_x=fullfile(bed_dir, 'data.nii.gz');
if ~exist(dti_u,'file'), error('Bedpost unable to find %s', dti_u); end;
copyfile(dti_u, dti_x);
[bvec, bval] = getBVec(dti);
dti_x=fullfile(bed_dir, 'bvecs');
copyfile(bvec, dti_x);
dti_x=fullfile(bed_dir, 'bvals');
copyfile(bval, dti_x);
dti_faThr=prepostfixSub('', '_FA_thr', dti);
dti_x=fullfile(bed_dir, 'nodif_brain_mask.nii.gz');
copyfile(dti_faThr, dti_x);
command=sprintf('bedpostx "%s" ', bed_dir);
fslParallelSub;
doFslCmd (command);
while ~exist(bed_done,'file')
    pause(1.0);
end
fprintf ('Bedpost took %f seconds to run.\n', toc(t_start) );
%end doDtiBedpostSub()

function [bvec, bval] = getBVec(dti)
%both.bval
dti = unGzNameSub(deblank(dti));
[p,n] = fileparts(dti);
dti = fullfile(p,n);
bvec = [dti 'both.bvec'];
bval = [dti 'both.bval'];
if exist(bvec,'file') && exist(bval,'file'), return; end;
bvec = [dti '.bvec'];
bval = [dti '.bval'];
if ~exist(bvec,'file') || ~exist(bval,'file'), error('Can not find files %s %s', bvec, bval); end;
%end getBVec()

function doDtiTractSub(dti)
pth = fileparts(dti);
bed_dir=fullfile(pth, 'bedpost');
bed_dirX=fullfile(pth, 'bedpost.bedpostX');
bed_merged=fullfile(bed_dirX, 'merged');
bed_mask=fullfile(bed_dirX, 'nodif_brain_mask');
if ~exist(bed_dir,'file') || ~exist(bed_dirX,'file')  
    fprintf('Please run bedpost to create files %s %s\n',bed_dir, bed_dirX);
    return;
end
template_roiW=prepostfixSub('', '_roi', dti);
dti_u=prepostfixSub('', 'u', dti);
dti_faThr=prepostfixSub('', '_FA_thr', dti);
if ~exist(template_roiW,'file') || ~exist(dti_u,'file') || ~exist(dti_faThr,'file')
    fprintf('Can not find %s or %s or %s\n',template_roiW, dti_u, dti_faThr);
    return;
end
template_roiWThr=prepostfixSub('', '_roi_thr', dti);
command=sprintf('fslmaths "%s" -mas "%s" "%s"',template_roiW, dti_faThr, template_roiWThr);
fprintf('Creating thresholded image %s\n', template_roiWThr);
doFslCmd (command);
fprintf('PROBTRACKX: Create seed data\n');
mask_dir=fullfile(pth, 'masks');
if ~exist(mask_dir, 'file'), mkdir(mask_dir); end;
nROI = 189; %666
%now run probtrackx
prob_dir=fullfile(pth, 'probtrackx');
if ~exist(prob_dir, 'file'), mkdir(prob_dir); end;
nPerm = 5000; %666
t_start=tic;
commands = [];
for i = 1: nROI
    maski=fullfile(mask_dir, [num2str(i),'.nii.gz']);
    command=sprintf('fslmaths "%s" -thr %d -uthr %d -bin "%s" -odt char', template_roiWThr, i,i, maski);
    doFslCmd (command, i == 1); %only show text for 1st region
    command=sprintf('fslstats "%s" -M',  maski);
    [~,cmdout] = doFslCmd (command, i == 1); %only show text for 1st region
    if str2num(cmdout) < 1.0  %#ok<ST2NM>
        delete(maski);
    else
        prob_diri=fullfile(prob_dir, num2str(i));
        if exist(prob_diri, 'file'), rmdir(prob_diri, 's'); end;
        mkdir(prob_diri);
        command=sprintf('probtrackx2 -x "%s" --dir="%s" --forcedir  -P %d -s "%s" -m "%s" --opd --pd -l -c 0.2 --distthresh=0', ...
            maski, prob_diri, nPerm ,bed_merged, bed_mask );
        commands = [commands {command}]; %#ok<AGROW>
    end %if voxels survive
end %for each region
if numel(commands) < 1
   fprintf('No regions survive thresholding with FA (poor normalization?) %s', dti); 
   return;
end
fprintf ('computing probtrackx for %d regions (this may take a while)\n',numel(commands) );   
doThreads(commands, prob_dir);
fprintf ('probtrackx2 took %f seconds to run.\n', toc(t_start) ); %t_start=tic;
%sum(nOK(:))

function doThreads(commands, out_dir)
fslParallelSub;
command_file=fullfile(out_dir, 'commands.txt');
log_dir=fullfile(out_dir, 'logs');
if exist(log_dir, 'file'), rmdir(log_dir, 's'); end;
mkdir(log_dir);
fid = fopen(command_file, 'wt');
for i = 1 : numel(commands)   
    fprintf(fid,'%s\n',commands{i});
end
fclose(fid);
command=sprintf('fsl_sub -l %s -N probtrackx -T 1 -t %s', log_dir, command_file );
doFslCmd (command, false);  
%end doThreads

function fslParallelSub
maxThreads = feature('numCores');
if maxThreads > 15, maxThreads = maxThreads - 1; end;
setenv('FSLPARALLEL', num2str(maxThreads));
%end fslEnvSub()

function doDtiWarpJhuSub(imgs)
if isempty(imgs.T1) || isempty(imgs.DTI), return; end; %required
T1 = prefixSub('wb',imgs.T1); %warped brain extracted image
FA = prepostfixSub('', '_FA', imgs.DTI);
MD = prepostfixSub('', '_MD', imgs.DTI);
if ~exist(T1,'file') || ~exist(FA,'file') || ~exist(MD,'file') , return; end; %required
FA = unGzSub (FA);
MD = unGzSub (MD);
nFA = rescaleSub(FA);
JHU = fullfile(fileparts(which(mfilename)) , 'jhu_spm.nii');
if ~exist(JHU,'file')
    error('Unable to find template %s', JHU);
end
[outhdr, outimg] = nii_reslice_target(JHU, '', T1, 0) ;
roiname = prepostfixSub('', '_roi', imgs.DTI);
if exist(roiname,'file') %e.g. FSL made .nii.gz version
    delete(roiname);
    roiname = prepostfixSub('', '_roi', imgs.DTI);
end
roiname = unGzNameSub(roiname);
outhdr.fname = roiname;
spm_write_vol(outhdr,outimg);
oldNormSub({T1, roiname}, nFA, 8, 10, 0 );
delete(roiname);
wroiname = prepostfixSub('w', '_roi', imgs.DTI);
movefile(wroiname, roiname);
%end doFaMdSub()

function imgname = unGzNameSub(imgname)
[p, n, x] = fileparts(imgname);
if strcmpi(deblank(x),'.gz') %.nii.gz
    imgname = fullfile(p,n);
end
%end unGzNameSub()

function oldNormSub(src, tar, smoref, reg, interp)
%coregister T2 to match T1 image, apply to lesion
if isempty(src) || isempty(tar), return; end;
if ~exist('smoref','var'), smoref = 0; end;
if ~exist('reg','var'), reg = 1; end;
if ~exist('interp','var'), interp = 1; end;
matlabbatch{1}.spm.tools.oldnorm.estwrite.subj.source = src(1);
matlabbatch{1}.spm.tools.oldnorm.estwrite.subj.wtsrc = '';
matlabbatch{1}.spm.tools.oldnorm.estwrite.subj.resample = src;
matlabbatch{1}.spm.tools.oldnorm.estwrite.eoptions.template = {tar};
matlabbatch{1}.spm.tools.oldnorm.estwrite.eoptions.weight = '';
matlabbatch{1}.spm.tools.oldnorm.estwrite.eoptions.smosrc = 8;
matlabbatch{1}.spm.tools.oldnorm.estwrite.eoptions.smoref = smoref; % <-- !!!
matlabbatch{1}.spm.tools.oldnorm.estwrite.eoptions.regtype = 'mni';
matlabbatch{1}.spm.tools.oldnorm.estwrite.eoptions.cutoff = 25;
matlabbatch{1}.spm.tools.oldnorm.estwrite.eoptions.nits = 16;
matlabbatch{1}.spm.tools.oldnorm.estwrite.eoptions.reg = reg;
matlabbatch{1}.spm.tools.oldnorm.estwrite.roptions.preserve = 0;
matlabbatch{1}.spm.tools.oldnorm.estwrite.roptions.bb = [nan nan nan; nan nan nan];%[-78 -112 -70; 78 76 85]; 
matlabbatch{1}.spm.tools.oldnorm.estwrite.roptions.vox = [nan nan nan];%[1 1 1];
matlabbatch{1}.spm.tools.oldnorm.estwrite.roptions.interp = interp;
matlabbatch{1}.spm.tools.oldnorm.estwrite.roptions.wrap = [0 0 0];
matlabbatch{1}.spm.tools.oldnorm.estwrite.roptions.prefix = 'w';
spm_jobman('run',matlabbatch);
%end oldNormSub()


function cleanupDtiDir(imgname)
FA = prepostfixSub('', '_FA', imgname);
MD = prepostfixSub('', '_MD', imgname);
deleteImg(FA);
deleteImg(MD);
[p, n, x] = fileparts(imgname);
bDir = fullfile(p,'bedpost.bedpostX');
if exist(bDir, 'file'), rmdir(bDir, 's'); end;
bDir = fullfile(p,'bedpost');
if exist(bDir, 'file'), rmdir(bDir, 's'); end;
bDir = fullfile(p,'masks');
if exist(bDir, 'file'), rmdir(bDir, 's'); end;
bDir = fullfile(p,'probtrackx');
if exist(bDir, 'file'), rmdir(bDir, 's'); end;
%end rmBedpostDir()

function fsldir = fslDirSub;
fsldir= '/usr/local/fsl/';
%end fslDirSub()

function fslEnvSub
fsldir = fslDirSub;
curpath = getenv('PATH');
k = strfind(curpath, fsldir);
if isempty(k)
	setenv('PATH',sprintf('%s:%s',fullfile(fsldir,'bin'),curpath));
end
%end fslEnvSub()

function command  = fslCmdSub (command)
fsldir = fslDirSub;
cmd=sprintf('sh -c ". %setc/fslconf/fsl.sh; ',fsldir);
command = [cmd command '"'];
%fslCmdSub

function [status,cmdout]  = doFslCmd (command, verbose)
if ~exist('verbose', 'var'), 
    verbose = true;
end;
fslEnvSub;
cmd = fslCmdSub(command);
if verbose
    fprintf('Running \n %s\n', cmd);
end
[status,cmdout]  = system(cmd);
%end doFslCmd()

% function [status,cmdout]  = doFslCmd (command, verbose)
% if ~exist('verbose', 'var'), 
%     verbose = true;
% end;
% fsldir= '/usr/local/fsl/';
% setenv('FSLDIR', fsldir);
% curpath = getenv('PATH');
% setenv('PATH',sprintf('%s:%s',fullfile(fsldir,'bin'),curpath));
% cmd=sprintf('sh -c ". %setc/fslconf/fsl.sh; ',fsldir);
% cmd = [cmd command '"'];
% if verbose
%     fprintf('Running \n %s\n', cmd);
% end
% [status,cmdout]  = system(cmd);
% setenv('PATH',sprintf('%s',curpath));
%end doFslCmd()


% function dtiSub(dtia)
% fsldir= '/usr/local/fsl/';
% if ~exist(fsldir,'dir'), error('%s: fsldir (%s) not found',mfilename,fsldir); end
% flirt = [fsldir 'bin/flirt'];
% if ~exist(flirt,'file')
% 	error('%s: fsl not installed (%s)',mfilename,flirt);
% end
% setenv('FSLDIR', fsldir);
% curpath = getenv('PATH');
% setenv('PATH',sprintf('%s:%s',fullfile(fsldir,'bin'),curpath));
% %./dti_travis.sh "DTIA_LM1001""
% command= [fileparts(which(mfilename)) filesep 'dti_travis2.sh'];
% command=sprintf('sh -c ". %setc/fslconf/fsl.sh; %s "%s""\n',fsldir,command, dtia);
% system(command);
% %end dtiSub()

function n = bvalCountSub(fnm)
[pth,nam] = nii_filepartsSub(fnm);
bnm = fullfile(pth, [nam, '.bval']);
vnm = fullfile(pth, [nam, '.bvec']);
if ~exist(bnm, 'file') || ~exist(vnm, 'file')
    n = 0;
    return;
end
fileID = fopen(bnm,'r');
[A, n] = fscanf(fileID,'%g'); %#ok<ASGLU>
fclose(fileID);
%end bvalCountSub()

function doRestSub(imgs, matName, TRsec, SliceOrder)
if isempty(imgs.T1) || isempty(imgs.Rest), return; end; %we need these images
%if isFieldSub(matName, 'rest_aal'), fprintf('Skipping Rest (already computed) %s\n', imgs.Rest); return; end;
if isFieldSub(matName, 'alf'), fprintf('Skipping Rest (already computed) %s\n', imgs.Rest); return; end;
nii_rest(imgs.Rest, imgs.T1, TRsec, SliceOrder);
nii_nii2mat (prefixSub('fdwa',imgs.Rest), 'rest', matName)
nii_nii2mat (prefixSub('alf_dwa',imgs.Rest), 'alf', matName)
%end doRestSub()

function idx = roiIndexSub(roiName)
[~, ~, idx] = nii_roi_list(roiName);
if idx < 1, error('Invalid roi name %s', roiName); end;
%end roiIndexSub()


function doAslSub(imgs, matName)
if isempty(imgs.T1) || isempty(imgs.ASL), return; end; %we need these images
if isFieldSub(matName, 'cbf'), fprintf('Skipping ASL (CBF already computed) %s\n', imgs.ASL); return; end;
nV = nVolSub (imgs.ASL) ;
[mx, ind] = max(nV);
if mx < 60, fprintf('not enough ASL volumes for %s\n', matName); end;
asl = imgs.ASL(ind,:);
if ~exist(prefixSub('b',imgs.T1),'file'), return; end; %required
if exist(prefixSub('wmeanCBF_0_src',asl),'file'), return; end; %already computed
[cbf, c1L, c1R, c2L, c2R] = nii_pasl12(asl, imgs.T1);
nii_nii2mat(cbf, 'cbf', matName); %2
stat = load(matName);
stat.cbf.nV = nV;
stat.cbf.c1L = c1L;
stat.cbf.c1R = c1R;
stat.cbf.c2L = c2L;
stat.cbf.c2R = c2R;
save(matName,'-struct', 'stat');
%end doAslSub()

function nVol = nVolSub (fnm)
%Report number of volumes
% v= nVols('img.nii');
% v= nVolS(strvcat('ASL_P005.nii', 'ASL_P005_1.nii'))
% v= nVolS({'ASL_P005.nii', 'ASL_P005_1.nii'})
nVol = [];
fnm = cellstr(fnm);
for v = 1 : numel(fnm)
    hdr = spm_vol(deblank(char(fnm(v,:))));
    nVol = [nVol numel(hdr)]; %#ok<AGROW>
end
%end nVol()

function doI3MSub(imgs, matName)
if isempty(imgs.T1), return; end; %we need these images
if isFieldSub(matName, 'i3mT1'), fprintf('Skipping i3m (already computed) %s\n', imgs.ASL); return; end;

w = prefixSub('w',imgs.T1); %warped T1
if  ~exist(w, 'file')  return; end; %exit: we require warped T1
i3m = prefixSub('zw',imgs.T1);
if isFieldSub(matName, 'i3mT1'), return; end;
%if exist(i3m,'file'), return; end; %i3m already computed

nii_i3m(w,'',0.5,10,0.25,1); %i3m T1 image
nii_nii2mat(i3m, 'i3mT1', matName); %4
%end doI3MSub()

function doT1Sub(imgs, matName)
if isempty(imgs.T1), return; end;
if size(imgs.T1,1) > 1 || size(imgs.T2,1) > 1 || size(imgs.Lesion,1) > 1
    error('Require no more than one image for these modalities: T1, T2, lesion');
end;
nii_enat_norm(imgs.T1,imgs.Lesion,imgs.T2);
if isempty(imgs.Lesion), return; end;
if isFieldSub(matName, 'lesion'), return; end;
wr = prefixSub('wr',imgs.Lesion);
if ~exist(wr,'file'), wr = prefixSub('ws',imgs.Lesion); end; %T1 but no T2
if ~exist(wr,'file'), wr = prefixSub('wsr',imgs.Lesion); end; %T1 and T2, smoothed
if ~exist(wr,'file'), error('Unable to find %s', wr); end;
nii_nii2mat(wr, 'lesion', matName); %1
%end doT1Sub()

function bt1 = normSub(t1)
template = fullfile(spm('Dir'),'tpm','TPM.nii');
if ~exist(template,'file')
    error('Unable to find template named %s',template);
end
if iscell(t1), t1 = char(t1); end;
matlabbatch{1}.spm.spatial.preproc.channel(1).vols = {t1};
matlabbatch{1}.spm.spatial.preproc.channel(1).biasreg = 0.001;
matlabbatch{1}.spm.spatial.preproc.channel(1).biasfwhm = 60;
matlabbatch{1}.spm.spatial.preproc.channel(1).write = [0 1];
matlabbatch{1}.spm.spatial.preproc.tissue(1).tpm = {[template ',1']};
matlabbatch{1}.spm.spatial.preproc.tissue(1).ngaus = 1;
matlabbatch{1}.spm.spatial.preproc.tissue(1).native = [1 1];
matlabbatch{1}.spm.spatial.preproc.tissue(1).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(2).tpm = {[template ',2']};
matlabbatch{1}.spm.spatial.preproc.tissue(2).ngaus = 1;
matlabbatch{1}.spm.spatial.preproc.tissue(2).native = [1 1];
matlabbatch{1}.spm.spatial.preproc.tissue(2).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(3).tpm = {[template ',3']};
matlabbatch{1}.spm.spatial.preproc.tissue(3).ngaus = 2;
matlabbatch{1}.spm.spatial.preproc.tissue(3).native = [1 0];
matlabbatch{1}.spm.spatial.preproc.tissue(3).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(4).tpm = {[template ',4']};
matlabbatch{1}.spm.spatial.preproc.tissue(4).ngaus = 3;
matlabbatch{1}.spm.spatial.preproc.tissue(4).native = [1 0];
matlabbatch{1}.spm.spatial.preproc.tissue(4).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(5).tpm = {[template ',5']};
matlabbatch{1}.spm.spatial.preproc.tissue(5).ngaus = 4;
matlabbatch{1}.spm.spatial.preproc.tissue(5).native = [1 0];
matlabbatch{1}.spm.spatial.preproc.tissue(5).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(6).tpm = {[template ',6']};
matlabbatch{1}.spm.spatial.preproc.tissue(6).ngaus = 2;
matlabbatch{1}.spm.spatial.preproc.tissue(6).native = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(6).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.warp.mrf = 1;
matlabbatch{1}.spm.spatial.preproc.warp.cleanup = 1;
matlabbatch{1}.spm.spatial.preproc.warp.reg = [0 0.001 0.5 0.05 0.2];
matlabbatch{1}.spm.spatial.preproc.warp.affreg = 'mni';
matlabbatch{1}.spm.spatial.preproc.warp.fwhm = 0;
matlabbatch{1}.spm.spatial.preproc.warp.samp = 3;
matlabbatch{1}.spm.spatial.preproc.warp.write = [1 1];
spm_jobman('run',matlabbatch);
bt1 = extractSub(0.005, t1, prefixSub('c1', t1), prefixSub('c2', t1));
%end normSub();

function t1Bet = extractSub(thresh, t1, c1, c2, c3)   
%subroutine to extract brain from surrounding scalp
% t1: anatomical scan to be extracted
% c1: gray matter map
% c2: white matter map
% c3: [optional] spinal fluid map
fprintf('Brain extraction of %s\n', t1);
[pth,nam,ext] = spm_fileparts(t1);
%load headers
mi = spm_vol(t1);%bias corrected T1
gi = spm_vol(c1);%Gray Matter map
wi = spm_vol(c2);%White Matter map
%load images
m = spm_read_vols(mi);
g = spm_read_vols(gi);
w = spm_read_vols(wi);
if nargin > 4 && ~isempty(c3)
   ci = spm_vol(c3);%CSF map
   c = spm_read_vols(ci);
   w = c+w; 
end;
w = g+w;
if thresh <= 0
    m=m.*w;
else
    mask= zeros(size(m));
    for px=1:length(w(:)),
      if w(px) >= thresh
        mask(px) = 255;
      end;
    end;
    spm_smooth(mask,mask,1); %feather the edges
    mask = mask / 255;
    m=m.*mask;
end;
mi.fname = fullfile(pth,['b',  nam, ext]);
mi.dt(1) = 4; %16-bit precision more than sufficient uint8=2; int16=4; int32=8; float32=16; float64=64
spm_write_vol(mi,m);
t1Bet = mi.fname;
%end extractSub()


function is = isFieldSub(matname, fieldname)
is = false;
if ~exist(matname, 'file'), return; end;
m = load(matname);
is = isfield(m, fieldname);
%end isFieldSub

function nam = prefixSub (pre, nam)
[p, n, x] = spm_fileparts(nam);
nam = fullfile(p, [pre, n, x]);
%end prefixSub()

function imgs = unGzAllSub(imgs)
imgs.ASL = unGzSub(imgs.ASL);
%imgs.DTI = unGzSub(imgs.DTI); %fsl is fine with gz
%imgs.DTIrev = unGzSub(imgs.DTIrev); %fsl is fine with gz
imgs.fMRI = unGzSub(imgs.fMRI);
imgs.Lesion = unGzSub(imgs.Lesion);
imgs.Rest = unGzSub(imgs.Rest);
imgs.T1 = unGzSub(imgs.T1);
imgs.T2 = unGzSub(imgs.T2);
%end subjStructSub()

function fnm = unGzSub (fnm)
if isempty(fnm), return; end;
innames = fnm;
fnm = [];
for i = 1: size(innames,1)
    f = deblank(innames(i,:));
    f = unGzCSub(f);
    fnm = strvcat(fnm, f);
end
%end unGzSub()

function deleteImg(fnm)
if isempty(fnm), return; end;
fnm = deblank(fnm);
[pth,nam,ext] = spm_fileparts(fnm);
if strcmpi(ext,'.gz') %.nii.gz
    fnm = fullfile(pth, nam);
    [pth,nam,ext] = spm_fileparts(fnm);
end
if strcmpi(ext,'.nii') %.nii.gz
    fnm = fullfile(pth, nam);
end
fnmN = [fnm, '.nii'];
if exist(fnmN,'file')
  delete(fnmN);  
end
fnmZ = [fnm, '.nii.gz'];
if exist(fnmZ,'file')
  delete(fnmZ);  
end
%
function [fnm, isGz] = unGzCSub (fnm)
if isempty(fnm), return; end;
fnm = deblank(fnm);
[pth,nam,ext] = spm_fileparts(fnm);
isGz = false;
if strcmpi(ext,'.nii') % fsl can not handle .nii coexisting with .nii.gz
    fnmz = fullfile(pth, [nam ext '.gz']);
    if exist(fnmz, 'file') && exist(fnm, 'file')
        delete(fnmz); %remove .nii.gz if .nii exists
    end
end
if strcmpi(ext,'.gz') %.nii.gz
    ofnm = fnm;
    fnm = char(gunzip(fnm));  
    isGz = true;
    delete(ofnm);
elseif strcmpi(ext,'.voi') %.voi -> 
    onam = char(gunzip(fnm));
    fnm = fullfile(pth, [nam '.nii']);
    movefile(onam,fnm);
    isGz = true;
end;  
%end unGzSub()


function imgs = subjStructSub(subjName)
imgs.name = subjName;
imgs.ASL = '';
imgs.DTI = '';
imgs.DTIrev = '';
imgs.fMRI = '';
imgs.Lesion = '';
imgs.Rest = '';
imgs.T1 = '';
imgs.T2 = '';
%end subjStructSub()

function imgName = imgfindSub (imgName, imgKey, inDir)
%look for a filename that includes imgKey in folder inDir or subfolders
% for example if imgKey is 'T1' then T1 must be in both folder and image name myFolder\T1\T1.nii
%if ~isempty(imgName), return; end;
[pth, nam] = fileparts(inDir); %#ok<ASGLU> %e.g. 'T1folder' for /test/T1folder
nameFiles = subFileSub(inDir);
nameFiles = sort(nameFiles); %take first file for multiimage sequences, e.g. ASL
for i=1:size(nameFiles,1)
    pos = isStringInKey (nameFiles(i), imgKey);
    if pos == 1 && isImgSub(char(nameFiles(i)))
        imgName = strvcat(imgName, [inDir, filesep, char(nameFiles(i))]);
        
    end; %do not worry about bvec/bval
end
if isempty(imgName), fprintf('WARNING: unable to find any "%s" images in folder %s\n',deblank(imgKey(1,:)), inDir); end;
if size(imgName) > 1, imgName = deblank(imgName(1,:)); end;
%end imgfindSub()

function isKey = isStringInKey (str, imgKey)
isKey = true;
for k = 1 : size(imgKey,1)
    key = deblank(imgKey(k,:));
    pos = strfind(lower(char(str)),lower(key));
    if ~isempty(pos), isKey = pos(1); return; end;
end
isKey = false;
%isStringInKey()
    
% function imgName = imgfindSub (imgName, imgKey, inDir, imgKey2)
% %look for a filename that includes imgKey in folder inDir or subfolders
% % for example if imgKey is 'T1' then T1 must be in both folder and image name myFolder\T1\T1.nii
% if ~exist('imgKey2','var'), imgKey2 = imgKey; end;
% if ~isempty(imgName), return; end;
% [pth, nam] = fileparts(inDir); %#ok<ASGLU> %e.g. 'T1folder' for /test/T1folder 
% if isempty(strfind(lower(char(nam)), lower(imgKey))), return; end;
% if exist([inDir,filesep, 'Native'],'file')
%     inDir = [inDir,filesep, 'Native'];
%     %fprintf('xxxx %s\n', inDir);
% end
% nameFiles = subFileSub(inDir);
% nameFiles = sort(nameFiles); %take first file for multiimage sequences, e.g. ASL
% for i=1:size(nameFiles,1)
%     pos = strfind(lower(char(nameFiles(i))),lower(imgKey));
%     if isempty(pos)
%         pos = strfind(lower(char(nameFiles(i))),lower(imgKey2));
%     end
%     if ~isempty(pos) && isImgSub(char(nameFiles(i)))
%         imgName = [inDir, filesep, char(nameFiles(i))];
%         return
%     end; %do not worry about bvec/bval
% end
% fprintf('WARNING: unable to find any "%s" images in folder %s\n',imgKey, inDir);
% %end imgfindSub

function nameFiles=subFileSub(pathFolder)
d = dir(pathFolder);
isub = ~[d(:).isdir];
nameFiles = {d(isub).name}';
%end subFileSub()

function nameFolds=subFolderSub(pathFolder)
d = dir(pathFolder);
isub = [d(:).isdir];
nameFolds = {d(isub).name}';
nameFolds(ismember(nameFolds,{'.','..'})) = [];
%end subFolderSub()

function isImg = isImgSub (fnm)
[pth,nam,ext] = spm_fileparts(fnm); %#ok<ASGLU>
isImg = false;
if strcmpi(ext,'.gz') || strcmpi(ext,'.voi') || strcmpi(ext,'.hdr') || strcmpi(ext,'.nii')
    isImg = true;
end;  
%end isImgSub()

function [pth,nam,ext,num] = nii_filepartsSub(fname)
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
%end nii_filepartsSub()