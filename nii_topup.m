function nii_topup(Vp,Vn, acqTimeSec, flipDim, viewResults)
%wrapper for DTI preprocessing using TOPUP - requires Matlab and FSL
%  Vp: 4D DTI data set(s) with positive blip
%  Vn: 4D DTI data set(s) with negative blip
%  acqTimeSec: time in seconds to acquire K-space (echoSpacing*Lines/iPAT)
%  flipDim: Dimension where polarity is reversed: 1=X(L/R), 2=Y(A/P), 3=Z(I/S)
%  viewResults: if true, FSLView is launched to display results
%Assumes b-value and b-vector files share name of Vp and Vn
% http://fsl.fmrib.ox.ac.uk/fsl/fslwiki/TOPUP/TopupUsersGuide
%
%Example
% nii_topup('a.nii.gz','b.nii.gz',0.03388,2); %assumes a.bval, a.bvec, b.bval, b.bvec
% nii_topup('a.nii','b.nii',0.03388,2); 
% nii_topup('aa.nii', '');
% nii_topup; %GUI

if nargin<1, 
   [files,pth] = uigetfile({'*.gz;*.nii;*.hdr;';'*.*'},'Choose positive DTI[s]', 'MultiSelect', 'on');
   Vp = strcat(pth,char(files));
end 
if isempty(Vp), disp('Error: no images provided'); return;  end; 
if nargin<2, 
    [files,pth] = uigetfile({'*.gz;*.nii;*.hdr;';'*.*'},'Choose negative DTI[s] (SAME ORDER)', 'MultiSelect', 'on'); 
    if ~ischar(files)
        Vn = {};
    else
        Vn = strcat(pth,char(files));
    end
end
if ~isempty(Vn) 
    if nargin<4 %get defaults for topup
        acqTimeSec = 0.03388; 
        answer = inputdlg({'Readout Time (sec):','Flipped dimension (1=X, 2=Y, 3=Z)'},'Input',1,{num2str(acqTimeSec), '2'});
        if isempty(answer), return; end;
        acqTimeSec = cellfun(@str2num, answer(1))
        flipDim = cellfun(@str2num, answer(2))
        if ((flipDim < 1) || (flipDim > 3)), disp('Error: flipDim must be a number 1..3.'); return; end;
        if ((acqTimeSec < 0) || (acqTimeSec > 0.5)), fprintf('Warning: acqTimeSec of %f seems implausible.',acqTimeSec); end;
    end
    fprintf('Assuming %f sec readout time\n',acqTimeSec);
    fprintf('Assuming polarity was reversed in dimension %d (1=X(L/R), 2=Y(A/P), 3=Z(I/S)\n',flipDim);
end
if ~exist('viewResults', 'var'), viewResults = true; end;
%check inputs
if ~isempty(Vn) && (size(Vp,1) ~= size(Vn,1)), disp('Error: requires the same number of positive and negative images (pairs).'); return; end;
setFslDirSub; %check FSL installation
for i=1:size(Vp,1) %process each image
    tStart = tic;
    nameP = deblank(Vp(i,:));
    if isempty(Vn)
       [imgEcc, imgMask] = eddySub(nameP); %only one phase encoding direction - cannot run topup
    else
       [imgEcc, imgMask] = topupSub(nameP, deblank(Vn(i,:)),acqTimeSec,flipDim); %run topup and eddy
    end
    fdtSub(imgEcc, imgMask);
    fprintf('Created tensor images %s in %f seconds.\n',nameP, toc(tStart) );
end;
if ~viewResults, return; end;
viewSub(imgEcc);
%end nii_topup()

function [imgEcc, imgMask] =  eddySub(imgName) 
%input: DTI data
%output: undistorted DTI (imgEcc) and brain mask (imgMask)
[pth nam ext] = fsl_filepartsSub(imgName);
imgEcc = fullfile(pth,['e' nam]); %name of eddy current
[bv0P] = getB0vols (imgName, imgEcc); %get b-values for 1st series
if isempty(bv0P), error('Unable to determine reference (B0) volume'); end;
fslCmd(sprintf('eddy_correct %s %s %d',imgName, imgEcc, bv0P(1)-1 )); %-1 since FSL volumes indexed from zero 
imgMask = fullfile(pth,['b' nam ]); %name of brain mask - fsl appends '_mask'
fslCmd(sprintf('bet %s %s  -f 0.2 -R -g 0 -n -m',imgName, imgMask));
imgMask = fullfile(pth,['b' nam '_mask']); %name of brain mask - fsl appends '_mask'
%end preprocSub

function [outNameV] = fdtSub(imgName, imgMask)
%process DTI data using the FMRIB's Diffusion Toolbox http://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FDT
[pth nam ext] = fsl_filepartsSub(imgName);
nameVal = fullfile(pth,[nam '.bval']); %name for b-values
nameVec = fullfile(pth,[nam  '.bvec']); %name for b-vectors
outNameV = fullfile(pth,[nam]); %name of vectors
fslCmd(sprintf('dtifit --data=%s --out=%s --mask=%s --bvecs=%s --bvals=%s',imgName,outNameV, imgMask, nameVec, nameVal));
%end fdtSub()

function viewSub(img) %display fractional anisotropy and vector images
[pth,nam] = fileparts(img);
if isempty(pth), pth = pwd; end;
faNam = fullfile(pth, [nam '_FA.nii.gz']); %Fractional Anisotropy
v1Nam = fullfile(pth, [nam '_V1.nii.gz']); %Primary vector direction
fslCmd(sprintf('fslview %s %s &',faNam,v1Nam));
%end dtiSub()

function nVol = fslInfoSub (nam)
nVol = 0;
[status,result] = fslCmd (sprintf('fslinfo %s',nam ));
if status ~= 0, return; end; %error reading image
t = textscan(result, '%s %s');
x = key2num(t, 'dim1');
y = key2num(t, 'dim2');
z = key2num(t, 'dim3');
nVol = key2num(t, 'dim4');
if mod(x,2) || mod(y,2) || mod(z,2)
   fprintf('Volume dimensions odd (%dx%dx%d), topup may fail due to "Subsampling levels incompatible with image data"\n', x,y,z);
end
%end fslInfoSub()

function val = key2num(t, key)
%check first column of cell array for key - if it exists return value in second column of the same row.
d = strfind(t{1}, key);
i = find(not(cellfun('isempty', d)));
val = nan;
if ~isempty(i)
   val = str2double(t{2}{i(1)}); 
end
%end key2num()

function setFslDirSub
%ensure fsl is installed and recent version (with topup)
fsldir = fslDirSub;
setenv('FSLDIR',fsldir);  % this to tell where FSL folder is
if ~exist(fsldir,'dir')
	error('%s: error fsldir (%s) not found',mfilename, fsldir);
end
setenv('PATH',sprintf('%s:%s',fullfile(fsldir,'bin'),getenv('PATH')));
topup = [fsldir '/bin/topup'];
if ~exist(topup)
	error('%s: error %s not found',mfilename,topup);
end
%end setFslDirSub()

function fsldir = fslDirSub
%return path to FSL
fsldir= '/usr/local/fsl/';
%end fslDirSub()

function [dti_u, dti_b] = topupSub (vol1,vol2,acqTimeSec,flipDim)
%input: two DTI images (vol1, vol2) with reverse phase encoding
%output: undistorted (dti_u) and brain-mask (dti_b)
nVol = fslInfoSub(vol1);
bv0P = getB0vols (vol1); %get b-values for 1st series
bv0N = getB0vols (vol2); %get b-values for 2nd series
if (bv0N ~= bv0P), fprintf('Error: different b-values: positive and negative images must be pairs %s %s \n',vol1,vol2); return; end;
both_b0 = spliceVolSub (vol1, vol2, bv0P); %make volume with just b0 images
nameAcqParams = writeAcqParamsSub(vol1,bv0P,acqTimeSec,flipDim);
[pth nam] = fsl_filepartsSub(vol1);
dti_u= fullfile(pth,[nam 'u']);
dti_t = fullfile(pth,['tp' nam]);
dti_b0s = fullfile(pth,['b0' nam]);
fslCmd(sprintf('topup --imain=%s --datain=%s --config=b02b0.cnf --out=%s --iout=%s',both_b0,nameAcqParams,dti_t, dti_b0s));
dti_tb0 = fullfile(pth,['mtp' nam]);
fslCmd(sprintf('fslmaths %s -Tmean %s', dti_b0s, dti_tb0));
dti_bx = fullfile(pth,['b' nam]);
fslCmd(sprintf('bet %s %s  -f 0.2 -R -n -m', dti_tb0, dti_bx));
dti_b = fullfile(pth,['b' nam '_mask']);
dti_txt2= fullfile(pth,[nam '_index.txt']); 
dlmwrite(dti_txt2,[ones(nVol,1)' 2*ones(nVol,1)'],'delimiter','\t');
dti_merge=fullfile(pth,[nam 'both']); 
fslCmd(sprintf('fslmerge -t %s %s %s', dti_merge, vol1, vol2));
[dti_bvecm, dti_bvalm] = bValVecMergeSub(vol1,vol2);
fslCmd(sprintf('eddy --imain=%s --mask=%s --acqp=%s --index=%s --bvecs=%s --bvals=%s --topup=%s --out=%s', ...
    dti_merge, dti_b, nameAcqParams, dti_txt2, dti_bvecm, dti_bvalm, dti_t, dti_u));
%end topupSub()  
  
function [nameAcqParams] = writeAcqParamsSub (vol1, vols, acqTimeSec, flipDim);
%create topup format text file that describes flipped dimension and acquisition time
[pth nam ext] = fsl_filepartsSub(vol1);
nameAcqParams = fullfile(pth,['b0' nam '.txt']);
M = zeros(length(vols)*2,4);
M(:,flipDim)=-1;
M(1:numel(vols),flipDim)=1;
M(:,4)=acqTimeSec;
fid=fopen(nameAcqParams,'wt');
dlmwrite(nameAcqParams, M, 'delimiter', '\t','precision', 6);
fclose(fid);
%end writeAcqParamsSub()

function [spliceName] = spliceVolSub (vol1, vol2, vols);
%splice B0 images from vol1 and vol2 into new volume
% example spliceVolSub ('a.nii', 'b.nii', [0,3]); creates volume with 1st and 3rd volumes from each a.nii and b.nii
tmpNames = '';
for vol=1:length(vols)
    outName = sprintf('data_tmp%d',vol);
    tmpNames = [tmpNames ' ' outName]; 
    fslCmd(sprintf('fslroi %s %s %d %d',vol1,outName,(vols(vol)-1),1 ) );
end;
for vol=1:length(vols)
    outName = sprintf('data_tmp%d',vol+length(vols));
    tmpNames = [tmpNames ' ' outName];
    fslCmd(sprintf('fslroi %s %s %d %d',vol2,outName,(vols(vol)-1),1 ) );
end;
[pth nam ext] = fsl_filepartsSub(vol1);
spliceName   = fullfile(pth,['b0' nam ext]);
fslCmd(sprintf('fslmerge -t %s %s',spliceName,tmpNames ));
for vol=1:(2*length(vols))
    delete(sprintf('data_tmp%d.nii.gz',vol));
end
%end spliceVolSub()

function [status, result] = fslCmd (Cmd)
%execute a fsl command, e.g. fslCmd('fslinfo a.nii');
command=sprintf('sh -c ". %setc/fslconf/fsl.sh; %sbin/%s"\n',fslDirSub, fslDirSub, Cmd);
fprintf('%s\n', command);
[status,result] = system(command);
%end fslCmd()

function [bv0] = getB0vols (imgName, copyName)
%read .bval file and return indices for B0 volumes
[pth nam] = fsl_filepartsSub(imgName);
nameVal = fullfile(pth,[nam '.bval']); %name for b-values
nameVec = fullfile(pth,[nam  '.bvec']); %name for b-vectors  
if (exist(imgName, 'file') == 0) , fprintf('Unable to find required image %s\n',imgName); return; end;  
if ( (exist(nameVal, 'file') == 0) || (exist(nameVec, 'file') == 0) ), fprintf('Unable to find required DTI files %s and %s\n',nameVec,nameVal); return; end;  
%bv = importdata(nameVal); %<- does not work with Matlab 2014b on Linux
fileID = fopen(nameVal);
bv = cell2mat( textscan(fileID,'%d'));
fclose(fileID);
bv0 = find(bv == 0);
if isempty(bv0), error('Unable to find any zero b-values in %s\n',nameVal); end;
if ~exist('copyName', 'var'), return; end;
%next - copy b-values and vectors to new name
[pth nam] = fsl_filepartsSub(copyName);
nameVal2 = fullfile(pth,[nam '.bval']); %name for b-values
nameVec2 = fullfile(pth,[nam  '.bvec']); %name for b-vectors 
copyfile(nameVal, nameVal2);
copyfile(nameVec, nameVec2);
%end getB0vols()

function [dti_bvecm, dti_bvalm] = bValVecMergeSub(vol1,vol2)
%combine DTI information for two volumes
[pth nam] = fsl_filepartsSub(vol1);
[pth2 nam2] = fsl_filepartsSub(vol2);
vnam = fullfile(pth,[nam '.bval']); %name for b-values
vnam2 = fullfile(pth2,[nam2  '.bval']); %name for b-vectors 
dti_bvalm = fullfile(pth,[nam 'u.bval']); %name for b-values
textMergeSub(vnam, vnam2, dti_bvalm);
vnam = fullfile(pth,[nam '.bvec']); %name for b-values
vnam2 = fullfile(pth2,[nam2  '.bvec']); %name for b-vectors 
dti_bvecm = fullfile(pth,[nam 'u.bvec']); %name for b-values
textMergeSub(vnam, vnam2, dti_bvecm, 3);
%end bValVecMerge()

function textMergeSub(nam1,nam2,outnam, nRow)
%combine two text files
fileID = fopen(nam1);
b1 = cell2mat( textscan(fileID,'%f'));
fclose(fileID);
%read 2nd file
fileID = fopen(nam2);
b2 = cell2mat( textscan(fileID,'%f'));
fclose(fileID);
if exist('nRow', 'var')
    b1 = reshape(b1, [numel(b1)/nRow nRow]);
    b2 = reshape(b2, [numel(b2)/nRow nRow]);
end
dlmwrite(outnam,[b1' b2'],'delimiter','\t');
%end textMergeSub

function [pth nam ext] = fsl_filepartsSub(fileName) 
% a.nii.gz has the extension ".nii.gz" not ".nii"
    [pth nam ext] = fileparts(fileName);
    if (length(ext)==3)  && min((ext=='.gz')==1)
       [pth nam ext2] = fileparts( fullfile(pth, nam)); %remove .nii for .nii.gz
       ext = [ext2 ext];
    end;
%end fsl_filepartsSub()
