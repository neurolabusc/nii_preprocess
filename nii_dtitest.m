function nii_dtitest (dtiBvec)
%test all possible vectors and polarities
%assumes angulations have been correctly adjusted
% dtiNii: name of bvec file(s), e.g. img.bvec
%Examples
% nii_dtibatch
% nii_dtibatch('DTI_axial_04.bvec');
% nii_dtibatch(strvcat('DTI_axial_04.bvec','DTI_30_deg_AP_05.bvec'));
%Chris Rorden 4 July 2014, BSD License
%Subsequently you can view the processed images from the command line
% fslview DTI_axial_04_FA.nii.gz DTI_axial_04_V1.nii.gz

if ~exist('dtiBvec','var') %file not specified
   [A,Apth] = uigetfile({'*.bvec';'*.*'},'Select b-vector file(s)', 'MultiSelect', 'on');
   dtiBvecNames = strcat(Apth,char(A));
end;
fsldir= '/usr/local/fsl';
for i=1:size(dtiBvecNames,1)
    dtiBvec = deblank(dtiBvecNames(i,:)); %positive image 
    [pth,nam,ext] = fileparts(dtiBvec); %#ok<NASGU>
    imgNam = fullfile(pth, [nam '.nii']); %file.nii
    if ~exist(imgNam,'file'), imgNam = fullfile(pth, [nam '.nii.gz']); end; %file.nii.gz
    exist(imgNam,'file');
    preprocSub(fsldir,imgNam,0);
    %next: permute all possible b-vector alternatives
    dtiSub(fsldir,dtiBvec)
end
viewSub(fsldir,dtiBvec) %compute vectors
%end main loop.... subroutines follow

function preprocSub(fsldir,imgNam,refVol) %eddy current correct and brain extract
setenv('FSLDIR', fsldir);
[pth,nam,ext] = fileparts(imgNam); %#ok<NASGU,ASGLU>
maskNam = fullfile(pth, 'dti'); %will generate image "dti_mask.nii.gz"
command=sprintf('sh -c ". ${FSLDIR}/etc/fslconf/fsl.sh; ${FSLDIR}/bin/bet %s %s -f 0.3 -g 0 -n -m"\n',imgNam,maskNam);
system(command);
eccNam = fullfile(pth, 'dti_eddy'); %will generate image "dti_eddy.nii.gz"
setenv('FSLDIR', fsldir);
setenv('PATH', [getenv('PATH') ':/usr/local/fsl/bin'])
command=sprintf('sh -c ". ${FSLDIR}/etc/fslconf/fsl.sh; ${FSLDIR}/bin/eddy_correct %s %s %d"\n',imgNam,eccNam, refVol);
system(command);
%end preprocSub

function dtiSub(fsldir,vNam) %compute vectors
%%/usr/local/fsl/bin/dtifit --data=/Users/rorden/Desktop/sliceOrder/dicom2/dtitest/dti_eddy.nii.gz --out=/Users/rorden/Desktop/sliceOrder/dicom2/dtitest/dti --mask=/Users/rorden/Desktop/sliceOrder/dicom2/dtitest/dti_mask.nii.gz --bvecs=/Users/rorden/Desktop/sliceOrder/dicom2/dtitest/s004a001.bvec --bvals=/Users/rorden/Desktop/sliceOrder/dicom2/dtitest/s003a001.bval
[pth,nam,ext] = fileparts(vNam); %#ok<NASGU>
bNam = fullfile(pth, [ nam '.bval'] ); %Eddy corrected data
eccNam = fullfile(pth, 'dti_eddy.nii.gz'); %Eddy corrected data
maskNam = fullfile(pth, 'dti_mask.nii.gz'); %Eddy corrected data
[pth,nam,ext] = fileparts(vNam); %#ok<NASGU>
outNam = fullfile(pth, nam); %Eddy corrected data
setenv('FSLDIR', fsldir);
setenv('PATH', [getenv('PATH') ':/usr/local/fsl/bin'])
command=sprintf('sh -c ". ${FSLDIR}/etc/fslconf/fsl.sh; ${FSLDIR}/bin/dtifit --data=%s --out=%s --mask=%s --bvecs=%s --bvals=%s"\n',eccNam, outNam, maskNam,vNam,bNam);
system(command);
delete(fullfile(pth, [nam '_V2.nii.gz']));
delete(fullfile(pth, [nam '_V3.nii.gz']));
delete(fullfile(pth, [nam '_L1.nii.gz']));
delete(fullfile(pth, [nam '_L2.nii.gz']));
delete(fullfile(pth, [nam '_L3.nii.gz']));
delete(fullfile(pth, [nam '_MO.nii.gz']));
delete(fullfile(pth, [nam '_MD.nii.gz']));
delete(fullfile(pth, [nam '_S0.nii.gz']));
%end dtiSub

function viewSub(fsldir,vNam) %compute vectors
%%/usr/local/fsl/bin/dtifit --data=/Users/rorden/Desktop/sliceOrder/dicom2/dtitest/dti_eddy.nii.gz --out=/Users/rorden/Desktop/sliceOrder/dicom2/dtitest/dti --mask=/Users/rorden/Desktop/sliceOrder/dicom2/dtitest/dti_mask.nii.gz --bvecs=/Users/rorden/Desktop/sliceOrder/dicom2/dtitest/s004a001.bvec --bvals=/Users/rorden/Desktop/sliceOrder/dicom2/dtitest/s003a001.bval
[pth,nam,ext] = fileparts(vNam); %#ok<NASGU>
faNam = fullfile(pth, [nam '_FA.nii.gz']); %Eddy corrected data
v1Nam = fullfile(pth, [nam '_V1.nii.gz']); %Eddy corrected data
setenv('FSLDIR', fsldir);
setenv('PATH', [getenv('PATH') ':/usr/local/fsl/bin'])
command=sprintf('sh -c ". ${FSLDIR}/etc/fslconf/fsl.sh; ${FSLDIR}/bin/fslview %s %s &"\n',faNam,v1Nam);
system(command);
%end dtiSub