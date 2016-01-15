function mat2pdf (matName)
%Show fidelity test for NiiStat mat files
% matName(s): files to analyze
%Example
%mat2pdf('P017.mat')
if isempty(which('spm')), error('Install SPM'); end;
if ~exist('matName','var')  %not specified
   [A,Apth] = uigetfile({'*.mat'},'Select first NiiStat file');
   matName = [Apth, A];
end;
mat = load(matName);
printSub(mat, matName);
%

function printSub(mat, matName)
spm_clf;
spm_figure('Clear', 'Graphics');
spm_orthviews('Reset');
fig = spm_figure('FindWin','Graphics');
ax  = axes('Position',[0.1 0.5 0.35 0.3],'Visible','off','Parent',fig);
[~, nam] = fileparts(matName);
text(0.0,-1.4, sprintf('normalized %s',nam),'Parent',ax, 'FontSize',24, 'fontn','Arial');
i3 = mat2niiSub(mat, 'i3mT1', [0.51 0.01 .4 .49]);
fm = mat2niiSub(mat, 'fmri', [0.01 0.32 .4 .49]);
md = mat2niiSub(mat, 'md',[0.51 0.32 .4 .49]);
rs = mat2niiSub(mat, 'alf',[0.01 0.64 .4 .49]);%%%
cb = mat2niiSub(mat, 'cbf',[0.51 0.64 .4 .49]);
ls = mat2niiSub(mat, 'lesion',[0.01 0.01 .4 .49]);
if ~isempty(ls)
     XYZmm = getCenterOfIntensitySub(ls);
     spm_orthviews('setcoords',XYZmm);
end;
spm_print(matName);
spm_clf;
rSub(i3); rSub(fm); rSub(md); rSub(rs); rSub(cb); rSub(ls);
%end printSub()

function rSub(fnm)
if isempty(fnm) || ~exist(fnm,'file'), return; end;
delete(fnm);
%end rSub()

function fnm = mat2niiSub(mat, modality, pos)
% e.g. mat2niiSub(mat, 'i3mT1', [0.51 0.01 .4 .49]
fnm = [];
if ~isfield(mat,modality) || ~isfield(mat.(modality),'hdr'), return; end; 
hdr = mat.(modality).hdr;
img = mat.(modality).dat;
fnm = [ 'temp_' modality '.nii'];
hdr.fname = fnm;
spm_write_vol(hdr,img);
spm_orthviews('Image',fnm,pos);
fig = spm_figure('FindWin','Graphics');
ax  = axes('Position',[0.0 0.1 1.0 1.0],'Visible','off','Parent',fig);
text(pos(1)+0.2,pos(2)+0.1, modality,'Parent',ax, 'FontSize',24, 'fontn','Arial', 'color','red');
%end mat2niiSub()

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
