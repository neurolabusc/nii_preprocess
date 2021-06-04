function dtinam = nii_dwidenoise (dtinam, degibbs)
% dtinam : name of dti image to denoise
% isDeGibbs : (optional, default=false) 0=no, 1=mrdegibbs, 2=unring
%             n.b. do not use if partial-fourier used by acquisition
%Examples
% nii_dwidenoise(('dti.nii')
% nii_dwidenoise(('dti.nii', 1) %use mrtrix to remove Gibbs artifacts
% nii_dwidenoise(('dti.nii', 1) %use mrtrix to remove Gibbs artifacts

if ~exist('dtinam','var')  %fnmFA not specified
   fprintf('Select DTI image\n');
   [A,Apth] = uigetfile({'*.nii;*.gz;*.hdr;';'*.*'},'Select DTI image');
   if isnumeric(A), return; end;
   dtinam = [Apth, A];
end
start = tic;
dtinam = dwidenoiseSub(dtinam);
if exist('degibbs','var') && (degibbs > 0)
    if (degibbs > 1)
        unringSub(dtinam);
    else
        degibbsSub(dtinam);
    end
end
fprintf('denoise/gibbs required %g seconds\n', toc(start) );
%end nii_dwidenoise()

function unringSub(fname)
%https://bitbucket.org/reisert/unring
% unring - tool for removal of the Gibbs ringing artefact
% Usage: outvol = unring(invol,params)
% Options: invol - input volume 
%          params - 3x1 array with [minW maxW nsh]
%                     nsh discretization of subpixel spaceing (default 20)
%                     minW  left border of window used for TV computation (default 1)
%                     maxW  right border of window used for TV computation (default 3)
%You may need to compile ringRm mex file for your computer ~
%   mex -I/usr/local/include -L/usr/local/lib -lfftw3 ringRm.cpp
% based on the algorithm in the publication 
% Kellner, E, Dhital B., Kiselev VG and Reisert, M. 
% Gibbs?ringing artifact removal based on local subvoxel?shifts. 
% Magnetic resonance in medicine, 76(5), 1574-1581.
if ~exist(fname, 'file'), warning( 'Unable to find %s', fname); return; end;
[hdr, img] = read_volsSub (fname);
fprintf('unringing %d volumes (make sure this acquisition did not use partial fourier)\n', size(img,4));
params = [1 3 20];
img = ringRm(double(img),params);
hdr = hdr(1);
%note we will overwrite as denoise created modified image
% This could be adpated, but make sure to handle bvec, bval and nii.gz
isDeleteTemp = true;
if isDeleteTemp
    delete(fname); %remove 
else
    [pth,nam,ext] = filepartsSub(fname);
    tempname = fullfile(pth,['denoise_', nam, ext ]);
	copyfile(fname, tempname);    
end
fname = hdr.fname; %handle .nii -> .nii.gz
for v = 1 : size(img,4) %smooth each volume
    hdr.n(1)=v;
    spm_write_vol(hdr,img(:,:,:,v));
end
%end unringSub()

function [hdr, img] = read_volsSub (fnm)
[fnm, isGz] = unGzSub (fnm); %convert FSL .nii.gz to .nii
hdr = spm_vol(fnm); %load header data
img = spm_read_vols(hdr); %load image data
if (isGz), delete(fnm); end; %remove .nii if we have .nii.gz
%end read_volsSub()

function [fnm, isGz] = unGzSub (fnm)
[pth,nam,ext] = spm_fileparts(fnm); %#ok<ASGLU>
isGz = false;
if strcmpi(ext,'.gz') %.nii.gz
    fnm = char(gunzip(fnm));
    isGz = true;
end;
%end unGzSub()

function degibbsSub(fname)
%http://mrtrix.readthedocs.io/en/latest/reference/commands/mrdegibbs.html
if ~exist(fname, 'file'), warning('Unable to find %s', fname); return; end;
exenam = fullfile(nii_mrtrix_pth, 'mrdegibbs');
if ~exist(exenam, 'file')
    error('unable to find %s', exenam); 
end
[pth,nam,ext] = filepartsSub(fname);
tempname = fullfile(pth,['temp', nam, ext]);
copyfile(fname, tempname);
cmd = '';
%optional: datatype - otherwise defaults to 32-bit float
if strcmpi(ext,'.nii') %ignore for .nii.gz files (alternatively ungzip to read header)
    hdr = spm_vol([tempname,',1']);
    if hdr.dt(1) == 4 %int16
        cmd = '-datatype int16';
    end
    if hdr.dt(1) == 512 %uint16
        cmd = '-datatype uint16';
    end  
end;
%run command
cmd = sprintf('%s %s -force -quiet "%s" "%s"', exenam, cmd, tempname, fname);
fprintf('Running (make sure this acquisition did not use partial fourier): %s\n', cmd);
status = systemSub(cmd);
if status ~= 0, error('unable to run command: %s', cmd); end;
delete(tempname);
%end degibbsSub()

function fname = dwidenoiseSub(fname)
if ~exist(fname, 'file'), warning('Unable to find %s', fname); return; end;
[pth,nam,ext] = filepartsSub(fname);
%copy bvec/bval
inname = fullfile(pth,nam);
outname = fullfile(pth,[nam, 'd']);
copyfile([inname '.bvec'], [outname '.bvec']);
copyfile([inname '.bval'], [outname '.bval']);
inname = [inname ext];
fname = [outname ext];
exenam = fullfile(nii_mrtrix_pth, 'dwidenoise');
if ~exist(exenam, 'file')
    error('unable to find %s', exenam); 
end
cmd = sprintf('%s -force -quiet "%s" "%s"', exenam, inname, fname);
fprintf('Running : %s\n', cmd);
status = systemSub(cmd);
if status ~= 0, error('unable to run command: %s', cmd); end;
%end dwidenoise()

function [pth,nam,ext,num] = filepartsSub(fname)
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
%end filepartsSub()

function status  = systemSub (cmd)
% Save library paths
MatlabPath = getenv('LD_LIBRARY_PATH');
% Make Matlab use system libraries
setenv('LD_LIBRARY_PATH',getenv('PATH'));
fprintf('%s\n',cmd);
status = system(cmd);
% Reassign old library paths
setenv('LD_LIBRARY_PATH',MatlabPath);
%systemSub