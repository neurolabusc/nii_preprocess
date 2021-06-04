function basil2mat(dir, matName)
%Copy BASIL results to NiiStat format Mat fie
% dir : base BASIL folder
% matName : NiiStat MatFile
%Examples
% basil2mat('/M2030/BASIL', '/M2030_T1_lime.mat')


%matName = '/home/coffee/Desktop/POLAR/T1_mprage_ns_sag_p2_iso_1mm_192_M2183_POLAR1033_session1_lime.mat';
%dir = '/home/coffee/Desktop/POLAR/BASIL'

sdir = fullfile(dir,'struc.anat');
c1 = fullfile(sdir, 'T1_fast_pve_1.nii.gz');
c2 = fullfile(sdir, 'T1_fast_pve_2.nii.gz');
sdir = fullfile(dir,'native_space');
native = fullfile(sdir, 'perfusion_calib.nii.gz'); 
sdir = fullfile(dir,'std_space');
std = fullfile(sdir, 'perfusion_calib.nii.gz');
if ~exist(c1,'file'), error('Unable to find %s', c1); end
if ~exist(c2,'file'), error('Unable to find %s', c2); end
if ~exist(std,'file'), error('Unable to find %s', std); end
if ~exist(native,'file'), error('Unable to find %s', native); end
if ~exist(matName, 'file'), error('Unable to find %s', matName); end;
addNiiSub(c1, 'basilC1', matName); %1
addNiiSub(c2, 'basilC2', matName); %1
addNiiSub(std, 'basilStd', matName); %1
addNiiSub(native, 'basilNative', matName); %1
%end basil2mat()

function addNiiSub(niiname, Voxfield, matName)
niiname = unGzSub (niiname);
[p,n,~] =spm_fileparts(niiname);
hdr = spm_vol(niiname);
img = spm_read_vols(hdr);
mn = min(img(:));
mx = max(img(:));
if mn == mx 
    error('No variability in image %s',niiname);
end
stat = [];
if ndims(img) == 3
    stat.(Voxfield).hdr = hdr;
    stat.(Voxfield).dat = img;
end
if exist(matName,'file')
    old = load(matName);
    stat = nii_mergestruct(stat,old); %#ok<NASGU>
end
save(matName,'-struct', 'stat');
%end addNiiSub()

function fnm = unGzSub (fnm)
[pth,nam,ext] = spm_fileparts(fnm);
if strcmpi(ext,'.gz') %.nii.gz
    onam = fnm;
    fnm = char(gunzip(fnm));
    delete(onam); % gunzip() does not delete the original file DPR 20200318   
elseif strcmpi(ext,'.voi') %.voi -> 
    onam = char(gunzip(fnm));
    fnm = fullfile(pth, [nam '.nii']);
    movefile(onam,fnm);
end  
%end unGzSub()


