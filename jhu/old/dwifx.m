function dwifx (isTrace)
%Process all the voi files in a folder
% isTrace: if true, single volume images are assumed to be TRACE, else single volumes assumed B0

voiDir = pwd;
if ~exist('isTrace','var'), isTrace = false; end;
if isempty(which('spm')) || ~strcmp(spm('Ver'),'SPM12'), error('SPM12 required'); end;
v = dir( fullfile(voiDir, '*.voi'));
vois = {v.name}';
missing = 0;
fail = 0;
for i = 1: size(vois,1) %: -1 : 5
    vname = char(deblank(vois(i,:)));
    if isempty(vname) || (vname(1) == '.'), continue; end;
    [~,n,x] = fileparts(vname);
    vname = fullfile(voiDir, [n,x]); %append path
    d = dir( fullfile(voiDir, [n,'*.nii']));
    dwi = {d.name}';
    if numel(dwi) < 1, missing = missing + 1; fprintf('%d\tmissing\t%s\n', missing, n); continue; end;
    dwiname = fullfile(voiDir, char(deblank(dwi(1,:))) );%append path
    %fprintf('%d normalize\t%s\t%s\n', i, dwiname, n);
    ok = normDwiSub(dwiname, vname, isTrace);
    if ~ok
        fail = fail + 1;
        fprintf('Error\t%d\t%s\n', fail, dwiname);
    end
end
%dwifx

function ok = normDwiSub(dwi, lesion, isTrace)
ok = true; %if mat exists assume success
if ~exist('isTrace','var'), isTrace = false; end;
[p,n] = fileparts(lesion);
matName = fullfile(p, [n, '.mat']);
if exist(matName,'file'), return; end;
ok = false; %assume failure
hdr = spm_vol(dwi);
trace = [];
if numel(hdr) > 1 %if two volumes, determine which is B=0 and which is trace
   img = spm_read_vols(hdr);
   img1 = img(:,:,:,1);
   img2 = img(:,:,:,2);
   if mean(img1(:)) < mean(img2(:))
    img2 = img(:,:,:,1);
    img1 = img(:,:,:,2);
   end
   hdr = hdr(1);
   [pth nm ext] = spm_fileparts(hdr.fname);
   dwi = fullfile(pth, ['b0_' nm ext]);
   trace = fullfile(pth, ['tr_' nm ext]);
   hdr.fname = dwi;
   spm_write_vol(hdr,img1);
   hdr.fname = trace;
   spm_write_vol(hdr,img2);
end
lesion = unGzSub(lesion); %.voi -> .nii
if isempty(trace)
    if exist('isTrace','var') && isTrace
        p = fileparts(which(mfilename));
        template = fullfile(p,'mean_TRACE_of_381_8mmFWHM.nii');
    else
        template  = fullfile(spm('Dir'),'toolbox','OldNorm','EPI.nii');
    end
    if ~exist(template,'file'), error('Can not find %s',template); end;
    [pth nm ext] = spm_fileparts(hdr(1).fname);
    %return; %to report broken files!
    if isTrace
        fprintf('%s expect 2 images (B0 and Trace). ASSUMING TRACE!',dwi);
        trace = fullfile(pth, ['tr_' nm ext]);
        copyfile(dwi,trace);
        dwi = trace;
        nii_setOrigin12({dwi,lesion}, 3, false); %T2 with yoked Lesion
        oldNormSub( {dwi,lesion}, template, 0);
        bo = [];        
    else
        fprintf('%s expect 2 images (B0 and Trace). ASSUMING B0!',dwi);
        b0 = fullfile(pth, ['b0_' nm ext]);
        copyfile(dwi,b0);
        dwi = b0;
        nii_setOrigin12({dwi,lesion}, 3, false); %T2 with yoked Lesion
        oldNormSub( {dwi,lesion}, template, 0);
        trace = [];
    end;
else
    template  = fullfile(spm('Dir'),'toolbox','OldNorm','EPI.nii');
    nii_setOrigin12({dwi,lesion,trace}, 3, false); %T2 with yoked Lesion
    oldNormSub( {dwi,lesion,trace}, template, 0);
end;
lesion = prefixSub('w', lesion);
nii_nii2mat(lesion, 'lesion', matName); %1
ok = true;
%save DWI data
if isempty(trace), return; end;
trace = prefixSub('w', trace);
if ~exist(trace,'file'), return; end;
hdr = spm_vol(trace);
img = spm_read_vols(hdr);
stat = [];
stat.DWI.hdr = hdr;
stat.DWI.dat = img;
old = load(matName);
stat = nii_mergestruct(stat,old); %#ok<NASGU>
save(matName,'-struct', 'stat');

%normDwiSub()

function nam = prefixSub(pre, nam)
if isempty(nam), return; end;
[p, n, x] = spm_fileparts(nam);
nam = fullfile(p, [pre, n, x]);
%end prefixSub()

function oldNormSub(src, tar, smoref, reg, interp)
%coregister T2 to match T1 image, apply to lesion
if isempty(src) || isempty(tar), return; end;
if ~exist('smoref','var'), smoref = 0; end;
if ~exist('reg','var'), reg = 1; end;
if ~exist('interp','var'), interp = 1; end;
matlabbatch{1}.spm.tools.oldnorm.estwrite.subj.source = src(1);
matlabbatch{1}.spm.tools.oldnorm.estwrite.subj.wtsrc = '';
matlabbatch{1}.spm.tools.oldnorm.estwrite.subj.resample = src(:);
matlabbatch{1}.spm.tools.oldnorm.estwrite.eoptions.template = {tar};
matlabbatch{1}.spm.tools.oldnorm.estwrite.eoptions.weight = '';
matlabbatch{1}.spm.tools.oldnorm.estwrite.eoptions.smosrc = 8;
matlabbatch{1}.spm.tools.oldnorm.estwrite.eoptions.smoref = smoref; % <-- !!!
matlabbatch{1}.spm.tools.oldnorm.estwrite.eoptions.regtype = 'mni';
matlabbatch{1}.spm.tools.oldnorm.estwrite.eoptions.cutoff = 25;
matlabbatch{1}.spm.tools.oldnorm.estwrite.eoptions.nits = 16;
matlabbatch{1}.spm.tools.oldnorm.estwrite.eoptions.reg = reg;
matlabbatch{1}.spm.tools.oldnorm.estwrite.roptions.preserve = 0;
%matlabbatch{1}.spm.tools.oldnorm.estwrite.roptions.bb = [nan nan nan; nan nan nan];%[-78 -112 -70; 78 76 85];
matlabbatch{1}.spm.tools.oldnorm.estwrite.roptions.bb = [-78 -112 -70; 78 76 85];
matlabbatch{1}.spm.tools.oldnorm.estwrite.roptions.vox = [1 1 1];
matlabbatch{1}.spm.tools.oldnorm.estwrite.roptions.interp = interp;
matlabbatch{1}.spm.tools.oldnorm.estwrite.roptions.wrap = [0 0 0];
matlabbatch{1}.spm.tools.oldnorm.estwrite.roptions.prefix = 'w';
spm_jobman('run',matlabbatch);
%end oldNormSub()

function fnm = unGzSub (fnm)
[pth,nam,ext] = spm_fileparts(fnm);
if strcmpi(ext,'.gz') %.nii.gz
    fnm = char(gunzip(fnm));
elseif strcmpi(ext,'.voi') %.voi ->
    onam = char(gunzip(fnm));
    fnm = fullfile(pth, ['LESION_' nam '.nii']);
    movefile(onam,fnm);
end;
%end unGzSub()