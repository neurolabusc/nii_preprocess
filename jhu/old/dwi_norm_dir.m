function dwi_norm_dir (voiDir)
%Process all the voi files in a folder

if isempty(vioDir), voiDir = pwd; end;
isTrace = false;
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
    dwi_norm(dwiname, vname, isTrace)
end
%dwi_norm_dir