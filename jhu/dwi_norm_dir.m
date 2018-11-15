function dwi_norm_dir (voiDir)
%Process all the voi files in a folder

if ~exist('voiDir', 'var'), 
    voiDir = pwd; 
    %voiDir = uigetdir;
end;
isTrace = false;
if isempty(which('spm')) || ~strcmp(spm('Ver'),'SPM12'), error('SPM12 required'); end;
v = dir( fullfile(voiDir, '*.voi'));
vois = {v.name}';
missing = 0;
for i = 1 : numel(vois)
    vname = char(deblank(vois(i,:)));
    if isempty(vname) || (vname(1) == '.'), continue; end;
    [~,n,x] = fileparts(vname);
    vname = fullfile(voiDir, [n,x]); %append path
    d = dir( fullfile(voiDir, [n,'*.nii']));
    dwi = {d.name}';
    if numel(dwi) < 1 %another JHU naming convention: AGS3790_DWI.voi is associated with AGS3790.nii.gz!!
        n = strrep(n, '_DWI', '');
        d = dir( fullfile(voiDir, [n,'*.nii.gz']));
        dwi = {d.name}';
        if numel(dwi) < 1
            d = dir( fullfile(voiDir, [n,'*.nii']));
            dwi = {d.name}';    
        end
    end
    if numel(dwi) < 1, missing = missing + 1; fprintf('%d\tmissing DWI named:\t%s\n', missing, n); continue; end;
    dwiname = fullfile(voiDir, char(deblank(dwi(1,:))) );%append path
    dwi_norm(dwiname, vname, isTrace)
end
%dwi_norm_dir