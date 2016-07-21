function nii_merge_travis

label = jhuLabelSub;

baseDir = '/Users/rorden/Desktop/matfiles/';
mergeDir = '/Users/rorden/desktop/LM_Mat_Files_20140902';
d = dir(baseDir);
isub = ~[d(:).isdir];
nameFiles = {d(isub).name}';
for s = 1: size(nameFiles,1)
    subj = [baseDir deblank(nameFiles{s})];
    [~, nam, ext] = fileparts(subj);
    if ~strcmpi(ext,'.mat'), continue; end; %only mat files
    tsubj = fullfile(mergeDir,[nam, 'density.mat']);
    if exist(tsubj,'file')
        %dsubj
        mergeSub(subj, tsubj, label, 'dti_jhu');
    end
    tsubj = fullfile(mergeDir,[nam, 'fiber_count.mat']);
    if exist(tsubj,'file')
        mergeSub(subj, tsubj, label, 'dtifc_jhu');
    end
end
%end nii_merge_travis()

function label = jhuLabelSub
pth = fileparts(which('NiiStat'));
if isempty(pth), error('Unable to find NiiStat'); end;
pth = [pth filesep 'roi' filesep 'jhu.txt'];
if ~exist(pth,'file'), error('Unable to find %s',pth); end;
fid = fopen(pth);  % Open file
label=[];
tline = fgetl(fid);
while ischar(tline)
    %disp(tline)
    label=strvcat(label,tline); %#ok<REMFF1>
    tline = fgetl(fid);
end
fclose(fid);
%end labelSub()

function mergeSub(matName, txtName, labels, statname)
%size(label,1)
m = spm_load(txtName);
if size(m,1) ~= size(m,2), error('Matrix not square (number of columns and rows differ)'); end;
if size(labels,1) ~= size(m,1), error('Wrong number of items in matrix'); end;
stat.(statname).label = labels;
stat.(statname).r = m;
if length(matName) < 1, return; end
if exist(matName,'file')
    old = load(matName);
    %old = rmfield(old,statname);
    if isfield(old,'rest_aal')
        if max(old.rest_aal.r(:)) == min(old.rest_aal.r(:))
            fprintf('WARNING: Please check resting state data of %s\n',matName);
            old = rmfield(old,'rest_aal');
            old = rmfield(old,'rest_aalcat');
            old = rmfield(old,'rest_bro');
            old = rmfield(old,'rest_cat');
            old = rmfield(old,'rest_fox');
            old = rmfield(old,'rest_jhu');
        end
      end
    stat = nii_mergestruct(stat,old);
end
save(matName, '-struct', 'stat');
%end mergeSub()



