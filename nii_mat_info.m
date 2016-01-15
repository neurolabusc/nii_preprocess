function nii_mat_info(baseDir)

if ~exist('baseDir','var')
    baseDir = pwd;
end
nameFiles = subFileSub(baseDir);
nameFiles = sort(nameFiles); %take first file for multiimage sequences, e.g. ASL
fprintf('Found %d mat files (subjects) in %s\n',size(nameFiles,1), baseDir);
fprintf('img\tlesionVolume\tASLvols\tLgm\tRgm\tLwm\tRwm\n');
for i=1:size(nameFiles,1)
    fnm = char(nameFiles(i));
    [p,n,x] = fileparts(fnm); %#ok<ASGLU>
    if strcmpi(x,'.mat')
        m = load(fnm);
        lesionVol = 0;
        if isfield(m,'lesion')
            m.lesion.dat(~isfinite(m.lesion.dat)) = 0;
            lesionVol = sum(m.lesion.dat(:)) ;
        end
        if isfield(m,'cbf')
            fprintf('%s\t%g\t%g\t%g\t%g\t%g\t%g\n', n, lesionVol, m.cbf.nV, m.cbf.c1L,  m.cbf.c1R, m.cbf.c2L,  m.cbf.c2R);
        end
     end
    
    
end
%end nii_mat_info

function nameFiles=subFileSub(pathFolder)
d = dir(pathFolder);
isub = ~[d(:).isdir];
nameFiles = {d(isub).name}';
%end subFileSub()

