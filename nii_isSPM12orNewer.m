function nii_isSPM12orNewer
%check that SPM is installed and is at least release 6225
if exist('spm','file') ~= 2, error('Please install SPM12 or later'); end;
[v,r] = spm('Ver','',1);
if ~contains(v,'SPM12'), error('%s: Only tested on SPM12', mfilename); end
r = str2double(r);
if (r < 7771), error('%s: please update SPM12', mfilename); end
if (r > 7592) && (r < 7965)
    fnm = which('spm_preproc8');
    txt = fileread(fnm);
    badStr = 'param(6)*scal^2;';
    if contains(txt,badStr)
       warning("see https://github.com/james-cole/brainageR/issues/3\n")
       error('change "%s" to "param(6)*scal;" in "%s"', badStr, fnm);
    end %bad string
end %r 7593..7964
%end isSPM12orNewer()