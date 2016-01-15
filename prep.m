function prep

f = dir('*.nii');
for i = 1 : numel(f)
    fname = f(i).name;
    nii_setOrigin12(fname, 2, true);
    
end;