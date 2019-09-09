function varargout = coreg_DKI(fn_source,fn_target,fn_other,dest_d)
%fn_target: target image (i.e. what you're registering to)
%fn_source: source image (i.e. image to register / get matrix from)
%other: matrix of strings listing all images to apply matrix to (including source if you want it outputted)
%dest_d: directory to store coregistered images in

%Format this function so it creates fn_other like the output from
%spm_select would be

%function coregister(fn_target, fn_source, folder_other)
%fn_other = spm_select('fplist', folder_other, fn_other_filt);
% coregistration and reslicing parameters

estflg.cost_fun = 'nmi';
estflg.sep      = [4 2];
estflg.tol      = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
estflg.fwhm     = [7 7];
wrtflg.interp   = 7;
wrtflg.wrap     = [0 0 0];
wrtflg.mask     = 0;
wrtflg.mean     = 0;
wrtflg.which    = 1;

hdr_trg = spm_vol(fn_target);
hdr_src = spm_vol(fn_source);

x  = spm_coreg(hdr_trg, hdr_src, estflg);

M  = inv(spm_matrix(x));

MM = zeros(4, 4, size(fn_other, 1));

for j=1:size(fn_other, 1)
    MM(:,:,j) = spm_get_space(strtrim(fn_other(j,:)));
end

for j = 1:size(fn_other, 1)
    spm_get_space(strtrim(fn_other(j,:)), M*MM(:,:,j));
end;

fn_other = str2mat(fn_target, fn_other);

%Check to see if user wants to ouput transformation matrix (do this before
%reslicing incase user over-writes original source)
if nargout>0
    ht = spm_vol(fn_other(1,:)); 
    hs = spm_vol(fn_other(2,:)); 
    varargout{1} = ht.mat\hs.mat; 
end

%Reslice
spm_reslice(fn_other, wrtflg);

%Rename images
for i = 2:size(fn_other,1)
    fn = strtrim(fn_other(i,:));
    idx = find(fn==filesep);
    source = fullfile(fn(1:idx(end)-1),['r' fn(idx(end)+1:end)]);
    name = fn(idx(end)+1:end-4);
    
    if~isdir(dest_d); mkdir(dest_d); end
    movefile(source,fullfile(dest_d,[name '.nii']),'f');
end

fprintf('\nRegistration Complete...\n\n')