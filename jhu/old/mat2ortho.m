function mat2ortho(path)
%Display all math files in folder path
% psnm : (optional) append data to post-script file
%Examples

if exist('path','var'), cd(path); end;
[~,n] = fileparts(pwd);

v = dir( '*.mat');
v ={v.name}';
for i = 1: size(v,1)
    nii_mat2ortho(char(deblank(v(i,:))), [n, '.ps']);
end