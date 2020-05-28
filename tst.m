function tst()


%imgs.ASL = '/home/chris/neuro/master/In/M2024/ASL_M2024_R01.nii';
%basilDir = '/home/chris/neuro/master/In/M2024/BASIL';
%matName = '/home/chris/neuro/master/In/M2024/m.mat';
%basil2mat(basilDir, matName)

fnm = '/media/chris/Chris5TB/Universe/master/M2029/Lesion_Mapping/T2_LM1072.nii.gz';
checkLesion(fnm);

function checkLesion(fnm);
[p,n,x] = fileparts(fnm);
les = dir(fullfile(p, 'LESION*'));
if isempty(les)
    error('Unable to find LESION file in: %s\n', p);
end
les = fullfile(p, les(1).name);
lhdr = nii_tool('hdr', les);
hdr = nii_tool('hdr', fnm);
if ((lhdr.dim(2) ~= hdr.dim(2)) || (lhdr.dim(3) ~= hdr.dim(3)) || (lhdr.dim(4) ~= hdr.dim(4)) )
    error('Dims do not match: %s %s\n', fnm, les);
end
dx = [lhdr.srow_x - hdr.srow_x; lhdr.srow_y - hdr.srow_y; lhdr.srow_z - hdr.srow_z];
dx = max(abs(dx(:)));
if (dx > 0.001)
    
    %fprintf( ['i srow_x: ', repmat('%g ', 1, numel(hdr.srow_x)), '\n'], hdr.srow_x);
    %fprintf( ['l srow_x: ', repmat('%g ', 1, numel(lhdr.srow_x)), '\n'], lhdr.srow_x);
    %fprintf( ['i srow_y: ', repmat('%g ', 1, numel(hdr.srow_y)), '\n'], hdr.srow_y);
    %fprintf( ['l srow_y: ', repmat('%g ', 1, numel(lhdr.srow_y)), '\n'], lhdr.srow_y);
    %fprintf( ['i srow_z: ', repmat('%g ', 1, numel(hdr.srow_z)), '\n'], hdr.srow_y);
    %fprintf( ['l srow_z: ', repmat('%g ', 1, numel(lhdr.srow_z)), '\n'], lhdr.srow_y);
    fprintf('SForm(%g):\t%s\n', les);
    [p,n,x] = fileparts(les);
    orig = fullfile(p, ['_',n,x]);
    copyfile(les, orig);
    nii = nii_tool('load', les);
    nii.hdr.srow_x = hdr.srow_x;
    nii.hdr.srow_y = hdr.srow_y;
    nii.hdr.srow_z = hdr.srow_z;
    nii_tool('save', nii, les);
    return;
end
fprintf('OK: %s\n', fnm);