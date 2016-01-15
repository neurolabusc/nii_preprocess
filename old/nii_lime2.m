function nii_lime2 (imgs)
%analyze images from stroke patients
% imgs is a structure with elements
%  imgs.name = names of subject (required)
%  imgs.t1 = name of T1 file
%  imgs.t2 = name of T2 file
%example
% imgs.name = 'P097'
% imgs.DTI=  '/Users/rorden/Desktop/pipeline/P097/DTI_P097_20140527/i.nii'
% nii_lime2(imgs)

prefolder = '/Users/rorden/Desktop/pre/';
if exist(prefolder, 'file') ~= 7
    error('Unable to find %s', pre);
 
end
if ~isfield(imgs,'name') || isempty(imgs.name)
	error('Field imgs.name required');
end
newdir = fullfile(prefolder, imgs.name);
if ~exist(newdir, 'file'), mkdir(newdir); end;
%next copy and rename files
%diary( fullfile(prefolder, 'notes.txt') );

if isfield(imgs,'DTI') && ~isempty(imgs.DTI)
    copyDtiSub(imgs.DTI, newdir, ['DTI_' imgs.name ]);
end
if isfield(imgs,'ASL') && ~isempty(imgs.ASL)
    copyImgSub(imgs.ASL, newdir, ['PASL_' imgs.name ]);
end

if isfield(imgs,'T1') && ~isempty(imgs.T1)
    copyImgSub(imgs.T1, newdir, ['T1_' imgs.name ]);
end
if isfield(imgs,'T2') && ~isempty(imgs.T2)
    copyImgSub(imgs.T2, newdir, ['T2_' imgs.name ]);
end

if isfield(imgs,'fMRI') && ~isempty(imgs.fMRI)
    copyImgSub(imgs.fMRI, newdir, ['fMRI_' imgs.name ]);
end
if isfield(imgs,'rest') && ~isempty(imgs.rest)
    copyImgSub(imgs.rest, newdir, ['REST_' imgs.name ]);
end
if isfield(imgs,'Lesion') && ~isempty(imgs.Lesion)
    copyImgSub(imgs.Lesion, newdir, ['LS_' imgs.name ]);
end
%diary OFF
%end nii_lime2


function copyDtiSub(innames, outdir, outnames)
if isempty(innames), return; end;
if ~exist('nameIncrement','var'), nameIncrement = 0; end;
for i = 1: size(innames,1)
	inname = deblank(innames(i,:)); 
    [pth, nam, ext] = spm_fileparts(inname); %#ok<ASGLU>
    nameFiles = subFileSub(pth);
    if isempty(inname), return; end;
    numDTI = 0;
    for n = 1: size(nameFiles,1) %for each participant
        nx = deblank(nameFiles{n}); 
        [p, n, x] = spm_fileparts(nx);
        if strcmpi(x,'.bval')
            fileID = fopen(fullfile(pth, nx),'r');
            formatSpec = '%f';
            A = fscanf(fileID,formatSpec);
            fclose(fileID);
            if numel(A) > 5
                numDTI = numDTI + 1;
                if strfind(lower(n),lower('PAs0')) %posterior->anterior
                    tag = 'PA';
                elseif strfind(lower(n),lower('APs0')) %anterior->posterior
                    tag = 'AP';
                elseif numDTI > 1
                    tag = num2str(numDTI-1);
                else
                    tag = '';
                end
                outname = outnames;
                if i > 1, outname = [outname '_' num2str(i-1)]; end; %#ok<AGROW>
                copyFileSub(fullfile(pth, [n '.bval']), outdir, [tag outname]); %bval
                copyFileSub(fullfile(pth, [n '.bvec']), outdir, [tag outname]);
                copyImgSub(fullfile(pth, [n '.nii']), outdir, [tag outname]);
            end
        end
    end
end %for each DTI
%copyDtiSub()

function nameFiles=subFileSub(pathFolder)
d = dir(pathFolder);
isub = ~[d(:).isdir];
nameFiles = {d(isub).name}';
%end subFileSub()

function copyImgSub(innames, outdir, newnames)
if isempty(innames), return; end;
if ~exist('nameIncrement','var'), nameIncrement = 0; end;
for i = 1: size(innames,1)
    inname = deblank(innames(i,:));
    newname = newnames;
    if i > 1, newname = [newname '_' num2str(i-1)]; end; %#ok<AGROW>
    [pth, nam, ext] = spm_fileparts(inname); 
    fprintf('Copy\t%s\t%s\n',inname, newname);
    if strcmpi(ext,'.gz') %remove both .nii and .gz from .nii.gz
        [pth2, nam, ext] = spm_fileparts(nam);   %#ok<NASGU,ASGLU>
    end
    filenamegz = fullfile(pth, [nam '.nii.gz']);
    filename = fullfile(pth, [nam '.nii']);
    filenamevoi = fullfile(pth, [nam '.voi']);
    %fprintf('Copying %s\n', newname);
    if exist(filenamegz,'file') 
        copyFileSub(filenamegz, outdir, newname);
    elseif exist(filenamevoi,'file') 
        outname = copyFileSub(filenamevoi, outdir, newname);
        [pth, nam] = spm_fileparts(outname); 
        movefile(outname, fullfile(pth, [nam '.nii.gz']) );
    else %.nii image
        if ~exist(filename,'file')
           error('Unable to find image named %s or %s\n',filename, filenamegz); 
        end
        outname = copyFileSub(filename, outdir, newname);
        gzip(outname);
        delete(outname);
    end;
end
%end copyImgSub()

function outname = copyFileSub(inname, outdir, newname)
if isempty(inname), return; end;
[pth, nam, ext] = spm_fileparts(inname); %#ok<ASGLU>
if strcmpi(ext,'.gz')
    ext = '.nii.gz';
end
outname = fullfile(outdir, [newname ext]);
copyfile(inname, outname);
%end copyFileSub()

% function nam = namSub(fullfilename)
% [~,nam] = fileparts(fullfilename);
% %end nameSub
% 
% function img = checkFilename (nam)
% img = nam;
% if ~isempty(img) && exist(img, 'file') ~= 2 
%     fprintf('%s warning: unable to find image named %s\n',mfilename,img);
%     img = '';
% end  
% %end checkFilename()