function naming_batch (baseDir)
%find relevant images for patients
% baseDir: root folder that holds images
%Data organized as follows
% ~/p1/p1.nii
% ~/rb/rb.nii
% ~/sv/sv.nii


if ~exist('baseDir','var') || isempty(baseDir)
    baseDir = pwd; %uigetdir('','Pick folder that contains all subjects');
end
nameFMRIs = subFolderSub(baseDir);
nameFMRIs = sort(nameFMRIs);
if numel(nameFMRIs) < 1, error('No images found'); end;
p.t1name = -1; %NONE
p.TRsec = 10; %repeat time is 1.92 seconds per volume
p.slice_order = -1; %SKIP
p.phase =  ''; %phase image from fieldmap
p.magn = ''; %magnitude image from fieldmap
%statistical information (optional: if not provided not statitics)
p.names{1} = 'naming';
p.names{2} = 'abstract';
%onsets for 1st session, onsets{session, cond}
p.onsets{1,1} = [16.734	26.031	40.609	51.641	70.141	88.625	109.469	128.516	138.75	169.813	186.813	201.141	208.859	219.875	250.672	259.875	266.359	278.953	307.234	328.422	337.219	349.531	358.922	367.844	377.547	388.234	401.109	410.891	421	440.406	451.469	468.563	477.625	506.578	520.594	541.813	549.75	581.125	587.406	597.563];
p.onsets{1,2} = [11.547	57.813	80.797	96.266	119.016	147.891	158	180.672	228.016	237.36	287.157	301.813	320.891	426.75	456.75	491.704	496.672	530.641	556.344	569.829];


p.onsets{1,1} = p.onsets{1,1} - (p.TRsec/2); %we will time from start not middle of acq
p.onsets{1,2} = p.onsets{1,2} - (p.TRsec/2); %we will time from start not middle of acq
%onsets for 2nd session, onsets{session, cond}
p.duration{1} = 3; %duration 1 for events, longer for blocks
p.mocoRegress = true; %should motion parameters be included in statistics?
startPwd = pwd;
% for i = 1 : numel(nameFMRIs)
%     fname = nameFMRIs{i};
%     [pth,nam,ext] = fileparts(fname);
%     matname = fullfile(pth, [nam, '.mat']);
%     if exist(matname, 'file') delete(matname); end;
%     cropName4D = nii_cropVolumes (fname, 0, 60);
%     if ~strcmp(cropName4D,fname)
%         
%         movefile(fname, fullfile(pth,[nam '120' ext]));
%         movefile(cropName4D, fname); 
%         fprintf('Created cropped image\n');
%     end
% end
for i = 1 : numel(nameFMRIs)
    fname = nameFMRIs{i};
    %fname = '/Users/rorden/Downloads/naming80_normals/p3/p3.nii';
    [pth,nam,ext] = fileparts(fname);
    cd(pth);
    statdir = fullfile(pth,nam);
    if exist(statdir, 'file'), rmdir(statdir, 's'); end;
    p.fmriname = fullfile('',[nam,ext]); 
    nii_batch12 (p);
end
cd(startPwd);


fnms = '';
tfnms = '';
cfnms = '';
for i = 1 : numel(nameFMRIs)
    fname = nameFMRIs{i};
    [pth,nam,ext] = fileparts(fname);
    fnms = strvcat(fnms, fullfile(pth,['wbmean', nam,ext])); %#ok<REMFF1>
    tfnms = strvcat(tfnms, fullfile(pth,nam,['spmT_0002.nii'])); %#ok<REMFF1>
    cfnms = strvcat(cfnms, fullfile(pth,nam,['con_0002.nii'])); %#ok<REMFF1>
    
end
nii_mean_stdev (fnms, true, 'intensity_');
nii_mean_stdev (tfnms, false, 't_');
nii_mean_stdev (cfnms, false, 'beta_');
%end naming_batch()

function nameFMRI=subFolderSub(pathFolder)
d = dir(pathFolder);
isub = [d(:).isdir];
nameFMRI = {d(isub).name}';
nameFMRI(ismember(nameFMRI,{'.','..'})) = [];
for i = numel(nameFMRI): -1 : 1
   nameFMRI{i} = fullfile(pathFolder, deblank(nameFMRI{i}), [deblank(nameFMRI{i}) '.nii']);
   
   if ~exist(nameFMRI{i},'file'), nameFMRI{i} = []; end;
   %fprintf('%s\n',nameFMRI{i});
   
end
%end subFolderSub()
