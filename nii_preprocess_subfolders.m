function nii_preprocess_subfolders(pth)
%find all limegui.mat files in subfolders and re-process them
% pth: parent folder
%Examples
% nii_preprocess_subfolders('~/a'); %would process ~/a/b/1_limegui.mat.mat, ~/a/c/2_limegui.mat
% nii_preprocess_subfolders; % search from current working directory
%Notes: will skip folders with "_" in name, e.g. will skip M2002:
% pth/M2000 pth/M2001 pth/M2002_dementia pth/M2003
if ~exist('pth','var'), pth = pwd; end;
f = subFolderSub(pth);
if isempty(f), error('No folders in parent folder %s', pth); end;
%global ForcefMRI;  ForcefMRI = true; warning('FORCED fMRI REPROCESSING'); %comment line for auto-processing
%global ForceRest;  ForceRest =  warning('FORCED REST REPROCESSING'); %comment line for auto-processing
%global ForceASL;  ForceASL = true; warning('FORCED ASL REPROCESSING'); %comment line for auto-processing
global ForceDTI; ForceDTI = [];
global ForcefMRI;  ForcefMRI = []; %comment line for auto-processing
global ForceRest;  ForceRest = [];   %comment line for auto-processing
global ForceASL;  ForceASL = []; %comment line for auto-processing
isPreprocess = false; %process data
isUpdateVers = true; %set LIME version to latest
isCopymat = true; %copy output files to M.vers and LIME.vers folders
isMakeModalityTable = true;

MasterID = {'M2025', 'M2071', 'M2005', 'M2036', 'M2006', 'M2002', 'M2069', 'M2007', 'M2020', 'M2141', 'M2142', 'NA', 'NA', 'M2046', 'M2014', 'M2074', 'M2061', 'M2078', 'M2143', 'M2072', 'M2144', 'M2145', 'M2146', 'M2075', 'M2059', 'M2040', 'M2076', 'M2016', 'M2147', 'M2148', 'M2082', 'M2079', 'M2070', 'M2149', 'M2150', 'M2151', 'M2152', 'M2085', 'M2153', 'M2086', 'M4138', 'M2154', 'M4148', 'M2088', 'M2044', 'M2087', 'M2155', 'M2084', 'M2030', 'M2156', 'M2157', 'M2158', 'M2159', 'M2160', 'M2094', 'M2161', 'M2162', 'M2031', 'M2017', 'M2163', 'M2103', 'M4180', 'M2164', 'M2106', 'M2165', 'M2110', 'M2109', 'M2114', 'M2166', 'M2112', 'M2111', 'M2029', 'M4209', 'M2119', 'M2125', 'M2126', 'M2124', 'M2121', 'M2117', 'M2130', 'M2127', 'M2133', 'M2136', 'M2129', 'M2123', 'M2135', 'M2134', 'M2139', 'M2138', 'M2170'};
LimeID = {'LM1001', 'LM1002', 'LM1003', 'LM1004', 'LM1005', 'LM1006', 'LM1007', 'LM1008', 'LM1009', 'LM1010', 'LM1011', 'LM1012', 'LM1013', 'LM1014', 'LM1015', 'LM1016', 'LM1017', 'LM1018', 'LM1019', 'LM1020', 'LM1021', 'LM1022', 'LM1023', 'LM1024', 'LM1025', 'LM1026', 'LM1027', 'LM1028', 'LM1029', 'LM1030', 'LM1031', 'LM1032', 'LM1033', 'LM1034', 'LM1035', 'LM1036', 'LM1037', 'LM1038', 'LM1039', 'LM1040', 'LM1041', 'LM1042', 'LM1043', 'LM1044', 'LM1045', 'LM1046', 'LM1047', 'LM1048', 'LM1049', 'LM1050', 'LM1051', 'LM1052', 'LM1053', 'LM1054', 'LM1055', 'LM1056', 'LM1057', 'LM1058', 'LM1059', 'LM1060', 'LM1061', 'LM1062', 'LM1063', 'LM1064', 'LM1065', 'LM1066', 'LM1067', 'LM1068', 'LM1069', 'LM1070 ', 'LM1071', 'LM1072', 'LM1073', 'LM1074', 'LM1075', 'LM1076', 'LM1077', 'LM1078', 'LM1079', 'LM1080', 'LM1081', 'LM1082', 'LM1083', 'LM1084', 'LM1085', 'LM1086', 'LM1087', 'LM1088', 'LM1089', 'LM1090'};
%Show Master -> Lime table
%for i = 1:numel(MasterID),
%    fprintf('%s\t%s\n', MasterID{i}, LimeID{i});
%end

vers = nii_matver;
ButtonName = questdlg(sprintf('Version %.4f?', vers.lime), 'Yes', 'No');
if ~strcmpi(ButtonName,'Yes'), return; end;


%f = {'M2025';}; %for a single folder
%f = {'M2001';'M4217';'M4218'};%
if isPreprocess
    t = tic;
    n = 0;
    for i = 1: numel(f) %change 1 to larger number to restart after failure
    %for i = numel(f): -1 : 1 %change 1 to larger number to restart after failure
       cpth = char(deblank(f(i))); %local child path
       fprintf('===\t%s participant %d/%d : %s\t===\n', mfilename, i, numel(f), cpth);
       %if ~isempty(strfind(cpth,'M2039')), error('all done'); end; %to stop at specific point
       if ~isempty(strfind(cpth,'_'))
          fprintf('Warning: "_" in folder name: skipping %s\n', char(cpth) );
          continue
       end
       cpth = char(fullfile(pth,cpth)); %full child path
       nii_preprocess_gui(cpth);
       n = n + 1;
    end
    fprintf('Processed %d *limegui.mat file in %gs\n', n, toc(t))
end %if isPreprocess

if isUpdateVers
    for i = 1: numel(f) %change 1 to larger number to restart after failure
       cpth = char(deblank(f(i))); %local child path
       if ~isempty(strfind(cpth,'_'))
          fprintf('Warning: "_" in folder name: skipping %s\n', char(cpth) );
          continue
       end
       inpth = fullfile(pth, cpth);
       mfile = dir(char(fullfile(inpth,'*lime.mat')));
       if isempty(mfile), continue; end;
       if numel(mfile) > 1, warning('Multiple lime.mat files in %s', inpth); end;
       inname = fullfile(inpth, mfile(1).name);
       m = load(inname);
       if ~isfield(m,'T1') || ~isfield(m.T1,'lime'), warning('No T1 field for %s\n',inname); continue; end;
       if (m.T1.lime == vers.lime), fprintf('%s\tOK\t->\t%.4f\n', cpth, vers.lime); continue; end; %already up to date
       fprintf('%s\t%.4f\t->\t%.4f\n', cpth, m.T1.lime, vers.lime);
       m.T1.lime = vers.lime;
       save(inname, '-struct', 'm');
    end
end


if isCopymat
    vers = nii_matver;
    vers = sprintf('%.4f', vers.lime);
    %fprintf('%s\n', vers);
    %vers = '2017.0307'; %%% GY
    outpth = fullfile(pth, ['M.', vers]);
    mkdir(outpth);
    %MasterID = {'M2025', 'M2071', 'M2005', 'M2036', 'M2006', 'M2002', 'M2069', 'M2007', 'M2020', 'M2141', 'M2142', 'M2046', 'M2014', 'M2074', 'M2061', 'M2078', 'M2143', 'M2072', 'M2144', 'M2145', 'M2146', 'M2075', 'M2059', 'M2040', 'M2076', 'M2016', 'M2147', 'M2148', 'M2082', 'M2079', 'M2070', 'M2149', 'M2150', 'M2151', 'M2152', 'M2085', 'M2153', 'M2086', 'M4138', 'M2154', 'M4148', 'M2088', 'M2044', 'M2087', 'M2155', 'M2084', 'M2030', 'M2156', 'M2157', 'M2158', 'M2159', 'M2160', 'M2094', 'M2161', 'M2162', 'M2031', 'M2017', 'M2163', 'M2103', 'M4180', 'M2164', 'M2106', 'M2165', 'M2110', 'M2109', 'M2114', 'M2166', 'M2112', 'M2111', 'M2029'};
    %LimeID = {'LM1001', 'LM1002', 'LM1003', 'LM1004', 'LM1005', 'LM1006', 'LM1007', 'LM1008', 'LM1009', 'LM1010', 'LM1011', 'LM1014', 'LM1015', 'LM1016', 'LM1017', 'LM1018', 'LM1019', 'LM1020', 'LM1021', 'LM1022', 'LM1023', 'LM1024', 'LM1025', 'LM1026', 'LM1027', 'LM1028', 'LM1029', 'LM1030', 'LM1031', 'LM1032', 'LM1033', 'LM1034', 'LM1035', 'LM1036', 'LM1037', 'LM1038', 'LM1039', 'LM1040', 'LM1041', 'LM1042', 'LM1043', 'LM1044', 'LM1045', 'LM1046', 'LM1047', 'LM1048', 'LM1049', 'LM1050', 'LM1051', 'LM1052', 'LM1053', 'LM1054', 'LM1055', 'LM1056', 'LM1057', 'LM1058', 'LM1059', 'LM1060', 'LM1061', 'LM1062', 'LM1063', 'LM1064', 'LM1065', 'LM1066', 'LM1067', 'LM1068', 'LM1069', 'LM1070', 'LM1071', 'LM1072'};
    %MasterID = {'M2025', 'M2071', 'M2005', 'M2036', 'M2006', 'M2002', 'M2069', 'M2007', 'M2020', 'M2141', 'M2142', 'M2046', 'M2014', 'M2074', 'M2061', 'M2078', 'M2143', 'M2072', 'M2144', 'M2145', 'M2146', 'M2075', 'M2059', 'M2040', 'M2076', 'M2016', 'M2147', 'M2148', 'M2082', 'M2079', 'M2070', 'M2149', 'M2150', 'M2151', 'M2152', 'M2085', 'M2153', 'M2086', 'M4138', 'M2154', 'M4148', 'M2088', 'M2044', 'M2087', 'M2155', 'M2084', 'M2030', 'M2156', 'M2157', 'M2158', 'M2159', 'M2160', 'M2094', 'M2161', 'M2162', 'M2031', 'M2017', 'M2163', 'M2103', 'M4180', 'M2164', 'M2106', 'M2165', 'M2110', 'M2109', 'M2114', 'M2166', 'M2112', 'M2111', 'M2029', 'M4209', 'M2119', 'M2125', 'M2126', 'M2124', 'M2121', 'M2117', 'M2130', 'M2127', 'M2133', 'M2136', 'M2129', 'M2135'};
    %LimeID = {'LM1001', 'LM1002', 'LM1003', 'LM1004', 'LM1005', 'LM1006', 'LM1007', 'LM1008', 'LM1009', 'LM1010', 'LM1011', 'LM1014', 'LM1015', 'LM1016', 'LM1017', 'LM1018', 'LM1019', 'LM1020', 'LM1021', 'LM1022', 'LM1023', 'LM1024', 'LM1025', 'LM1026', 'LM1027', 'LM1028', 'LM1029', 'LM1030', 'LM1031', 'LM1032', 'LM1033', 'LM1034', 'LM1035', 'LM1036', 'LM1037', 'LM1038', 'LM1039', 'LM1040', 'LM1041', 'LM1042', 'LM1043', 'LM1044', 'LM1045', 'LM1046', 'LM1047', 'LM1048', 'LM1049', 'LM1050', 'LM1051', 'LM1052', 'LM1053', 'LM1054', 'LM1055', 'LM1056', 'LM1057', 'LM1058', 'LM1059', 'LM1060', 'LM1061', 'LM1062', 'LM1063', 'LM1064', 'LM1065', 'LM1066', 'LM1067', 'LM1068', 'LM1069', 'LM1070', 'LM1071', 'LM1072', 'LM1073', 'LM1074', 'LM1075', 'LM1076 ', 'LM1077', 'LM1078', 'LM1079', 'LM1080', 'LM1081', 'LM1082', 'LM1083', 'LM1084', 'LM1086'};
	limepth = fullfile(pth, ['LIME.', vers]);
    mkdir(limepth);
    n = 0;
    nLime = 0;
    for i = 1: numel(f) %change 1 to larger number to restart after failure
       cpth = char(deblank(f(i))); %local child path
       %fprintf('===\t%s participant %d/%d : %s\t===\n', mfilename, i, numel(f), cpth);
       %if ~isempty(strfind(cpth,'M2082')), error('all done'); end; %to stop at specific point
       if ~isempty(strfind(cpth,'_'))
          fprintf('Warning: "_" in folder name: skipping %s\n', char(cpth) );
          continue
       end
       %cpth = char(fullfile(pth,cpth)); %full child path
       inpth = fullfile(pth, cpth);
       mfile = dir(char(fullfile(inpth,'*lime.mat')));
       if isempty(mfile), continue; end;
       if numel(mfile) > 1, warning('Multiple lime.mat files in %s', inpth); end;
       outname = fullfile(outpth, [cpth, '.mat']);
       inname = fullfile(inpth, mfile(1).name);
       fprintf('%s -> %s\n', inname, outname);
       copyfile(inname,outname);
       %nii_preprocess_gui(cpth);
       n = n + 1;
       i = find(strncmp(MasterID,cpth,numel(cpth)));
       if isempty(i), continue; end;
       %fprintf('%s\n', LimeID{i(1)} );
       limename = fullfile(limepth,[LimeID{i(1)},'.mat']);
       fprintf('%s -> %s\n', inname, limename);
       copyfile(inname, limename,'f');
       nLime = nLime + 1;
    end
    fprintf('Copied %d *lime.mat files to %s\n', n, outpth);
    fprintf('Copied %d *lime.mat files to %s\n', nLime, limepth);
    if numel(f) < 2, return; end; %just touching up a single person
    if isMakeModalityTable
        nii_modality_table(outpth);
        nii_modality_table(limepth);
    end
end %if isCopymat
%end nii_preprocess_subfolders()

function nameFolds=subFolderSub(pathFolder)
d = dir(pathFolder);
isub = [d(:).isdir];
nameFolds = {d(isub).name}';
nameFolds(ismember(nameFolds,{'.','..'})) = [];
%end subFolderSub()