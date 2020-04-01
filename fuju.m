function fugu
pth = fileparts(which(mfilename))
numVolASLImg = 74;
jsonFile = 'xx';


if ~exist(jsonFile,'file')
    warning('Guessing ASL parameters: missing JSON "%s"\n', jsonFile)
    p = fileparts(which(mfilename));
    p = fullfile(p,'json');
    switch numVolASLImg

        case 101 %handles POLAR and legacy pasl sequences
            jsonFile = fullfile(p,'POLAR_dummy.json');
        case 74
            jsonFile = fullfile(p,'LEGACY_PCASL_dummy.json');
        case 97
            jsonFile = fullfile(p,'LARC_dummy.json'); 
            [pth,nam,ext] = spm_fileparts(imgs.ASL);
            ASLrev = imgs.ASLrev;   % [pth '/' nam 'rev',ext];
        case 60
            jsonFile = fullfile(p,'SEN_dummy.json'); 
        otherwise
            error('Can''t guess sequence type');
    end
    if ~exist(jsonFile,'file')
        error('Can''t find json %s', jsonFile);
    end
    
end
jsonFile