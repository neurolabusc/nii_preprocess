function exitCode = nii_bas(asl, t1, aslRev, inJSON, inCalScan, anatDir, dryRun, outDir)

%warning('in nii_basil');
% asl='/work/reddydp/neuro/ABC_MASTER_IN/ABC1001/ASL_ABC1001_ABC.nii';
% t1='/work/reddydp/neuro/ABC_MASTER_IN/ABC1001/T1_ABC1001_ABC.nii';
% aslRev='/work/reddydp/neuro/ABC_MASTER_IN/ABC1001/ASLrev_ABC1001_ABC.nii';
% inJSON='/work/reddydp/neuro/Neuro/nii_preprocess/ASL_DUMMY_FILES/LARC_dummy.json';
% inCalScan='';
% anatDir='';
% dryRun=false;
% outDir='/work/reddydp/neuro/ABC_MASTER_IN/ABC1001/BASIL';


%made by Ryan Joseph
%asl = .nii file representing ASL data
%t1 = structural image for that participant
%aslRev = reverse phase encoded ASL data (can be very small)
%inCalScan = calibration scan, usually the first image if a file contains of an odd number
%of volumes.
%anatDir = directory for anat processed by FSL may already exist
%dryRun = show command only, dont actually process
%overwrite = bool to do something...
%outDir = place to store all results

exitCode = 123;
if isempty(which('nii_tool')), error('nii_tool required (https://github.com/xiangruili/dicm2nii)'); end;

% TESTING
TESTING = ~exist('asl','var');

% v = version('-release');
% disp(['MATLAB VERSION: ' v]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% LINKS:
%   https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/BASIL/UserGuide
%   https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/BASIL/Tutorial
%   https://users.fmrib.ox.ac.uk/~chappell/asl_primer/data/index.html (ASL Data for Examples)

% pASL (Siemens): "PulseSequenceDetails" = "ep2d_pasl",
%   https://www.jiscmail.ac.uk/cgi-bin/webadmin?A2=ind1404&L=FSL&D=0&P=185303
%   set "labeling" to pASL
%   Bolus duration: inversion time 1
%   TIs: inversion time 2

% pcASL (Oxford): "PulseSequenceDetails" = "ep2d_VEPCASL",
%   set "labeling" to cASL/pcASL
%   Bolus duration: BolusDuration
%   PLDs: InitialPostLabelDelay
%

% pcASL (U Southern California/Danny Wang): "PulseSequenceDetails" = "ep2d_pcasl_UI_PHC",
%   set "labeling" to cASL/pcASL
%   Bolus duration: NumRFBlocks*0.0185*1000 %0.0185ms, convert to sec %https://www.mccauslandcenter.sc.edu/crnl/tools/asl
%   PLD: PostLabelDelay
%   No calibration image

% pASL (FAIREST): "PulseSequenceDetails": "ep2d_fairest_UI_iPAT"
%   http://www.pubmed.com/21606572,
%   set "labeling" to cASL/pcASL
%   Bolus duration: PostInversionDelay
%   PLD: PostLabelDelay
%   Creates 2 series: M0 Calibration Scan (Proton Density) and label+tagged pairs
%
%   TODO: SliceTiming: there is no slice timing! (use from another
%   sequence)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% GLOBALS

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SEQUENCE_PASL_SIEMENS = 'ep2d_pasl';               
SEQUENCE_PCASL_OXFORD = 'ep2d_VEPCASL';            %LARC sequence
SEQUENCE_PCASL_USC_WANG = 'ep2d_pcasl_UI_PHC';     %legacy USC pCASL
SEQUENCE_PASL_FAIREST = 'ep2d_fairest_UI_iPAT';
SEQUENCE_PCASL_TGSE = 'jw_tgse_VEPCASL';
SEQUENCE_PASL_TGSE = 'tgse_pasl';
SEQUENCE_PCASL_2D_VE11C = 'ep2d_pcasl_ve11c';

% TODO: replace if statement in readJSON with map
% SEQUENCE_ALIASES = containers.Map({'to_ep2d_VEPCASL'},[SEQUENCE_PCASL_OXFORD]);

% Grouping order:
%   repeats: ibf=tis
%   label/control pairs: ibf=rpt
GROUPING_ORDER_REPEATS = 'tis';
GROUPING_ORDER_LABEL_CONTROL_PAIRS = 'rpt';

% Label/Control pairs: --iaf
%   label then control: (iaf=tc)
%   control then label: (iaf=ct)
PAIRS_LABEL_THEN_CONTROL = 'tc';
PAIRS_CONTROL_THEN_LABEL = 'ct';

% TODO:
% - are these hard coded values correct for all sequences?
% '--t1b 1.65 ',           ...
% '--alpha 0.85 ',         ...

% EXAMPLE SEQUNCES FOR TESTING
if TESTING
    disp ' ';
    disp '========= TESTING =========';
    disp ' ';
    
    % chris (SEQUENCE_PCASL_OXFORD)
    % sequence = SEQUENCE_PCASL_OXFORD;
    % outDir = '/Volumes/Ryan/aslRorden/out';
    % dir = '/Volumes/Ryan/aslRorden';
    % asl = [dir '/24_asl_ep2d_PCASL.nii'];
    % t1 = [dir '/6_anat_T1w.nii'];
    
    % william (SEQUENCE_PCASL_OXFORD)
    % outDir = '/Volumes/Ryan/aslRorden/PCASL_20190208135617';
    % dir = outDir;
    % asl = [dir '/PCASL_20190208135617_to_ep2d_PCASL_20190208135617_11.nii.gz'];
    % aslRev = [dir '/PCASL_20190208135617_to_ep2d_PCASL_PA_20190208135617_15.nii.gz'];
    % t1 = [dir '/PCASL_20190208135617_T1_mprage_ns_sag_p2_iso_1mm_192_20190208135617_2.nii.gz'];
    
    % SEQUENCE_PCASL_USC_WANG
    % outDir = '/Volumes/Ryan/aslRorden/pasl_pcasl_old/pcasl/20120607135306';
    % dir = outDir;
    % asl = [dir '/20120607135306_pCASL_1200ms_3.5sTR_10.nii.gz'];
    % t1 = [dir '/20120607135306_T1_mprage_ns_sag_p2_iso_1.0mm_192_7.nii.gz'];
    
    % SEQUENCE_PASL_SIEMENS
    % outDir = '/Volumes/Ryan/aslRorden/pasl_pcasl_old/pasl';
    % dir = outDir;
    % asl = [dir '/4_ep2d_tra_pasl_p2.nii'];
    % t1 = [dir '/2_Anatomical_1x1x1_ipat=2_nw.nii'];
    
    % SEQUENCE_PASL_FAIREST
    % outDir = '/Volumes/Ryan/aslRorden/pasl_fairest_new';
    % dir = outDir;
    % asl = [dir '/20091031130655_PASL_fairest_16x5_no-fs1.5_7.nii.gz'];
    % t1 = [dir '/20091031130655_Anatomical_1x1x1_ipat=2_nw_4.nii.gz'];
    % inCalScan = [dir '/20091031130655_PASL_fairest_16x5_M0_nofs_5.nii.gz'];
    
    % sen pcasl
    % outDir = '/Volumes/Ryan/aslRorden/sen';
    % dir = outDir;
    % asl = [dir '/16_ep2d_pcasl_UI_PHC_p4.nii'];
    % t1 = [dir '/4_T1_MPRage.nii'];
    
    %outDir = '/home/coffee/Desktop/RyGuy/pasl_3d_11_MATLAB';
    %sourceDir = '/home/coffee/Desktop/dcm_qa_asl/Ref';
    %asl = [sourceDir '/Sen.nii'];
    %aslRev = [];
    %inJSON = '';
    %t1 = '';
   % inCalScan = [sourceDir 'pasl_3d_m0_13.nii'];
    %inCalScan = [];
    %aslOutput = [outDir '/oxford_asl_output'];
    %anatDir = '';
    %dryRun = true;
    
else
    % TODO: is this ok for m0 files that are made?
    aslOutput = outDir;
end

if dryRun
    disp '========= DRY RUN =========';
    disp ' ';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('aslRev','var')
    aslRev = [];
end

% verify inputs
if ~isempty(t1) && ~exist(t1,'file'), error('Requires T1 scan %s\n', t1); end

if ~exist(asl,'file'), error('Requires actual asl scan  %s\n', asl); end
if ~isempty(aslRev) && ~exist(aslRev,'file'), error('Requires reverse phase encoded file %s\n', aslRev); end

if ~exist('outDir','var') outDir = pwd; end
if ~exist(outDir,'dir') mkdir(outDir); end
if ~exist(aslOutput,'dir') mkdir(aslOutput); end

% will create:
aslX = [];        %asl data with M0 removed
aslRevX = [];     %reverse phased first volulme
m0X = [];         %M0 calibration scan
if exist(aslRev,'file')
    nii = nii_tool('load', aslRev);
    [~,n,x] = nii_fileparts(aslRev);
    niiRev.hdr = nii.hdr;    % outDir = '/Volumes/Ryan/aslRorden/sen';

    niiRev.hdr.dim(1) = 3; %3d not 4D data
    niiRev.hdr.dim(5) = 1; %single volume
    aslRevX = fullfile(outDir, ['m0', n, x]);
    niiRev.img = nii.img(:,:,:,1);
    nii_tool('save', niiRev, aslRevX);
end

% load ASL
if ~exist(asl,'file'), error('Unable to find %s\n', asl); end
nii = nii_tool('load', asl);
[TR, PLDs, BD, SecPerSlice, EffectiveEchoSpacing, sequence] = readJson(asl, inJSON);

% check bolus duration is within a rangeread
if BD > 20
    warning('Abnormally high Bolus Duration: %s detected, dividing by 1000 assuming it was entered as ms...',BD);
    BD = BD /1000;
end

% additonal steps for scans with calibration images
switch sequence
    case {SEQUENCE_PASL_SIEMENS, SEQUENCE_PCASL_OXFORD}
        ok = check_volumesOK(nii, true);
        if ~ok, warning('Unable to process ASL dataset: wrong number of volumes'); return; end
        %save calibration image
        [~,n,x] = nii_fileparts(asl);
        niiM0.hdr = nii.hdr;
        niiM0.hdr.dim(1) = 3; %3d not 4D data
        niiM0.hdr.dim(5) = 1; %single volume
        m0X = fullfile(outDir, ['m0_', n, x]);
        niiM0.img = nii.img(:,:,:,1);
        nii_tool('save', niiM0, m0X);
        %save subsequent volumes (ASL images)
        [~,n,x] = nii_fileparts(asl);
        niiAsl.hdr = nii.hdr;
        niiAsl.hdr.dim(5) = nii.hdr.dim(5) - 1; %single volume
        aslX = fullfile(outDir, ['asl_', n, x]);
        niiAsl.img = nii.img(:,:,:,2:end);
        nii_tool('save', niiAsl, aslX);
    case {SEQUENCE_PCASL_USC_WANG,SEQUENCE_PCASL_2D_VE11C}
        check_volumes(nii, false)
        aslX = asl;
    case {SEQUENCE_PCASL_TGSE, SEQUENCE_PASL_TGSE}
        check_volumes(nii, false)
        aslX = asl;
        m0X = inCalScan;
    case SEQUENCE_PASL_FAIREST
        check_volumes(nii, false)
        aslX = asl;
        if ~exist(inCalScan, 'file')
            warning('SEQUENCE_PASL_FAIREST requires a calibration image to be passed in directly');
            return;
        end
        m0X = inCalScan;

    otherwise
        warning('invalid sequence %s', sequence)
        return;
        
end

% sequence specific options
switch sequence
    case SEQUENCE_PASL_SIEMENS
        options.cASL = false;
        options.groupingOrder = GROUPING_ORDER_LABEL_CONTROL_PAIRS;
        options.labelControlPairs = PAIRS_LABEL_THEN_CONTROL;
    case {SEQUENCE_PCASL_OXFORD, SEQUENCE_PCASL_USC_WANG,SEQUENCE_PCASL_2D_VE11C}
        options.cASL = true;
        options.groupingOrder = GROUPING_ORDER_LABEL_CONTROL_PAIRS;
        options.labelControlPairs = PAIRS_LABEL_THEN_CONTROL;
    case SEQUENCE_PASL_FAIREST
        options.cASL = false;
        % TODO: what grouping order is FAIREST? chris seems to think
        % label/control pairs
        options.groupingOrder = GROUPING_ORDER_REPEATS;
        options.labelControlPairs = PAIRS_LABEL_THEN_CONTROL;
    case SEQUENCE_PCASL_TGSE
        options.cASL = true;
        options.groupingOrder = GROUPING_ORDER_REPEATS;
        options.labelControlPairs = PAIRS_LABEL_THEN_CONTROL;
     case SEQUENCE_PASL_TGSE
        options.cASL = false;
        options.groupingOrder = GROUPING_ORDER_LABEL_CONTROL_PAIRS;
        options.labelControlPairs = PAIRS_CONTROL_THEN_LABEL;
end


controlBeforeLabel = aslOrder(aslX);
if ~controlBeforeLabel && strcmp(options.labelControlPairs,PAIRS_CONTROL_THEN_LABEL)
	controlBeforeLabel = aslOrder(aslX, true);
    warning('ASL has label-then-control order, but specified as control-then-label')
    return;
end
if controlBeforeLabel && strcmp(options.labelControlPairs,PAIRS_LABEL_THEN_CONTROL)
	controlBeforeLabel = aslOrder(aslX, true);
	warning('ASL has control-then-label order, but specified as label-then-control')
    return;
end
options.output = aslOutput;
options.inputImage = aslX;
options.calibrationImage = m0X;
options.reversePhaseEncodedImage = aslRevX;
options.T1Image = t1;
options.TR = TR;
options.PLDs = PLDs;
options.BD = BD;
options.secPerSlice = SecPerSlice;
options.effectiveEchoSpacing = EffectiveEchoSpacing;
options.partialVolumeCorrection = true;

if ~isempty(options.T1Image)
    % generate anat .struct file
    if isempty(anatDir)
        anatDir = options.output;
    else
        [p,n,x] = fileparts(anatDir);
        if strcmpi([n,x], 'struc.anat')
            anatDir = p;
        end
    end
    options.anatFile = [anatDir '/struc.anat'];
    if ~exist(anatDir, 'dir')
        error('Missing anatDir named "%s"\', anatDir);
    end

    if ~exist(options.anatFile, 'dir')
        command = ['fsl_anat -i "' options.T1Image '" -o ' anatDir '/struc'];
        if ~dryRun
            exitCode = fslCmd(command);
            if exitCode ~= 0
                return;
            end
        else
            % show command for dry runs
            disp(command);
        end
    else %exist(options.anatFile, 'dir')
        if false
            fprintf('Using existing files %s\n', options.anatFile); 
        else
            %20200325 RNN and DR suggest we always recreated files
            command = ['fsl_anat -i "' options.T1Image '" --clobber -o ' anatDir '/struc'];
            if ~dryRun
                exitCode = fslCmd(command);
                if exitCode ~= 0
                    return;
                end
            else
                % show command for dry runs
                disp(command);
            end
        end
    end
else
    options.anatFile = '';
    options.partialVolumeCorrection = false;
end

command = makeCommand(options);

% run the command
if ~isempty(command) && ~dryRun
    exitCode = fslCmd(command);
else
    disp(command);
    exitCode = 0;
end


% Local functions follow ============================================================================================================
% check number of volumes
    function check_volumes(nii, expectsCalibImage)
        vols = nii.hdr.dim(5);
        if expectsCalibImage && (mod(vols, 2) == 0)
            error('You said it has a calibration image, but there is an even number of volumes');
        elseif ~expectsCalibImage && (mod(vols, 2) ~= 0)
            error('You said it doesn''t have a calibration image, but there is an odd number of volumes');
        end
    end

    function ok = check_volumesOK(nii, expectsCalibImage)
        vols = nii.hdr.dim(5);
        ok = true;
        if expectsCalibImage && (mod(vols, 2) == 0)
            ok = false;
        elseif ~expectsCalibImage && (mod(vols, 2) ~= 0)
            ok = false;
        end
    end

    function command = makeCommand(options)
        
        % generate tis
        tis = '';
        for n = 1:numel(options.PLDs)
            % TODO: is this true for all sequences?
            label = options.BD + options.PLDs(n);
            tis = strcat(tis, num2str(label));
            if n < 6
                tis = strcat(tis, ',');
            end
        end
        
        % program
        part1 = 'oxford_asl';
        
        % input data
        % tis = PostLabelDelay+BolusTime
        % casl = labeling cASL/pcASL
        if options.cASL
            cASLFlag = '--casl';
        else
            cASLFlag = '';
        end
        
        if options.secPerSlice
            slicedtFlag = ['--slicedt ' num2str(options.secPerSlice)];
        else
            slicedtFlag = '';
        end
        
        part2 = [
            '-i "', options.inputImage, '" ',             ...  % input image
            cASLFlag, ' ',                                ...  % Labeling (pASL or cASL/pcASL)
            '--ibf=', options.groupingOrder, ' ',         ...  % (Grouping order): Input block format. Specifically for multi-delay (multi-PLD) ASL data to identify whther the individual delays/PLDs are groups togther or by repeats of the same sequence of PLDs.
            '--iaf=', options.labelControlPairs, ' ',     ...  % (Label/Control pairs): Input ASL format: specifies if the data has already been label-control subtracted (diff, default), or is in the form of label(tag)-control pairs (tc or ct depending upon if label/tag is first).
            '--tis ', tis, ' ',                           ...  % TIs/PLDs depending on labeling (pASL/pcASL)
            '--bolus ', num2str(options.BD), ' ',         ...  % Bolus duration
            slicedtFlag                                        % Time per slice (ms)
            ];
        
        % structural
        if isempty(options.anatFile)
            part3 = []; %no structural file
        else
            part3 = ['--fslanat="', options.anatFile, '"'];
        end
        
        % calibration
        if isempty(options.calibrationImage)
            part4 = []; %no calibration image
        else
            part4 = [
                '-c "', options.calibrationImage, '" ',   ...
                '--tr ', num2str(options.TR), ' ',        ...
                '--cgain 1.00 ',                          ...
                '--cmethod voxel'
                ];
        end
        
        % distortion correction
        if isempty(options.reversePhaseEncodedImage)
            part5 = []; %no phase reversed image
        else
            part5 = [
                '--cblip="',  options.reversePhaseEncodedImage ,'" ',                   ...
                '--echospacing=', num2str(options.effectiveEchoSpacing * 1000) ,' '     ...
                '--pedir=y'
                ];
        end
        
        % analysis
        if options.partialVolumeCorrection
            pvcorrFlag = '--pvcorr ';
        else
            pvcorrFlag = '';
        end
        
        part6 = [
            '--wp ',                 ...  % analysis conforms to white paper
            '--t1b 1.65 ',           ...
            '--alpha 0.85 ',         ...  % inversion efficiency
            '--fixbolus ',           ...  % fixed label duration
            '--spatial ',            ...  % adaptive spatial regularization
            '--artoff ',             ...
            pvcorrFlag,              ...  % partial volume correction
            '--mc ',                 ...  % motion correction
            '-o "', options.output, '"'
            ];
        
        command = [part1 ' ' part2 ' ' part3 ' ' part4 ' ' part5 ' ' part6];
    end

    function SecPerSlice = calculateSliceTiming(SliceTiming)
        if iscell(SliceTiming)
            SecPerSlice = (max(cell2mat(SliceTiming(:))) -min(cell2mat(SliceTiming(:))))/(numel(SliceTiming)-1);
            if (  min(cell2mat(SliceTiming(:))) ~= 0)
                warning('Warning: Non-Ascending Slice Timing Detected');
            end
        else
            SecPerSlice = (max(SliceTiming(:)) -min(SliceTiming(:)))/(numel(SliceTiming)-1);
            if (  min(SliceTiming(:)) ~= 0)
                warning('Warning: Non-Ascending Slice Timing Detected');
            end
        end
    end      

    function [TR, PLDs, BD, SecPerSlice, EffectiveEchoSpacing, sequence] = readJson(niiName, json)
        [p,n,x] = nii_fileparts(niiName);
        % if not json file is provided then infer from .nii file name
        if isempty(json)
            json = fullfile(p,[n,'.json']);
        end
        if ~exist(json,'file'), error('Unable to find %s\n', json); end
        js = loadJson(json);
        if isfield(js, 'MultibandAccelerationFactor'), error('Update to support multiband'); end
        if ~isfield(js, 'PulseSequenceDetails'), error('Require PulseSequenceDetails'); end
        % get sequence name
        result = regexp(js.PulseSequenceDetails, '^%(CustomerSeq|SiemensSeq)+%_(\w+)$', 'tokens');
        if isempty(result)
           result = regexp(js.PulseSequenceDetails, '^%(CustomerSeq|SiemensSeq)+%\\(\w+)$', 'tokens'); 
        end
        if isempty(result)  
            js.PulseSequenceDetails
            error('Cant parse pulse sequence details'); 
        
        end
        prefix = result{1}{1};
        sequence = result{1}{2};
        %sequence = js.PulseSequenceDetails;
        % TODO: make reverse lookup table 
        
        
        if strcmpi(sequence,'to_ep2d_VEPCASL')
            sequence = SEQUENCE_PCASL_OXFORD;
        end
        if strcmpi(sequence,'tgse_pcasl_ve11c') %aka 'jw_tgse_VEPCASL'
            sequence = SEQUENCE_PCASL_TGSE;
        end
        switch sequence
           
            case SEQUENCE_PCASL_OXFORD
                %if ~isfield(js, 'InversionTime'), error('Require InversionTime'); end
                if ~isfield(js, 'BolusDuration'), error('Require BolusDuration'); end
                if ~isfield(js, 'InitialPostLabelDelay'), error('Require InitialPostLabelDelay'); end
                if ~isfield(js, 'RepetitionTime'), error('Require RepetitionTime'); end
                if ~isfield(js, 'SliceTiming'), error('Require SliceTiming'); end
                if ~isfield(js, 'EffectiveEchoSpacing'), error('Require EffectiveEchoSpacing'); end
                SecPerSlice = calculateSliceTiming(js.SliceTiming);
                
                BD = js.BolusDuration;
                PLDs = js.InitialPostLabelDelay;
                TR = js.RepetitionTime;
                EffectiveEchoSpacing = js.EffectiveEchoSpacing;

            case SEQUENCE_PCASL_TGSE
                %if ~isfield(js, 'InversionTime'), error('Require InversionTime'); end
                if ~isfield(js, 'InitialPostLabelDelay'), error('Require InitialPostLabelDelay'); end
                if ~isfield(js, 'RepetitionTime'), error('Require RepetitionTime'); end
                if isfield(js, 'MaximumT1Opt')
                    BD = js.MaximumT1Opt;
                elseif isfield(js, 'T1')
                    BD = T1;
                else
                    error('Require MaximumT1Opt');
                end
                PLDs = js.InitialPostLabelDelay;
                TR = js.RepetitionTime;
                SecPerSlice = '';
                EffectiveEchoSpacing = '';                  

                
            case SEQUENCE_PASL_SIEMENS
                if ~isfield(js, 'InversionTime') && isfield(js, 'RepetitionTime') && (abs(2.5 - js.RepetitionTime) < 0.05) && isfield(js, 'SliceTiming') && (numel(js.SliceTiming) == 14)
                    js.InversionTime = 0.8;
                    warning('Unable to read ASL values from JSON: guessing');
                    %error('x')
                end

                if ~isfield(js, 'InversionTime'), error('Require InversionTime'); end
                if ~isfield(js, 'RepetitionTime'), error('Require RepetitionTime'); end
                if ~isfield(js, 'SliceTiming'), error('Require SliceTiming'); end
                if ~isfield(js, 'EffectiveEchoSpacing'), error('Require EffectiveEchoSpacing'); end
                SecPerSlice = calculateSliceTiming(js.SliceTiming);
                
                BD = js.InversionTime;
                % TODO: what should this be??
                PLDs = [1];
                TR = js.RepetitionTime;
                EffectiveEchoSpacing = js.EffectiveEchoSpacing;
                  
            case {SEQUENCE_PCASL_USC_WANG,SEQUENCE_PCASL_2D_VE11C}
                if ~isfield(js, 'NumRFBlocks') && isfield(js, 'RepetitionTime') && (abs(3.5 - js.RepetitionTime) < 0.05) && isfield(js, 'SliceTiming') && (numel(js.SliceTiming) == 17)
                    js.NumRFBlocks = 80;
                    js.PostLabelDelay = 1.2;
                    warning('Unable to read pCASL values from JSON: guessing');
                end
                if ~isfield(js, 'NumRFBlocks'),
                    error('Require NumRFBlocks'); 
                end
                if ~isfield(js, 'PostLabelDelay'), error('Require PostLabelDelay'); end
                if ~isfield(js, 'RepetitionTime'), error('Require RepetitionTime'); end
                if ~isfield(js, 'SliceTiming'), error('Require SliceTiming'); end
                if ~isfield(js, 'EffectiveEchoSpacing'), error('Require EffectiveEchoSpacing'); end
                SecPerSlice = calculateSliceTiming(js.SliceTiming);
                
                BD = js.NumRFBlocks * 0.0185;
                PLDs = [js.PostLabelDelay];
                TR = js.RepetitionTime;
                EffectiveEchoSpacing = js.EffectiveEchoSpacing;
                
            case SEQUENCE_PASL_FAIREST
                if ~isfield(js, 'PostInversionDelay'), error('Require PostInversionDelay'); end
                if ~isfield(js, 'PostLabelDelay'), error('Require PostLabelDelay'); end
                if ~isfield(js, 'RepetitionTime'), error('Require RepetitionTime'); end
                if ~isfield(js, 'EffectiveEchoSpacing'), error('Require EffectiveEchoSpacing'); end
                %   TODO: SliceTiming is not available, hardcode something meaningful
                SecPerSlice = 1;
                BD = js.PostInversionDelay;
                PLDs = [js.PostLabelDelay];
                TR = js.RepetitionTime;
                EffectiveEchoSpacing = js.EffectiveEchoSpacing;
            
            case SEQUENCE_PASL_TGSE
                if ~isfield(js, 'BolusDuration'), error('Require BolusDuration'); end
                if ~isfield(js, 'InversionTime'), error('Require InversionTime'); end
                if ~isfield(js, 'RepetitionTime'), error('Require RepetitionTime'); end
                BD = js.BolusDuration;
                PLDs = [js.InversionTime];
                TR = js.RepetitionTime;
                SecPerSlice = '';
                EffectiveEchoSpacing = '';    
                
            otherwise
                error('Invalid sequence %s', sequence);
        end
    end

    function [json, fnm] = loadJson(fnm)
        p = fileparts(which(mfilename));
        if ~exist('fnm','var'), fnm = fullfile(p, 'config.json'); end
        if ~exist(fnm,'file')
            fprintf('Creating empty JSON: Unable to find %s\n', fnm);
            json = [];
            return;
        end
        fid = fopen(fnm);
        raw = fread(fid,inf);
        str = char(raw');
        fclose(fid);
        if exist('jsondecode', 'builtin')
            json = jsondecode(str);
        else
            % https://www.mathworks.com/matlabcentral/fileexchange/20565-json-parser
            js = parse_json(str);
            json = js{1};
        end
    end

    function [pth,nam,ext,num] = nii_fileparts(fname)
        % extends John Ashburner's spm_fileparts.m to include '.nii.gz' as ext
        num = '';
        if ~ispc, fname = strrep(fname,'\',filesep); end
        [pth,nam,ext] = fileparts(fname);
        ind = find(ext==',');
        if ~isempty(ind)
            num = ext(ind(1):end);
            ext = ext(1:(ind(1)-1));
        end
        if strcmpi(ext,'.gz')
           [pth nam ext] = fileparts(fullfile(pth, nam));
           ext = [ext, '.gz'];
        end
    end

    function path = findBestPath(paths)
        path = nan;
        for i = 1:numel(paths)
            if exist(paths{i},'dir')
                path =  paths{i};
                return
            end
        end
    end

    function exitCode = fslCmd(Cmd)
        % find fsl path
        fsldir = findBestPath({ '/usr/local/fsl/',
            '/Volumes/Ryan/fsl/',
            '/usr/share/fsl/6.0/',
            getenv('FSLDIR')
            });
        if ~exist(fsldir,'dir')
            error('%s: fsldir not found',mfilename);
        end
        setenv('FSLDIR', fsldir);
        flirt = [fsldir '/bin/flirt'];
        if ~exist(flirt,'file')
            error('%s: flirt (%s) not found',mfilename,flirt);
        end
        opts = '';
        global noGz % <- your need this!
        if ~isempty(noGz) && noGz, opts = 'export FSLOUTPUTTYPE=NIFTI'; end
        command=sprintf('sh -c ". %s/etc/fslconf/fsl.sh; %s %s/bin/%s"\n',fsldir, opts, fsldir, Cmd);
        %command=sprintf('sh -c ". %s/etc/fslconf/fsl.sh; %s/bin/%s" >log.txt\n',fsldir,fsldir, Cmd);
        
        fprintf(command);
        exitCode = system(command);
    end

    function controlBeforeLabel = aslOrder(fnm, verbose)
        %determine whether label-control order is lclclc or clclcl
        % fnm : name of NIfTI image to inspect
        % verbose : report text 
        %return
        % controlBeforeLabel == true: order is clclcl
        % controlBeforeLabel == false: order is lclcl

        %Example
        % aslOrder('pasl_3d_11.nii')
        % aslOrder %use GUI
        controlBeforeLabel = false;
        if ~exist('fnm', 'var')
           [n,p] = uigetfile('*.nii;*nii.gz'); 
           fnm = fullfile(p,n);
        end
        if ~exist('verbose', 'var')
           verbose = false; 
        end
        if ~exist(fnm, 'file')
           error('unable to find %s', fnm); 
        end
        img = nii_tool('img', fnm);
        if size(img,4) < 3
           error('asl files should have at least 3 volumes: %s', fnm);
        end
        if verbose
            fprintf('%s report for %s\n', mfilename, fnm);
        end
        v1 = 1;
        if mod(size(img,4),2) ~= 0 %odd number of volumes!
            warning('odd number of volumes, assuming first volume is M0 image that will be discarded: %s', fnm);
            v1 = 2; %skip first volume
            %determine if
            % mlclclc
            % mclclcl
        end
        sum = zeros(2,1);
        for i = v1 : size(img,4)
            im = img(:,:,:,i);
            sum(2 - mod(i,2)) = sum(2 - mod(i,2)) + mean(im(:));
            %sum(2 - mod(i,2)) = sum(2 - mod(i,2)) + median(im(:));
            if verbose 
                fprintf('%d %d %g %g\n', i, 2 - mod(i,2), mean(im(:)), median(im(:)) );
            end
        end
        %control scans are brighter than labelled scans
        controlBeforeLabel = (sum(1) > sum(2)); 
        prefix = '';
        if v1 == 2
            controlBeforeLabel = ~controlBeforeLabel;
            prefix = 'm';
        end
        if ~verbose, return; end
        fprintf('volumes = %d sum odd = %g sum even = %g\n', size(img,4), sum(1), sum(2));
        if controlBeforeLabel
            fprintf('order %sclclcl...\n', prefix); 
        else
            fprintf('order %slclclc...\n', prefix);
        end
    end    
    
end

