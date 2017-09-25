function imgs = nii_preprocess_gui(loadprev, ignoreImgPaths)
%Allow user to select files for preprocessing
% loadprev : (optional) if 1 then user is prompted to select existing mat file
%                       if string, that mat file is loaded
%                       otherwise, user creates mat file by selecting images
%
% ignoreImgPaths : (optional) if true images MUST exist in mat file's folder
%Examples
% nii_preprocess_gui %use gui
% nii_preprocess_gui(1) %use gui to select previous mat file
% nii_preprocess_gui('T1_M2094_limegui.mat') %specify previous mat file
% nii_preprocess_gui(pwd) %open mat file if it exists or auto-create new one
if nargin > 0
    if isdir(loadprev) %if user passes folder, than either process limegui file or auto-generate and process limegui file
        mfile = dir(char(fullfile(loadprev,'*limegui.mat')));
        if isempty(mfile) %limegui does not exist: create one
            img = getfiles_x (loadprev);
            if isempty(img.T1)
               warning('Unable to find T1 scan in folder %s\n', loadprev);
               return;
            end
            [p,n] = fileparts(img.T1);
            loadprev = fullfile(loadprev, [n, '_limegui.mat']);
            save(loadprev,'-struct', 'img');
        else
            loadprev = char(fullfile(loadprev, mfile(1).name));
            if numel(mfile) >1
                fprintf('Warning: multiple limegui files in folder. Assume you wanted %s\n', loadprev);
            end
        end
    end %isdir
    if ischar(loadprev) && exist(loadprev, 'file')
        [p, n, x] = fileparts(loadprev);
        f = [n, x];
    else
        [f, p] = uigetfile('*limegui.mat', 'Choose name mat file');
    end
    matName = fullfile(p,f);
    imgs = load(matName);
    if exist('ignoreImgPaths', 'var') && ignoreImgPaths
        imgs = deletePathSub (imgs);
    end
    if ~isfield(imgs,'T1') || isempty(imgs.T1)
        fprintf('%s error: no T1 image for %s\n', mfilename, matName);
        return;
    end
    imgs = checkPathSub (imgs,matName);
    imgs = stripGzSub (imgs);
    nii_preprocess(imgs);
    return;
end
%Graphical interface for running nii_preprocess
imgs = [];
go = false;
Xpressed = false; %user closes window without pressing button
while ~go && ~Xpressed
    [img, go, Xpressed] = getfiles_x;
    if isempty(img)
        fprintf('Items cleared\n');
    elseif isempty(img.T1)
        fprintf('Error: T1 scan required\n');
    else
        imgs = [imgs, img]; %#ok<AGROW>
    end;
end
if isempty(imgs), fprintf('Aborting: no images selected'); return; end;
if Xpressed, fprintf('Aborting: user never pressed go'); return; end;
for i = 1: numel(imgs)
    img = imgs(i);
    %save filenames that will be processed -
    % to rerun in future: img = load('name_limegui.mat'); nii_preprocess(img);
    [p,n] = fileparts(img.T1);
    matName = fullfile(p, [n, '_limegui.mat']);
    save(matName,'-struct', 'img');
    %process the data
    nii_preprocess(img);
end
end %nii_preprocess_gui()

function [img, go, Xpressed] = getfiles_x (autoLoadDir)
%if autoLoadDir is a folder, then files automatically filled
go = false;
Xpressed = false;
alreadyWarned = false;
modalities = {'T1', 'T2', 'ASL', 'fMRI', 'Rest', 'DTI', 'DTIrev', 'Lesion'};
modalitiesVerbose = {'T1', 'T2', 'ASL', 'fMRI', 'Resting State', 'DTI', 'DTI (reversed phase)', 'Lesion'};
modalitiesDTI = [false, false, false, false, false, true, true, false];
for i = 1: numel(modalities)
	img.(char(modalities(i))) = [];
end;
if exist('autoLoadDir','var') && isdir(autoLoadDir)
    autoload(autoLoadDir)
    return;
end
btns = zeros(numel(modalities),1);
renameCheck = 0;
btnID = 0;
d = nestedfx;
uiwait(d);

    function d = nestedfx
      d = dialog('Position',[300 300 620 360],'Name','image selection','CloseRequestFcn',@close_callback);
      topMostPos = 330;
      pos = 10:30:topMostPos;
      btnID = uicontrol('Parent',d,...
                   'Style','Edit',...
                   'Position',[10 pos(end) 600 25],...
                   'String','Participant ID');
      btnTxt = uicontrol('Parent',d,...
                   'Style','text',...
                   'Position',[10 pos(end-1) 600 25],...
                   'String','Select the desired scans');
        for i = 1: numel(modalities)
            btns(i) = uicontrol('Parent',d,...
                    'Position',[10 pos(end- i - 1) 600 25],...
                    'String', char(modalitiesVerbose(i)),...
                    'Callback',{@filename_callback,char(modalities(i)), i, modalitiesDTI(i)});
        end;
        btnGo = uicontrol('Parent',d,...
                   'Position',[10 10 120 25],...
                   'String','Go',...
                   'Callback',@go_callback);
        btnNext = uicontrol('Parent',d,...
                   'Position',[140 10 120 25],...
                   'String','Save',...
                   'Callback',@save_callback);
        btnLoad  = uicontrol('Parent',d,...
                   'Position',[270 10 120 25],...
                   'String','Load',...
                   'Callback',@load_callback);
        btnAuto  = uicontrol('Parent',d,...
                   'Position',[400 10 120 25],...
                   'String','Auto',...
                   'Callback',@auto_callback);
        renameCheck = uicontrol('Parent',d,...
                   'style','checkbox',...
                   'units','pixels',...
                   'position',[520 15 100 15],...
                   'string','Rename files');

    end
    function go_callback(~,~,~)
        go = true;
        delete(gcf);
    end
    function autoload(imgDir)
        fnms = dir(fullfile(imgDir,'*.nii'));
        if numel(fnms) < 1, disp('Unable to find .nii images'); return; end;

        for i = 1: numel(modalities)
            m = char(modalities(i));
            img.(m) = [];
            modalityKey = [m,'_'];
            for f = 1 : numel(fnms)
                if strncmpi(fnms(f).name,modalityKey,numel(modalityKey))
                    img.(m) = fullfile(imgDir, fnms(f).name);
                end
            end
        end
        
    end %autoload()
    function auto_callback(~, ~, ~) %filename_callback(obj,callbackdata,myField)
        title = 'Select folder images (T1_*.nii, ASL_*.nii, etc)\n';
        fprintf(title);
        imgDir = uigetdir(pwd, title);
        autoload(imgDir);
        for i = 1: numel(modalities)
            m = char(modalities(i));
            if isempty(img.(m))
                set(btns(i),'string', char(modalitiesVerbose(i)));
            else
                set(btns(i),'string', [m, ': ',img.(m)]);
            end
        end
        [p,n] = fileparts(imgDir);
        set(btnID,'String',n);
    end %end auto_callback()
    function load_callback(~, ~, ~) %filename_callback(obj,callbackdata,myField)
        [f, p] = uigetfile('*limegui.mat', 'Select a mat file');
        matName = fullfile(p,f);
        imgNew = load(matName);
        if ~isfield(imgNew,'T1'), fprintf('Not a valid limegui.mat file (no T1): %s\n', matName); return; end;
        if ~ischar(imgNew.T1), fprintf('Not a valid limegui.mat file (T1 is not a file name): %s\n', matName); return; end;
        imgNew = checkPathSub (imgNew,matName);
        for i = 1: numel(modalities)
            m = char(modalities(i));
            img.(m) = [];
            if (isfield(imgNew, m)) && (~isempty(imgNew.(m)))
                img.(m) = imgNew.(m);
            end
            m
            if isempty(img.(m))
               set(btns(i),'string', char(modalitiesVerbose(i)));
            else
                set(btns(i),'string', [m, ': ',img.(m)]);
            end;
        end
    end %end load_callback()
    function save_callback(~, ~, ~) %filename_callback(obj,callbackdata,myField)
        if isempty(img.T1)
           warning('Please specify T1 scan');
           return;
        end
        [p,n] = fileparts(img.T1);
        matName = fullfile(p, [n, '_limegui.mat']);
        %'*limegui.mat'
        [f,p] = uiputfile(matName,'Save file name');
        matName = fullfile(p,f);
        save(matName,'-struct', 'img');
        img = [];
        delete(gcf); %create blank document
    end %end next_callback()
    function close_callback(~,~,~)
        delete(gcf);
        Xpressed = true;
    end
%     function check_callback(hObject, ~, ~) %filename_callback(obj,callbackdata,myField)
%         checkValue = get(hObject,'Value');
%         disp(checkValue);
%     end %end check_callback()
    function filename_callback(obj,~,myField, myIndex, isDti) %filename_callback(obj,callbackdata,myField)
        dorename = false;
        checkValue = get(renameCheck,'Value');
        if checkValue > 0
            if ~alreadyWarned
                m = msgbox('Warning: Your image will be renamed to a shortened version','Rename Warning');
                alreadyWarned = true;
                uiwait(m);
            end
            dorename = true;
        end
        if exist('isDti','var') && isDti
            [A,Apth] = uigetfile({'*.bvec';'*.*'},['Select ', myField, ' vectors']);
            if isnumeric(A), return; end; %user pressed cancel
            [Apth,A,ext] = fileparts(fullfile(Apth,A));
            imgName = fullfile(Apth,[A '.nii']);
            if exist(imgName, 'file'),
                ext = '.nii';
            else
                A = [A, '.nii'];
                ext = '.gz';
                imgName = fullfile(Apth,[A, ext]);
                if ~exist(imgName, 'file')
                    error('Unable to find %s', imgName);
                end
            end;

        else
            [A,Apth] = uigetfile({'*.nii;*.gz;*.hdr;';'*.*'},['Select ', myField, ' image']);
            if isnumeric(A), return; end;
            [Apth,A,ext] = fileparts(fullfile(Apth,A));
        end;

        if strfind(A,'nii') %if .nii still exists in filename 'A'
            newext = ['.nii' ext];
            realname = A(1:end-4);
        else
            newext = ext;
            realname = A;
        end

        if dorename
            ID = get(btnID,'String');
            if strcmpi(ID,'Participant ID') %user did not provide ID so don't use the field value
                ID = '';
            end
            Ashort = [myField '_' ID];
            copyfile(fullfile(Apth,[A ext]),fullfile(Apth,[Ashort newext]));
            img.(myField) = fullfile(Apth,[Ashort newext]);
            set(btns(myIndex),'string', [myField, ': ',fullfile(Apth,[Ashort newext])]);
            %if myIndex == 6 || myIndex == 7 %rename bvec and bval files too for DTI
            if exist('isDti','var') && isDti %rename bvec and bval files too for DTI
                origbvecFile = fullfile(Apth, [realname '.bvec']);
                origbvalFile = fullfile(Apth, [realname '.bval']);
                newbvecFile = fullfile(Apth, [myField '_' ID '.bvec']);
                newbvalFile = fullfile(Apth, [myField '_' ID '.bval']);
                copyfile(origbvecFile,newbvecFile);
                copyfile(origbvalFile,newbvalFile);
            end
        else
            img.(myField) = fullfile(Apth,[A ext]);
            set(btns(myIndex),'string', [myField, ': ', fullfile(Apth,[A ext])]);
        end
        %disp(Apth);
        if strcmpi(deblank(Apth),filesep)
            cd(pwd)
        else
            cd(Apth); %change directory for next search
        end
        if ~nii_check_dims({img.DTI; img.DTIrev}), warning('Fix DTI'); end;
        
    end %end filename_callback()
end %end nestedfx()

function imgs = deletePathSub (imgs)
    fields = fieldnames(imgs);
    for i = 1:numel(fields)
        if ~ischar(imgs.(fields{i})) || isempty(imgs.(fields{i})), continue; end;
        fnm = imgs.(fields{i});
        [~,n,x] = spm_fileparts(fnm);
        fnmPwd = fullfile(pwd,[n,x]);
        imgs.(fields{i}) = fnmPwd;
    end
end%deletePathSub()

function imgs = checkPathSub (imgs, matName)
    fields = fieldnames(imgs);
    for i = 1:numel(fields)
        if ~ischar(imgs.(fields{i})) || isempty(imgs.(fields{i})), continue; end;
        fnm = imgs.(fields{i});
        if exist(fnm, 'file'), continue; end;
        [~,n,x] = spm_fileparts(fnm);
        fnmPwd = fullfile(pwd,[n,x]);
        if exist(fnmPwd,'file'),
            imgs.(fields{i}) = fnmPwd;
            fprintf('Note: unable to find %s, will use %s\n', fnm, fnmPwd);
            continue;
        end;
        if strcmpi(x,'.gz') %.nii.gz
            fnmPwd = fullfile(pwd,n);
            if exist(fnmPwd,'file'),
                imgs.(fields{i}) = fnmPwd;
                fprintf('Note: unable to find %s, will use %s\n', fnm, fnmPwd);
                continue;
            end;
        end
        %unable to find files in working directory: check mat files directory
        pth = fileparts(matName);
        fnmPth = fullfile(pth,[n,x]);
        if exist(fnmPth,'file'),
            imgs.(fields{i}) = fnmPth;
            fprintf('Note: unable to find %s, will use %s\n', fnm, fnmPth);
            continue;
        end;
        if strcmpi(x,'.gz') %.nii.gz
            fnmPth = fullfile(pth,n);
            if exist(fnmPth,'file'),
                imgs.(fields{i}) = fnmPth;
                fprintf('Note: unable to find %s, will use %s\n', fnm, fnmPth);
                continue;
            end;
        end

        fprintf('WARNING: unable to find %s\n', fnm);
    end
end%checkPathSub()

function imgs = stripGzSub (imgs)
    fields = fieldnames(imgs);
    for i = 1:numel(fields)
        if ~ischar(imgs.(fields{i})) || isempty(imgs.(fields{i})), continue; end;
        fnm = imgs.(fields{i});
        [p,n,ext] = fileparts(imgs.T1);
        if strcmpi(ext,'.gz')
            fprintf('Note: removed gz from %s\n', fnm);
            imgs.(fields{i}) = fullfile(p,n);
        end
    end
end%stripGzSub()
