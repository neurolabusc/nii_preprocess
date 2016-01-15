function imgs = nii_lime_gui(loadprev)
%loadprev is true or false to load saved img name file
% 
%Examples
% nii_lime_gui %use gui
% nii_lime_gui(1) %use gui to select previous mat file
% nii_lime_gui('T1_M2094_limegui.mat') %specify previous mat file
if nargin < 1
    loadprev = [];
end
if ~isempty(loadprev)
    if ischar(loadprev) && exist(loadprev, 'file')
        [p, n, x] = fileparts(loadprev);
        f = [n, x];
    else
        [f, p] = uigetfile('*limegui.mat', 'Choose name mat file');
    end
    imgs = load(fullfile(p,f));
    imgs = checkPathSub (imgs);
    imgs = stripGzSub (imgs);
    nii_lime(imgs);
    return;
end
%Graphical interface for running nii_lime processing
imgs = [];
go = false;
Xpressed = false;
while ~go && ~Xpressed
    [img, go, Xpressed] = getfiles_x;
    if isempty(img.T1)
        fprintf('Error: T1 scan required\n');
    else
        imgs = [imgs, img]; %#ok<AGROW>
    end;
end
if isempty(imgs), return; end;
if Xpressed, return; end;
for i = 1: numel(imgs)
    img = imgs(i);
    %save filenames that will be processed - 
    % to rerun in future: img = load('name_limegui.mat'); nii_lime(img);
    [p,n] = fileparts(img.T1);
    matName = fullfile(p, [n, '_limegui.mat']);
    save(matName,'-struct', 'img');
    %process the data
    nii_lime(img);
end
%end %nii_lime_gui()

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
%stripGzSub()

function imgs = checkPathSub (imgs)
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
    fprintf('WARNING: unable to find %s\n', fnm);
end
%checkPathSub

function [img, go, Xpressed] = getfiles_x
go = true;
Xpressed = false;
alreadyWarned = false;
btns = zeros(10,1);
img.T1 = [];
img.T2 = [];
img.ASL = [];
img.fMRI = [];
img.Rest = [];
img.DTI = [];
img.DTIrev = [];
img.Lesion = [];
renameCheck = 0;
ed = 0;
d = nestedfx;
uiwait(d);

    function d = nestedfx
      d = dialog('Position',[300 300 620 360],'Name','LIME image selection','CloseRequestFcn',@my_closereq);
      topMostPos = 330;
      pos = 10:30:topMostPos;
      ed = uicontrol('Parent',d,...
                   'Style','Edit',...
                   'Position',[10 pos(end) 600 25],...
                   'String','Participant ID');
      txt = uicontrol('Parent',d,...
                   'Style','text',...
                   'Position',[10 pos(end-1) 600 25],...
                   'String','Select the desired scans');
      btns(1) = uicontrol('Parent',d,...
                   'Position',[10 pos(end-2) 600 25],...
                   'String','T1',...
                   'Callback', {@filename_callback,'T1', 1});
      btns(2) = uicontrol('Parent',d,...
                   'Position',[10 pos(end-3) 600 25],...
                   'String','T2',...
                   'Callback',{@filename_callback,'T2', 2});
      btns(3) = uicontrol('Parent',d,...
                   'Position',[10 pos(end-4) 600 25],...
                   'String','ASL',...
                   'Callback',{@filename_callback,'ASL', 3});
      btns(4) = uicontrol('Parent',d,...
                   'Position',[10 pos(end-5) 600 25],...
                   'String','fMRI',...
                   'Callback',{@filename_callback,'fMRI', 4});
      btns(5) = uicontrol('Parent',d,...
                   'Position',[10 pos(end-6) 600 25],...
                   'String','Resting State',...
                   'Callback',{@filename_callback,'Rest', 5});
      btns(6) = uicontrol('Parent',d,...
                   'Position',[10 pos(end-7) 600 25],...
                   'String','DTI',...
                   'Callback',{@filename_callback,'DTI', 6, true});
      btns(7) = uicontrol('Parent',d,...
                   'Position',[10 pos(end-8) 600 25],...
                   'String','DTI (reversed phase)',...
                   'Callback',{@filename_callback,'DTIrev', 7, true});
      btns(8) = uicontrol('Parent',d,...
                   'Position',[10 pos(end-9) 600 25],...
                   'String','Lesion',...
                   'Callback',{@filename_callback,'Lesion', 8});
               
        btnGo = uicontrol('Parent',d,...
                   'Position',[300 10 160 25],...
                   'String','Go',...
                   'Callback','delete(gcf)');
        btnNext = uicontrol('Parent',d,...
                   'Position',[10 10 160 25],...
                   'String','Next subject',...
                   'Callback',@next_callback);
        renameCheck = uicontrol('Parent',d,...
                   'style','checkbox',...
                   'units','pixels',...
                   'position',[500 10 100 15],...
                   'string','Rename files');
        
    %end nestedfx()
    function next_callback(~, ~, ~) %filename_callback(obj,callbackdata,myField)
        go = false;
        delete(gcf);
    %end %end next_callback()
    function my_closereq(~,~,~)
        go = false;
        delete(gcf);
        Xpressed = true;
    %end my_closereq()
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
            ID = get(ed,'String');
            if strcmpi(ID,'Participant ID') %user did not provide ID so don't use the field value
                ID = '';
            end
            Ashort = [myField '_' ID];
            copyfile(fullfile(Apth,[A ext]),fullfile(Apth,[Ashort newext]));
            img.(myField) = fullfile(Apth,[Ashort newext]);
            set(btns(myIndex),'string', [myField, ': ',fullfile(Apth,[Ashort newext])]);
            if myIndex == 6 || myIndex == 7 %rename bvec and bval files too for DTI
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

        cd(Apth); %change directory for next search
    %end %end filename_callback()
%end %end nestedfx()