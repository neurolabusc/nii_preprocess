function oxford_asl_test(fnm, fsldir)
%text if oxford_asl bugs patched
%https://github.com/ibme-qubic/oxford_asl/issues/8
% fnm : (optional) location of oxford_asl
% fsldir: (optional) location of main FSL installation
%Examples
% oxford_asl_test %default install
% oxford_asl_test('/home/chris/oxford_asl/oxford_asl')
% oxford_asl_test('~/oxford_asl/oxford_asl', '/usr/local/fsl')


if ~exist('fnm', 'var') || ~exist(fnm, 'file') || exist(fnm, 'dir')
    pth = getenv('FSLDIR');
    if isempty(pth)
       pth = '/usr/local/fsl';
    end
    fnm = fullfile(fullfile(pth,'bin'), 'oxford_asl');
    if ~exist(fnm,'file')
        error('%s unable to find %s\n', mfilename, fnm); 
    end
    fprintf('Testing %s\n', fnm);
end

allOK = true;
txt = fileread(fnm);
bad = 'asl2struc.mat';
good = 'asl2struct.mat';
if contains(txt, bad)
    warning('%s: change **%s** to read **%s**\n"', fnm, bad, good);
    allOK = false;
end
bad = 'if [ -f $tempdir/struc_bet.* -a -z $mask ];';
good = 'if [ -f $tempdir/struc_bet.* -a -z "$mask" ];';
if contains(txt, bad)
    warning('%s: change **%s** to read **%s**\n"', fnm, bad, good);
    allOK = false;
end

if ~exist('fsldir', 'var') || ~exist(fsldir, 'file')
    fsldir = getenv('FSLDIR');
end
if isempty(fsldir)
    fsldir = '/usr/local/fsl';
end
if exist(fsldir, 'dir')
    %by default, Ubuntu 20.04 does not define "python" (python3)
    fnms = {'/bin/fsl_abspath', '/python/asl/reg.py', '/python/asl/fslwrap.py', '/python/asl/fslhelpers.py'};
    for i = 1 : numel(fnms)
        pynm = fullfile(fsldir, fnms{i});
        if ~exist(pynm, 'file')
            warning('unable to find **%s**\n"', pynm);
            allOK = false;
            continue;
        end
        txt = fileread(pynm);
        if startsWith(txt, '#!/bin/env python') && ~startsWith(txt, '#!/bin/env python3')
            warning('%s: change **python** to **fslpythonw**\n"', pynm);
            allOK = false;
            continue;            
        end
    end
else
	warning('%s: change **python** to **fslpythonw**\n"', pynm);
	allOK = false;    
end

if allOK
    fprintf('%s success!\n', mfilename);
else
   error('%s failed: see warnings\n', mfilename); 
end