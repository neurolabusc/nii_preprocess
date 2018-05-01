function roi2dwi

rootDir = '/Volumes/FAT500/OTA/OtaTask/ota_img/OK';
outDir = '/Volumes/FAT500/OTA/OtaTask/ota_img/Process';

[filenames, files] =  dirSub(rootDir);
subDirs = filenames([files.isdir]);
for s = 1:length(subDirs)
  subDir = subDirs{s};
  if subDir(1)=='.', continue; end;
  %subDir = fullfile(rootDir, subDir);
  voi = dirSub(fullfile(rootDir, subDir, '*.voi'));
  if isempty(voi), fprintf('Unable to find VOI for %s\n', subDir); continue; end;
  if numel(voi) > 2, error('Unable to choose: Multiple VOIs %s\n', subDir); end;
  voi = voi{1};
  isGz = false;
  nii = dirSub(fullfile(rootDir, subDir, '*.nii'));
  if isempty(nii), nii = dirSub(fullfile(rootDir, subDir, '*.nii.gz')); isGz = true; end;
  if isempty(nii), fprintf('Unable to find image for %s\n', voi); continue; end;
  if numel(nii) > 2, error('Unable to choose: Multiple NIIs %s\n', subDir); end;
  nii = nii{1};
  
  % process subdirectory
  fprintf('%s + %s\n',voi, nii);   % for example
  copyfile(fullfile(rootDir, subDir, voi), fullfile(outDir, [subDir,'.voi']) );
  if isGz
    outname = fullfile(outDir, [subDir,'_dwi.nii.gz']);
    copyfile(fullfile(rootDir, subDir, nii), outname );
    gunzip(outname);
    delete(outname);
    
  else
    copyfile(fullfile(rootDir, subDir, nii), fullfile(outDir, [subDir,'_dwi.nii']) );
  end;
end
%end function

function [filenames, files] =  dirSub(pth)
%return visible files in path
files = dir(pth);   % assume starting from current directory
files = files(arrayfun(@(x) x.name(1), files) ~= '.'); %remove hidden files
filenames = {files.name};
%dirSub()