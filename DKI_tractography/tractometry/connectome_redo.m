%% connectome redo 

basedir='/Volumes/Data2/Aphasia_Project/'
folders=dir(basedir);
index_ROI=[1:2:93,96,97,99,113,115,167,182:2:188]; % left gray matter
combinations=nchoosek(index_ROI,2);
scalar_maps={'fa','md','dax','drad','mk','kax','krad','kfa'}; % calculate along tract stats of these metrics 
mkdir([basedir '/connectomes_redo'])
for i=1:length(folders) 

    if exist([basedir '/' folders(i).name '/output_redo/BL/along_tract_metrics'],'dir')
    subdir=[basedir '/' folders(i).name '/output_redo/BL/'];
    load([subdir 'along_tract_metrics/Scalar_connectome.mat']);
    subconnectome=squeeze(median(connectome(:,:,5,:),4));
    save(['/Users/Emilie/Box Sync/PhD_Projects/Recovery_Aphasia/connectomes_redo/BL_' folders(i).name 'median_MK.mat'],'subconnectome')
    end
 
    


end

basedir='/Volumes/Data2/Aphasia_Project/'
folders=dir(basedir);
index_ROI=[1:2:93,96,97,99,113,115,167,182:2:188]; % left gray matter
combinations=nchoosek(index_ROI,2);
scalar_maps={'fa','md','dax','drad','mk','kax','krad','kfa'}; % calculate along tract stats of these metrics 
%mkdir([basedir '/connectomes_redo'])
for i=1:length(folders) 

    if exist([basedir '/' folders(i).name '/output_redo/output_redo/1W/along_tract_metrics'],'dir')
    subdir=[basedir '/' folders(i).name '/output_redo/output_redo/1W/'];
    load([subdir 'along_tract_metrics/Scalar_connectome.mat']);
    subconnectome=squeeze(median(connectome(:,:,5,:),4));
    save(['/Users/Emilie/Box Sync/PhD_Projects/Recovery_Aphasia/connectomes_redo/1W_' folders(i).name 'median_MK.mat'],'subconnectome')
    end
 
    


end

%% 
subj_id={'119','123','125','134','137','138','141','148','150','171','175','177','180','187','188','189','191','192','193','195','197','199','201','205','208','209','211','214','215','217','218','220','223','224'};

basedir='/Volumes/Data2/Aphasia_Project/'
folders=dir(basedir);
index_ROI=[1:2:93,96,97,99,113,115,167,182:2:188]; % left gray matter
combinations=nchoosek(index_ROI,2);
scalar_maps={'fa','md','dax','drad','mk','kax','krad','kfa'}; % calculate along tract stats of these metrics 

scalar=5;

for i=1:length(subj_id) 

    if exist([basedir '/' subj_id{i} '/output_redo/output_redo/BL/along_tract_metrics/Scalar_connectome.mat'],'file')
    subdir=[basedir '/'  subj_id{i} '/output_redo/output_redo/BL/'];
    load([subdir 'along_tract_metrics/Scalar_connectome.mat']);
    subconnectome(:,:,:,i)=squeeze(connectome(:,:,scalar,:));
    %save(['/Users/Emilie/Box Sync/PhD_Projects/Recovery_Aphasia/connectomes_redo/BL_' folders(i).name 'median_MK.mat'],'subconnectome')
    end

end

subconnectome(:,:,:,29)=zeros(189,189,100); % problem with P215
subconnectome(subconnectome==0)=NaN;



test=squeeze(nanmedian(subconnectome(combinations(928,1),combinations(928,2),:,:),4));
test_sd=squeeze(nanstd(subconnectome(combinations(928,1),combinations(928,2),:,:),'',4));

test=squeeze(nanmedian(subconnectome(combinations(928,1),combinations(928,2),1:25,:),3));
[res(1,1),res(1,2)]=corr(semantics(~isnan(test),3),test(~isnan(test)));

test=squeeze(nanmedian(subconnectome(combinations(928,1),combinations(928,2),26:50,:),3));
[res(2,1),res(2,2)]=corr(semantics(~isnan(test),3),test(~isnan(test)));

test=squeeze(nanmedian(subconnectome(combinations(928,1),combinations(928,2),51:75,:),3));
[res(3,1),res(3,2)]=corr(semantics(~isnan(test),3),test(~isnan(test)));

test=squeeze(nanmedian(subconnectome(combinations(928,1),combinations(928,2),76:100,:),3));
[res(4,1),res(4,2)]=corr(semantics(~isnan(test),3),test(~isnan(test)));


figure, hold on
plot(test(~isnan(test)), 'k')
plot(test(~isnan(test))+test_sd(~isnan(test)), 'k--')
plot(test(~isnan(test))-test_sd(~isnan(test)), 'k--')

%Make the plot prettier
hold off, box off
xlim([0 100]), ylim([0 3])
title('\bf{Mean Diffusivity along track}')
xlabel('Distance along track (%)')
ylabel('Mean Diffusivity (MD)')

subconnectome_FU=zeros(189,189,100,length(subj_id));
for i=1:length(subj_id) 

    if exist([basedir '/' subj_id{i} '/output_redo/output_redo/1W/along_tract_metrics/Scalar_connectome.mat'],'file')
    subdir=[basedir '/'  subj_id{i} '/output_redo/output_redo/1W/'];
    load([subdir 'along_tract_metrics/Scalar_connectome.mat']);
    subconnectome_FU(:,:,:,i)=squeeze(connectome(:,:,scalar,:));
    %save(['/Users/Emilie/Box Sync/PhD_Projects/Recovery_Aphasia/connectomes_redo/BL_' folders(i).name 'median_MK.mat'],'subconnectome')
    end

end

subconnectome(isnan(subconnectome))=0;

int=(subconnectome_FU>0)& (subconnectome>0);
change=(subconnectome_FU.*int)-(subconnectome.*int);
change(change==0)=NaN;


test=squeeze(nanmedian(change(combinations(928,1),combinations(928,2),1:25,:),3));
[res(1,1),res(1,2)]=corr(semantics(~isnan(test),7),test(~isnan(test)));

test=squeeze(nanmedian(change(combinations(928,1),combinations(928,2),26:50,:),3));
[res(2,1),res(2,2)]=corr(semantics(~isnan(test),7),test(~isnan(test)));

test=squeeze(nanmedian(change(combinations(928,1),combinations(928,2),51:75,:),3));
[res(3,1),res(3,2)]=corr(semantics(~isnan(test),7),test(~isnan(test)));

test=squeeze(nanmedian(change(combinations(928,1),combinations(928,2),76:100,:),3));
[res(4,1),res(4,2)]=corr(semantics(~isnan(test),7),test(~isnan(test)));


