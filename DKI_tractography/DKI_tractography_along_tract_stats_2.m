%% calculate along tract measurements of MD/FA/MK between all pairs of a provided atlas (in diffusion space) 
function DKI_tractography_along_tract_stats(imgs,atlas,nb_nodes,scalar_maps)
global dwi_name
mask_lesion=1;
[p, n, x] = fileparts(imgs.T1);
[p_dki, n_dki , ~] = fileparts(imgs.DKI);
if exist(imgs.Lesion,'var'), [p_lesion, ~ , ~] = fileparts(imgs.Lesion);  else mask_lesion=0; end
atlas_path='/Users/Emilie/Box Sync/PhD_Projects/Recovery_Aphasia/atlas/Catani_atlas'; % make this a variable somewhere 


if exist([p '/scalars_mean_' atlas '.mat'],'file') || exist([p '/scalars_mean_' atlas '_excl.mat'],'file')
    fprintf('Skipping DKI tractography with atlas: %s\n', atlas); 
    return
end

%create gm-wm interface for seeding purposes 
if ~exist([p '/5tt_' n x],'file')
    fprintf('Running mrtrix 5ttgen and 5tt2gmwmi to %s\n', n); 
command=['5ttgen fsl ' p '/wb' n x ' ' p '/w5tt_' n x ' -force -quiet'];
system(command)
command=['5tt2gmwmi ' p '/w5tt_' n x ' ' p '/wgmwmi_' n x ' -force -quiet'];
system(command)
nfa=[p '/n' dwi_name '_fa_dki' x];
oldNormSub({[ p '/wb' n x],[p '/w5tt_' n x],[p '/wgmwmi_' n x]},nfa,8,10,0);
movefile([p '/ww5tt_' n x],[p '/5tt_' n x]);
movefile([p '/wwgmwmi_' n x],[p '/gmwmi_' n x]);
else
    fprintf('Skipping generating WM-GM interface: found %s\n', ['gmwmi_' n x]);
end

%%
% threshold gmwm interface at 50% chance 
gmwm=spm_read_vols(spm_vol([p '/gmwmi_' n x]));
gmwm=gmwm>0; % threshold was picked arbitrary, could likely be optimized 

if strcmpi(atlas,'jhu')
    atlasext = '_roi';
    hdr=spm_vol([ p_dki '/' n_dki  atlasext x]);
    ROI=spm_read_vols(hdr);  % read in atlas ROIs in native diffusion space
    %index_ROI=[1:2:93,96,97,99,113,115,167,182:2:188]; % left gray matter
    index_ROI=[7 11 15 31 35 37 39 184 186 1 9 13 25 27 29 49 69 71 41 43]; % uncomment for only language specific and domain general ROIs  
    %index_ROI=[1:max(ROI(:))];
else
    atlasext = ['_roi_' atlas];
    hdr=spm_vol([ p_dki '/' n_dki  atlasext x]);
    ROI=spm_read_vols(hdr);  % read in atlas ROIs in native diffusion space
    index_ROI=(1:max(ROI(:))); % if JHU do only language specific and domain general ROIs     
end

if mask_lesion && exist([p_lesion '/exclude' x],file)
    nfa=[p '/n' dwi_name '_fa_dki' x];
    copyfile([p '/wb' n x],[p '/Twb' n x]);
    oldNormSub({[ p '/b' n x],[p_lesion '/exclude' x]},nfa,8,10,0);
    copyfile([p '/Twb' n x],[p '/wb' n x]); % do not overwrite original wbT1
    movefile([p '/wexclude' x],[p '/DKI_exclude' x]);
    lesion = ['-exclude ' p_dki '/DKI_exclude.nii'];
    ext='_excl';
else 
    lesion='';
    ext='';
    fprintf('Warning: no voxels are excluded from tractography. If this is unexpected, make sure exclude.nii is in the same folder as imgs.Lesion')
end

%delete empty ROIs
count=1;
for i=1:length(index_ROI)
if sum(sum(sum((ROI==index_ROI(i)) & (gmwm>0))))
   temp(count)=i;
   count=count+1;
end
end
index_ROI=index_ROI(temp);
combinations=nchoosek(index_ROI,2);
clear temp count 
%%

% initialize 
if mod(nb_nodes,2)==0, nb_nodes=nb_nodes+1;end % when using tie at center you need uneven number of points


total_tracks=zeros(length(combinations),1);
total_tracks_final=zeros(max(ROI(:)),max(ROI(:)));

%track_mrtrix=cell(length(combinations),10000);

pixel_size = abs(diag(hdr.mat))';
pixel_size = pixel_size(1:3); 
transform=hdr.mat; % transformation matrix used to calculate tranformation from world coordinate space to image space 

dim=hdr.dim;

mkdir([p '/temp']);

parfor i=1:length(combinations)
fprintf('%s: ROI %d to ROI %d\t',atlas,combinations(i,1),combinations(i,2));
hdr_loop=spm_vol([ p_dki '/' n_dki atlasext x]);  

%% Create seed and include ROI 
seed=(ROI==combinations(i,1)) & (gmwm>0); hdr_loop.fname=[p '/temp/seed_' num2str(combinations(i,1)) '_' num2str(combinations(i,2)) '.nii']; spm_write_vol(hdr_loop,seed>0); % seed ROI 
include=(ROI==combinations(i,2))& (gmwm>0); hdr_loop.fname=[p '/temp/include_' num2str(combinations(i,1)) '_' num2str(combinations(i,2)) '.nii']; spm_write_vol(hdr_loop,include>0); % end ROI
%% Run MRtrix tractography - deterministic, dODF based, angular threshold 60, FA>0.1
% seed to end
%command=['tckgen -algorithm SD_STREAM -seed_image ' p '/temp/seed_' num2str(combinations(i,1)) '_' num2str(combinations(i,2)) '.nii -cutoff 0.1 -include ' p '/temp/include_' num2str(combinations(i,1)) '_' num2str(combinations(i,2)) '.nii ' lesion ' ' p '/SH_coeff.nii ' p '/temp/seed_include_' num2str(combinations(i,1)) '_' num2str(combinations(i,2)) '.tck -seeds 10000 -select 10000 -angle 60 -stop -force -seed_unidirectional -quiet'];
%system(command);
command=['tckgen -algorithm SD_STREAM -seed_image ' p '/temp/seed_' num2str(combinations(i,1)) '_' num2str(combinations(i,2)) '.nii -cutoff 0 -include ' p '/temp/include_' num2str(combinations(i,1)) '_' num2str(combinations(i,2)) '.nii ' lesion ' ' p '/SH_coeff.nii ' p '/temp/seed_include_' num2str(combinations(i,1)) '_' num2str(combinations(i,2)) '.tck -mask ' p '/mask.nii -seeds 10000 -select 10000 -angle 60 -stop -force -seed_unidirectional -quiet'];
system(command);
% end to seed 
 %command=['tckgen -algorithm SD_STREAM -seed_image ' p '/temp/include_' num2str(combinations(i,1)) '_' num2str(combinations(i,2)) '.nii -cutoff 0.1 -include ' p '/temp/seed_' num2str(combinations(i,1)) '_' num2str(combinations(i,2)) '.nii ' lesion ' ' p '/SH_coeff.nii ' p '/temp/include_seed_' num2str(combinations(i,1)) '_' num2str(combinations(i,2)) '.tck -seeds 10000 -select 10000 -angle 60 -stop -force -seed_unidirectional -quiet'];
% system(command);
command=['tckgen -algorithm SD_STREAM -seed_image ' p '/temp/include_' num2str(combinations(i,1)) '_' num2str(combinations(i,2)) '.nii -cutoff 0 -include ' p '/temp/seed_' num2str(combinations(i,1)) '_' num2str(combinations(i,2)) '.nii ' lesion ' ' p '/SH_coeff.nii ' p '/temp/include_seed_' num2str(combinations(i,1)) '_' num2str(combinations(i,2)) '.tck -mask ' p '/mask.nii -seeds 10000 -select 10000 -angle 60 -stop -force -seed_unidirectional -quiet'];
system(command);
% combine both results into one .tck 
command=['tckedit ' p '/temp/seed_include_' num2str(combinations(i,1)) '_' num2str(combinations(i,2)) '.tck ' p '/temp/include_seed_' num2str(combinations(i,1)) '_' num2str(combinations(i,2)) '.tck ' p '/temp/combined_' num2str(combinations(i,1)) '_' num2str(combinations(i,2)) '.tck -force -quiet'];
system(command);
        
%% Calculate along tract metrics (cite John Colby work: https://www.sciencedirect.com/science/article/pii/S1053811911012833?via%3Dihub and mrDiffusion) 

tracks = read_mrtrix_tracks ([ p '/temp/combined_' num2str(combinations(i,1)) '_' num2str(combinations(i,2)) '.tck']); % tracks from mrtrix are in world coordinates (voxels) 
delete([p '/temp/seed_' num2str(combinations(i,1)) '_' num2str(combinations(i,2)) '.nii'],[p '/temp/include_' num2str(combinations(i,1)) '_' num2str(combinations(i,2)) '.nii'],[p '/temp/seed_include_' num2str(combinations(i,1)) '_' num2str(combinations(i,2)) '.tck'],[p '/temp/include_seed_' num2str(combinations(i,1)) '_' num2str(combinations(i,2)) '.tck'],[p '/temp/combined_' num2str(combinations(i,1)) '_' num2str(combinations(i,2)) '.tck']);
total_tracks(i)=length(tracks.data); % calculate how many tracts were found between the two ROIs 
fprintf('Total Tracks: %d\n', total_tracks(i));
track_mrtrix(i)=tracks';
end
   
%reindex because parfor didnt let me index it properly from the beginning
for i=1:length(combinations)
    total_tracks_final(combinations(i,1),combinations(i,2))=total_tracks(i);
end

save([p '/total_tracks_' atlas ext '.mat'],'total_tracks_final');
track_mrtrix_tck=track_mrtrix(1);
track_mrtrix_tsf=track_mrtrix(1);

counter=0;
for i=1:length(combinations)
        for node=1:total_tracks(i)
        track_mrtrix_tck.data{node+counter}=track_mrtrix(i).data{node};
        track_mrtrix_tsf.data{node+counter}=i*ones(nb_nodes,1);
        end
       counter=counter+total_tracks(i);
end
track_mrtrix_tck.count=num2str(3);
track_mrtrix_tck.total_count=num2str(counter);
write_mrtrix_tracks(track_mrtrix_tck,[p '/all_' atlas ext '.tck']);

track_mrtrix_tsf.count=num2str(1);
track_mrtrix_tsf.total_count=num2str(counter);
write_mrtrix_tsf(track_mrtrix_tsf,[p '/all_' atlas ext '.tsf']);

%% segmentation 
if exist('nfa','var'), else nfa=[p '/n' dwi_name '_fa_dki' x]; end
    
oldNormSub({[ p '/wb' n x],[atlas_path '/Inferior_Occipito_Frontal_Fasciculus_Left' x],[atlas_path '/Inferior_Longitudinal_Fasciculus_Left' x],[atlas_path '/Anterior_Segment_Left' x],[atlas_path '/Posterior_Segment_Left' x],[atlas_path '/Long_Segment_Left' x],[atlas_path '/Uncinate_Left' x]},nfa,8,10,0);
movefile([atlas_path '/wInferior_Occipito_Frontal_Fasciculus_Left' x],[p '/IFOF' x]);
movefile([atlas_path '/wInferior_Longitudinal_Fasciculus_Left' x],[p '/ILF' x]);
movefile([atlas_path '/wAnterior_Segment_Left' x],[p '/Arc1' x]);
movefile([atlas_path '/wLong_Segment_Left' x],[p '/Arc2' x]);
movefile([atlas_path '/wPosterior_Segment_Left' x],[p '/Arc3' x]);
movefile([atlas_path '/wUncinate_Left' x],[p '/UNC' x]);
atlas_maps={'IFOF','ILF','Arc1','Arc2','Arc3','UNC'};

tracks_mm=track_mrtrix_tck; % initialize the matrix needed to transform tracks from voxels to mm 
total_tracks_int=length(tracks_mm.data); % calculate how many tracts were found between the two ROIs 

header=struct;
header.voxel_size=pixel_size;
header.dim=dim;
header.n_scalars=0;

        for nb_tracks=1:length(track_mrtrix_tck.data)
        int=transform\[track_mrtrix_tck.data{nb_tracks} ones(length(track_mrtrix_tck.data{nb_tracks}),1)]'; % transform to image space
        tracks_mm.data{nb_tracks}=int(1:3,:)'.*repmat(pixel_size,length(track_mrtrix_tck.data{nb_tracks}),1);  % convert to mm 
        end

        for nb_tracks=1:length(track_mrtrix_tck.data) % change from .tck format to .trk format 
        tracks_trk(nb_tracks).matrix=tracks_mm.data{nb_tracks};
        tracks_trk(nb_tracks).nPoints=length(tracks_mm.data{nb_tracks});
        end
        
%  interpolated_trk=trk_interp(tracks_trk,100,[],0); % interpolate tracts into 100 nodes, will make the number of nodes uneven because of the tie at center option 
%  tracks_interp     = trk_flip(header, interpolated_trk, [1 1 1]);
%  tracks_interp_str = trk_restruc(tracks_interp);
% 
hdr=spm_vol([p '/' atlas_maps{1} '.nii']);
map=double(spm_read_vols(hdr)>0);
m_test(:,1) = trk_add_sc_Emilie(header, tracks_interp_str, map, atlas_maps{1});
for sc_map=2:length(atlas_maps) 
        hdr=spm_vol([p '/' atlas_maps{sc_map} '.nii']);
        map=double(spm_read_vols(hdr)>0);
        m_test(:,sc_map) = trk_add_sc_Emilie(header, tracks_interp_str, map, atlas_maps{sc_map});
end 

[m,index]=max(m_test,[],2);


index_help=(1:1:total_tracks_int);

% clean IFOF 
        IFOF_index=(m~=0)&(index==1) & (m>0.5);
        tracks_IFOF=tracks_trk(IFOF_index);
        interpolated_trk=trk_interp(tracks_IFOF,101,[],1); % interpolate tracts into 100 nodes, will make the number of nodes uneven because of the tie at center option 
        tracks_interp     = trk_flip(header, interpolated_trk, [1 1 1]);
        track_mean=mean(tracks_interp,3);
        clear weights weights_param        
        for node=1:101
        int=tril(cov(permute(tracks_interp(node,:,:),[3 2 1])));
        cov_matrix=int([1:3 5:6 9]');
        sigma=[cov_matrix(1:3)'
               0 cov_matrix(4:5)'
               0 0 cov_matrix(6)'];
        sigma=sigma+sigma'-diag(diag(sigma));
        mu=track_mean(node,:);
        d=bsxfun(@minus,squeeze(tracks_interp(node,:,:))',mu);
        % mahalanobis distance of each point on each fiber from the tract core  
        weights(node,:)=sqrt(dot(d/(sigma),d,2))';
        weights_param(node,:)=mvnpdf(squeeze(tracks_interp(node,:,:))',mu,sigma)';     
        end
        
        index_keep=(sum(weights<3)==101); % remove tracks that are more than 4 standard deviations away from tract mean 
        lengths=trk_length(tracks_IFOF);
        
       int=index_help(IFOF_index);
       tracks_IFOF=track_mrtrix_tck;
       tracks_IFOF.data=track_mrtrix_tck.data(int( (index_keep) & (lengths>prctile(lengths,50)) ));
       tracks_IFOF.total_count=sum(IFOF_index);
      write_mrtrix_tracks(tracks_IFOF,[p '/IFOF.tck']);
       
       %trk
       int2=zeros(1,total_tracks_int);
       int2(int( (index_keep) & (lengths>prctile(lengths,50)) ))=1;
       IFOF_index=int2;
%ILF       
       ILF_index=(m~=0)&(index==2) & (m>0.5);
       tracks_ILF=tracks_trk(ILF_index);
        interpolated_trk=trk_interp(tracks_ILF,101,[],1); % interpolate tracts into 100 nodes, will make the number of nodes uneven because of the tie at center option 
        tracks_interp     = trk_flip(header, interpolated_trk, [1 1 1]);
        track_mean=mean(tracks_interp,3);
        clear weights weights_param        
        for node=1:101
        int=tril(cov(permute(tracks_interp(node,:,:),[3 2 1])));
        cov_matrix=int([1:3 5:6 9]');
        sigma=[cov_matrix(1:3)'
               0 cov_matrix(4:5)'
               0 0 cov_matrix(6)'];
        sigma=sigma+sigma'-diag(diag(sigma));
        mu=track_mean(node,:);
        d=bsxfun(@minus,squeeze(tracks_interp(node,:,:))',mu);
        % mahalanobis distance of each point on each fiber from the tract core  
        weights(node,:)=sqrt(dot(d/(sigma),d,2))';
        weights_param(node,:)=mvnpdf(squeeze(tracks_interp(node,:,:))',mu,sigma)';     
        end
        
        index_keep=(sum(weights<3)==101); % remove tracks that are more than 4 standard deviations away from tract mean 
        lengths=trk_length(tracks_ILF);
        
       int=index_help(ILF_index);
       tracks_ILF=track_mrtrix_tck;
       tracks_ILF.data=track_mrtrix_tck.data(int( (index_keep) & (lengths>prctile(lengths,50)) ));
       tracks_ILF.total_count=sum(ILF_index);
      write_mrtrix_tracks(tracks_ILF,[p '/ILF.tck']);
       %trk
       int2=zeros(1,total_tracks_int);
       int2(int( (index_keep) & (lengths>prctile(lengths,50)) ))=1;
       ILF_index=int2;
       
%ARC1       
       Arc1_index=(m~=0)&(index==3) & (m>0.5);
       tracks_Arc1=tracks_trk(Arc1_index);
        interpolated_trk=trk_interp(tracks_Arc1,101,[],1); % interpolate tracts into 100 nodes, will make the number of nodes uneven because of the tie at center option
        tracks_interp     = trk_flip(header, interpolated_trk, [1 1 1]);
        track_mean=mean(tracks_interp,3);
        clear weights weights_param
        for node=1:101
        int=tril(cov(permute(tracks_interp(node,:,:),[3 2 1])));
        cov_matrix=int([1:3 5:6 9]');
        sigma=[cov_matrix(1:3)'
               0 cov_matrix(4:5)'
               0 0 cov_matrix(6)'];
        sigma=sigma+sigma'-diag(diag(sigma));
        mu=track_mean(node,:);
        d=bsxfun(@minus,squeeze(tracks_interp(node,:,:))',mu);
        % mahalanobis distance of each point on each fiber from the tract core
        weights(node,:)=sqrt(dot(d/(sigma),d,2))';
        weights_param(node,:)=mvnpdf(squeeze(tracks_interp(node,:,:))',mu,sigma)';
        end

        index_keep=(sum(weights<3)==101); % remove tracks that are more than 4 standard deviations away from tract mean
        lengths=trk_length(tracks_Arc1);

       int=index_help(Arc1_index);
       tracks_Arc1=track_mrtrix_tck;
       tracks_Arc1.data=track_mrtrix_tck.data(int( (index_keep) & (lengths>prctile(lengths,50)) ));
      tracks_Arc1.total_count=sum(Arc1_index);
      write_mrtrix_tracks(tracks_Arc1,[p '/Arc1.tck']);
              %trk
       int2=zeros(1,total_tracks_int);
       int2(int( (index_keep) & (lengths>prctile(lengths,50)) ))=1;
       Arc1_index=int2;
%ARC2       
       Arc2_index=(m~=0)&(index==4) & (m>0.5);
       tracks_Arc2=tracks_trk(Arc2_index);
        interpolated_trk=trk_interp(tracks_Arc2,101,[],1); % interpolate tracts into 100 nodes, will make the number of nodes uneven because of the tie at center option 
        tracks_interp     = trk_flip(header, interpolated_trk, [1 1 1]);
        track_mean=mean(tracks_interp,3);
        clear weights weights_param        
        for node=1:101
        int=tril(cov(permute(tracks_interp(node,:,:),[3 2 1])));
        cov_matrix=int([1:3 5:6 9]');
        sigma=[cov_matrix(1:3)'
               0 cov_matrix(4:5)'
               0 0 cov_matrix(6)'];
        sigma=sigma+sigma'-diag(diag(sigma));
        mu=track_mean(node,:);
        d=bsxfun(@minus,squeeze(tracks_interp(node,:,:))',mu);
        % mahalanobis distance of each point on each fiber from the tract core  
        weights(node,:)=sqrt(dot(d/(sigma),d,2))';
        weights_param(node,:)=mvnpdf(squeeze(tracks_interp(node,:,:))',mu,sigma)';     
        end
        
        index_keep=(sum(weights<3)==101); % remove tracks that are more than 4 standard deviations away from tract mean 
        lengths=trk_length(tracks_Arc2);
        
       int=index_help(Arc2_index);
       tracks_Arc2=track_mrtrix_tck;
       tracks_Arc2.data=track_mrtrix_tck.data(int( (index_keep) & (lengths>prctile(lengths,50)) ));
      tracks_Arc2.total_count=sum(Arc2_index);
      write_mrtrix_tracks(tracks_Arc2,[p '/Arc2.tck']);
       
              %trk
       int2=zeros(1,total_tracks_int);
       int2(int( (index_keep) & (lengths>prctile(lengths,50)) ))=1;
       Arc2_index=int2;
       
%ARC3       
       Arc3_index=(m~=0)&(index==5) & (m>0.5);
       tracks_Arc3=tracks_trk(Arc3_index);
        interpolated_trk=trk_interp(tracks_Arc3,101,[],1); % interpolate tracts into 100 nodes, will make the number of nodes uneven because of the tie at center option 
        tracks_interp     = trk_flip(header, interpolated_trk, [1 1 1]);
        track_mean=mean(tracks_interp,3);
        clear weights weights_param        
        for node=1:101
        int=tril(cov(permute(tracks_interp(node,:,:),[3 2 1])));
        cov_matrix=int([1:3 5:6 9]');
        sigma=[cov_matrix(1:3)'
               0 cov_matrix(4:5)'
               0 0 cov_matrix(6)'];
        sigma=sigma+sigma'-diag(diag(sigma));
        mu=track_mean(node,:);
        d=bsxfun(@minus,squeeze(tracks_interp(node,:,:))',mu);
        % mahalanobis distance of each point on each fiber from the tract core  
        weights(node,:)=sqrt(dot(d/(sigma),d,2))';
        weights_param(node,:)=mvnpdf(squeeze(tracks_interp(node,:,:))',mu,sigma)';     
        end
        
        index_keep=(sum(weights<3)==101); % remove tracks that are more than 4 standard deviations away from tract mean 
        lengths=trk_length(tracks_Arc3);
        
       int=index_help(Arc3_index);
       tracks_Arc3=track_mrtrix_tck;
       tracks_Arc3.data=track_mrtrix_tck.data(int( (index_keep) & (lengths>prctile(lengths,50)) ));
       tracks_Arc3.total_count=sum(Arc3_index);
       write_mrtrix_tracks(tracks_Arc3,[p '/Arc3.tck']);  
              %trk
       int2=zeros(1,total_tracks_int);
       int2(int( (index_keep) & (lengths>prctile(lengths,50)) ))=1;
       Arc3_index=int2;
%UNC       
      UNC_index=(m~=0)&(index==6) & (m>0.5);
       tracks_UNC=tracks_trk(UNC_index);
        interpolated_trk=trk_interp(tracks_UNC,101,[],1); % interpolate tracts into 100 nodes, will make the number of nodes uneven because of the tie at center option 
        tracks_interp     = trk_flip(header, interpolated_trk, [1 1 1]);
        track_mean=mean(tracks_interp,3);
        clear weights weights_param        
        for node=1:101
        int=tril(cov(permute(tracks_interp(node,:,:),[3 2 1])));
        cov_matrix=int([1:3 5:6 9]');
        sigma=[cov_matrix(1:3)'
               0 cov_matrix(4:5)'
               0 0 cov_matrix(6)'];
        sigma=sigma+sigma'-diag(diag(sigma));
        mu=track_mean(node,:);
        d=bsxfun(@minus,squeeze(tracks_interp(node,:,:))',mu);
        % mahalanobis distance of each point on each fiber from the tract core  
        weights(node,:)=sqrt(dot(d/(sigma),d,2))';
        weights_param(node,:)=mvnpdf(squeeze(tracks_interp(node,:,:))',mu,sigma)';     
        end
        
        index_keep=(sum(weights<3)==101); % remove tracks that are more than 4 standard deviations away from tract mean 
        lengths=trk_length(tracks_UNC);
        
       int=index_help(UNC_index);
       tracks_UNC=track_mrtrix_tck;
       tracks_UNC.data=track_mrtrix_tck.data(int( (index_keep) & (lengths>prctile(lengths,50)) ));
      tracks_UNC.total_count=sum(UNC_index);
      write_mrtrix_tracks(tracks_UNC,[p '/UNC.tck']);         
      
        %trk
       int2=zeros(1,total_tracks_int);
       int2(int( (index_keep) & (lengths>prctile(lengths,50)) ))=1;
       UNC_index=int2;

track_density_IFOF=zeros(dim);
track_density_ILF=zeros(dim);
track_density_Arc1=zeros(dim);
track_density_Arc2=zeros(dim);
track_density_Arc3=zeros(dim);
track_density_UNC=zeros(dim);

% save trk        
for i=1:length(tracks_trk)
if IFOF_index(i)==1
tracks_trk(i).matrix(:,4)=ones(1,length(tracks_trk(i).matrix));
vox=ceil(tracks_trk(i).matrix(:,1:3)./repmat(header.voxel_size, tracks_trk(i).nPoints,1));
[~,ind]=unique(vox,'rows');
track_density_IFOF(sub2ind(dim,vox(ind,1),vox(ind,2),vox(ind,3)))=track_density_IFOF(sub2ind(dim,vox(ind,1),vox(ind,2),vox(ind,3)))+1;
elseif ILF_index(i)==1
tracks_trk(i).matrix(:,4)=2*ones(1,length(tracks_trk(i).matrix));    
vox=ceil(tracks_trk(i).matrix(:,1:3)./repmat(header.voxel_size, tracks_trk(i).nPoints,1));
[~,ind]=unique(vox,'rows');
track_density_ILF(sub2ind(dim,vox(ind,1),vox(ind,2),vox(ind,3)))=track_density_ILF(sub2ind(dim,vox(ind,1),vox(ind,2),vox(ind,3)))+1;
elseif Arc1_index(i)==1
tracks_trk(i).matrix(:,4)=3*ones(1,length(tracks_trk(i).matrix));
vox=ceil(tracks_trk(i).matrix(:,1:3)./repmat(header.voxel_size, tracks_trk(i).nPoints,1));
[~,ind]=unique(vox,'rows');
track_density_Arc1(sub2ind(dim,vox(ind,1),vox(ind,2),vox(ind,3)))=track_density_Arc1(sub2ind(dim,vox(ind,1),vox(ind,2),vox(ind,3)))+1;
elseif Arc2_index(i)==1
tracks_trk(i).matrix(:,4)=4*ones(1,length(tracks_trk(i).matrix));   
vox=ceil(tracks_trk(i).matrix(:,1:3)./repmat(header.voxel_size, tracks_trk(i).nPoints,1));
[~,ind]=unique(vox,'rows');
track_density_Arc2(sub2ind(dim,vox(ind,1),vox(ind,2),vox(ind,3)))=track_density_Arc2(sub2ind(dim,vox(ind,1),vox(ind,2),vox(ind,3)))+1;
elseif Arc3_index(i)==1
 tracks_trk(i).matrix(:,4)=5*ones(1,length(tracks_trk(i).matrix));  
 [~,ind]=unique(vox,'rows');
 track_density_Arc3(sub2ind(dim,vox(ind,1),vox(ind,2),vox(ind,3)))=track_density_Arc3(sub2ind(dim,vox(ind,1),vox(ind,2),vox(ind,3)))+1;

elseif UNC_index(i)==1
tracks_trk(i).matrix(:,4)=6*ones(1,length(tracks_trk(i).matrix));
vox=ceil(tracks_trk(i).matrix(:,1:3)./repmat(header.voxel_size, tracks_trk(i).nPoints,1));
[~,ind]=unique(vox,'rows');
track_density_UNC(sub2ind(dim,vox(ind,1),vox(ind,2),vox(ind,3)))=track_density_UNC(sub2ind(dim,vox(ind,1),vox(ind,2),vox(ind,3)))+1;
else
delete(i)=1;         
end
end
tracks_trk(delete>0)=[];

hdr.dt=[8 0];
hdr.fname=[p '/IFOF_TD.nii'];
spm_write_vol(hdr,track_density_IFOF);
hdr.fname=[p '/ILF_TD.nii'];
spm_write_vol(hdr,track_density_ILF);
hdr.fname=[p '/Arc1_TD.nii'];
spm_write_vol(hdr,track_density_Arc1);
hdr.fname=[p '/Arc2_TD.nii'];
spm_write_vol(hdr,track_density_Arc2);
hdr.fname=[p '/Arc3_TD.nii'];
spm_write_vol(hdr,track_density_Arc3);
hdr.fname=[p '/UNC_TD.nii'];
spm_write_vol(hdr,track_density_UNC);


oldNormstring=cell(1,length(atlas_maps)+1);
oldNormstring{1}=nfa;
for par=1:length(atlas_maps)
DKI_par=[p '/' atlas_maps{par} '_TD.nii'];
oldNormstring(par+1)={DKI_par};
end
oldNormSub(oldNormstring,[ p '/wb' n x], 8, 8 );

mkdir([p '/niistat_inputs/'])
for par=1:length(scalar_maps)
hdr=spm_vol([p '/w' dwi_name '_' scalar_maps{par} '_dki' x]);
vol=spm_read_vols(hdr);
 for td=1:length(atlas_maps)
 hdr=spm_vol([p '/w'  atlas_maps{td} '_TD' x]);
 track=spm_read_vols(hdr);
 combined=(track>0).*vol;
 hdr.fname=[p '/w' atlas_maps{td} '_' scalar_maps{par} x ];
 spm_write_vol(hdr,combined);
 end
end

copyfile(nfa,[p '/int.nii'])
m_coreg=coreg_DKI(nfa,[p '/wb' n '.nii'],[p '/int.nii'],p);    
transform2=spm_vol([p '/wb' n '.nii']);
transform2=transform2.mat;

header_sc.id_string='TRACK ';
header_sc.origin=[0 0 0];
header_sc.scalar_name=repmat(' ',[10,20]);
header_sc.n_properties=0;
header_sc.property_name=repmat(' ',[10,20]);
header_sc.vox_to_ras=transform2*m_coreg;
header_sc.reserved=repmat(' ',[444,1]);
header_sc.voxel_order='LPS ';
header_sc.pad2='LPS ';
header_sc.image_orientation_patient=[1 0 0 0 1 0];
header_sc.pad1='  ';
header_sc.invert_x=0;
header_sc.invert_y=0;
header_sc.invert_z=0;
header_sc.swap_xy=0;
header_sc.swap_yz=0;
header_sc.swap_zx=0;
header_sc.n_count=length(tracks_trk);
header_sc.version=2;
header_sc.hdr_size=1000; 
header_sc.n_scalars=1; 
header_sc.dim= header.dim;
header_sc.voxel_size=header.voxel_size;
trk_write(header_sc,tracks_trk,[p '/Language_tracks_MNI.trk'])
 
function oldNormSub(src, tar, smoref, reg, interp)
%coregister T2 to match T1 image, apply to lesion
if isempty(src) || isempty(tar), return; end
if ~exist('smoref','var'), smoref = 0; end
if ~exist('reg','var'), reg = 1; end
if ~exist('interp','var'), interp = 1; end
matlabbatch{1}.spm.tools.oldnorm.estwrite.subj.source = src(1);
matlabbatch{1}.spm.tools.oldnorm.estwrite.subj.wtsrc = '';
matlabbatch{1}.spm.tools.oldnorm.estwrite.subj.resample = src(:);
matlabbatch{1}.spm.tools.oldnorm.estwrite.eoptions.template = {tar};
matlabbatch{1}.spm.tools.oldnorm.estwrite.eoptions.weight = '';
matlabbatch{1}.spm.tools.oldnorm.estwrite.eoptions.smosrc = 8;
matlabbatch{1}.spm.tools.oldnorm.estwrite.eoptions.smoref = smoref; % <-- !!!
matlabbatch{1}.spm.tools.oldnorm.estwrite.eoptions.regtype = 'mni';
matlabbatch{1}.spm.tools.oldnorm.estwrite.eoptions.cutoff = 25;
matlabbatch{1}.spm.tools.oldnorm.estwrite.eoptions.nits = 16;
matlabbatch{1}.spm.tools.oldnorm.estwrite.eoptions.reg = reg;
matlabbatch{1}.spm.tools.oldnorm.estwrite.roptions.preserve = 0;
matlabbatch{1}.spm.tools.oldnorm.estwrite.roptions.bb = [nan nan nan; nan nan nan];%[-78 -112 -70; 78 76 85];
matlabbatch{1}.spm.tools.oldnorm.estwrite.roptions.vox = [nan nan nan];%[1 1 1];
matlabbatch{1}.spm.tools.oldnorm.estwrite.roptions.interp = interp;
matlabbatch{1}.spm.tools.oldnorm.estwrite.roptions.wrap = [0 0 0];
matlabbatch{1}.spm.tools.oldnorm.estwrite.roptions.prefix = 'w';
spm_jobman('run',matlabbatch);

