%% calculate along tract measurements of MD/FA/MK between all pairs of a provided atlas (in diffusion space) 
function DKI_tractography_along_tract_stats(imgs,atlas,nb_nodes,scalar_maps,atlas_maps)
global dwi_name
if nargin<5
atlas_maps={'Inferior_Occipito_Frontal_Fasciculus_Left','Inferior_Longitudinal_Fasciculus_Left','Long_Segment_Left','Anterior_Segment_Left','Posterior_Segment_Left','Uncinate_Left'};
end

mask_lesion=1;

[p, n, x] = fileparts(imgs.T1);
[p_dki, n_dki , ~] = fileparts(imgs.DKI);
if exist(imgs.Lesion), [p_lesion, ~ , ~] = fileparts(imgs.Lesion);  else mask_lesion=0; end
atlas_path = fullfile(fileparts(which('nii_preprocess')), 'catani_atlas');
% 
% if exist([p '/scalars_mean_' atlas '.mat'],'file') || exist([p '/scalars_mean_' atlas '_excl.mat'],'file')
%     fprintf('Skipping DKI tractography with atlas: %s\n', atlas); 
%     return
% end

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
% threshold gmwm interface at 0% chance 
gmwm=spm_read_vols(spm_vol([p '/gmwmi_' n x]));
gmwm=gmwm>0; % threshold was picked arbitrary, could likely be optimized 

if strcmpi(atlas,'jhu')
    atlasext = '_roi';
    hdr=spm_vol([ p_dki '/' n_dki  atlasext x]);
    ROI=spm_read_vols(hdr);  % read in atlas ROIs in native diffusion space
    index_ROI=[1:2:93,96,97,99,113,115,167,182:2:188]; % left gray matter
    %index_ROI=[7 11 15 31 35 37 39 184 186 1 9 13 25 27 29 49 69 71 41 43]; % uncomment for only language specific and domain general ROIs  
    %index_ROI=[1:max(ROI(:))];
else
    atlasext = ['_roi_' atlas];
    hdr=spm_vol([ p_dki '/' n_dki  atlasext x]);
    ROI=spm_read_vols(hdr);  % read in atlas ROIs in native diffusion space
    index_ROI=(1:max(ROI(:))); % if JHU do only language specific and domain general ROIs     
end
%made exclude mask in T1 space , warp exclude.nii into native diffusion
%space. 
if mask_lesion && exist([p_lesion '/exclude' x],'file')
    nfa=[p '/n' dwi_name '_fa_dki' x];
    copyfile([p '/wb' n x],[p '/Twb' n x]);
    oldNormSub({[ p '/b' n x],[p_lesion '/exclude' x]},nfa,8,10,0);
    copyfile([p '/Twb' n x],[p '/wb' n x]); % step included to not overwrite original wbT1
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

total_tracks=zeros(length(combinations),1);
total_tracks_final=zeros(max(ROI(:)),max(ROI(:)));

pixel_size = abs(diag(hdr.mat))';
pixel_size = pixel_size(1:3); 
transform=hdr.mat; % transformation matrix used to calculate tranformation from world coordinate space to image space 
dim=hdr.dim;

mkdir([p '/temp']);

if ~exist([p '/all_' atlas ext '.tck'],'file') % skip making of tck if already exists 
parfor i=1:length(combinations)
fprintf('%s: ROI %d to ROI %d\t',atlas,combinations(i,1),combinations(i,2));
hdr_loop=spm_vol([ p_dki '/' n_dki atlasext x]);  

%% Create seed and include ROI 
seed=(ROI==combinations(i,1)) & (gmwm>0); hdr_loop.fname=[p '/temp/seed_' num2str(combinations(i,1)) '_' num2str(combinations(i,2)) '.nii']; spm_write_vol(hdr_loop,seed>0); % seed ROI 
include=(ROI==combinations(i,2))& (gmwm>0); hdr_loop.fname=[p '/temp/include_' num2str(combinations(i,1)) '_' num2str(combinations(i,2)) '.nii']; spm_write_vol(hdr_loop,include>0); % end ROI
%% Run MRtrix tractography - deterministic, dODF based, angular threshold 60, FA>0.1
% seed to end
command=['tckgen -algorithm SD_STREAM -seed_image ' p '/temp/seed_' num2str(combinations(i,1)) '_' num2str(combinations(i,2)) '.nii -cutoff 0.1 -include ' p '/temp/include_' num2str(combinations(i,1)) '_' num2str(combinations(i,2)) '.nii ' lesion ' ' p '/SH_coeff.nii ' p '/temp/seed_include_' num2str(combinations(i,1)) '_' num2str(combinations(i,2)) '.tck -seeds 10000 -select 10000 -angle 60 -stop -force -seed_unidirectional -quiet'];
system(command);
% end to seed 
command=['tckgen -algorithm SD_STREAM -seed_image ' p '/temp/include_' num2str(combinations(i,1)) '_' num2str(combinations(i,2)) '.nii -cutoff 0.1 -include ' p '/temp/seed_' num2str(combinations(i,1)) '_' num2str(combinations(i,2)) '.nii ' lesion ' ' p '/SH_coeff.nii ' p '/temp/include_seed_' num2str(combinations(i,1)) '_' num2str(combinations(i,2)) '.tck -seeds 10000 -select 10000 -angle 60 -stop -force -seed_unidirectional -quiet'];
system(command);
% combine both results into one .tck 
command=['tckedit ' p '/temp/seed_include_' num2str(combinations(i,1)) '_' num2str(combinations(i,2)) '.tck ' p '/temp/include_seed_' num2str(combinations(i,1)) '_' num2str(combinations(i,2)) '.tck ' p '/temp/combined_' num2str(combinations(i,1)) '_' num2str(combinations(i,2)) '.tck -force -quiet'];
system(command);
        
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
        track_mrtrix_tsf.data{node+counter}=i*ones(size(track_mrtrix(i).data{node},1),1);
        end
       counter=counter+total_tracks(i);
end
track_mrtrix_tck.count=num2str(3);
track_mrtrix_tck.total_count=num2str(counter);
write_mrtrix_tracks(track_mrtrix_tck,[p '/all_' atlas ext '.tck']);

track_mrtrix_tsf.count=num2str(1);
track_mrtrix_tsf.total_count=num2str(counter);
write_mrtrix_tsf(track_mrtrix_tsf,[p '/all_' atlas ext '.tsf']);
else
       fprintf('Skipping tractography: found %s\n', ['all_' atlas ext '.tck']); % dont do tractography if .tck is already created, but load in outputs for the rest of the script
       load([p '/total_tracks_' atlas ext '.mat']);
       track_mrtrix_tck=read_mrtrix_tracks([p '/all_' atlas ext '.tck']);
       total_tracks=zeros(1,length(combinations)+1); 

         % offset total_tracks with 0 
        for cnt=1:length(combinations)
            total_tracks(cnt+1)=total_tracks_final(combinations(cnt,1),combinations(cnt,2));
        end

end
%% segmentation in catani atlas ROIs 
if ~exist([p '/Language_tracks_MNI.trk'],'file')
    
if exist('nfa','var'), else nfa=[p '/n' dwi_name '_fa_dki' x]; end

% transform atlas from mni space to native diffusion space 
oldNormstring=cell(1,length(atlas_maps)+1);
oldNormstring{1}=[ p '/wb' n x];
for par=1:length(atlas_maps)
atlas_par=[atlas_path '/' atlas_maps{par} x];
oldNormstring(par+1)={atlas_par};
end
oldNormSub(oldNormstring,nfa, 8, 10,0 );

for par=1:length(atlas_maps)
movefile([atlas_path '/w' atlas_maps{par} x],[p '/' atlas_maps{par} x]);
end
if ~exist('track_mrtrix_tck','var'), track_mrtrix_tck=read_mrtrix_tracks([p '/all_' atlas ext '.tck']);, end
tracks_mm=track_mrtrix_tck; % initialize the matrix needed to transform tracks from voxels to mm 
total_tracks_int=length(tracks_mm.data); % calculate how many tracts were found between the two ROIs 

header=struct;
header.voxel_size=pixel_size;
header.dim=dim;
header.n_scalars=0;

tracks_trk=repmat(struct,1,length(track_mrtrix_tck.data)); %initialize

        for nb_tracks=1:length(track_mrtrix_tck.data)
        int=transform\[track_mrtrix_tck.data{nb_tracks} ones(length(track_mrtrix_tck.data{nb_tracks}),1)]'; % transform tracks to image space
        tracks_mm.data{nb_tracks}=int(1:3,:)'.*repmat(pixel_size,length(track_mrtrix_tck.data{nb_tracks}),1);  % convert to mm 
        end

        for nb_tracks=1:length(track_mrtrix_tck.data) % change from .tck format to .trk format 
        tracks_trk(nb_tracks).matrix=tracks_mm.data{nb_tracks};
        tracks_trk(nb_tracks).nPoints=length(tracks_mm.data{nb_tracks});
        end
        
% Calculate overlap between atlas and tracks: sample the atlas values at the coordinates of the tracks.         
for sc_map=1:length(atlas_maps) 
        hdr=spm_vol([p '/' atlas_maps{sc_map} '.nii']);
        map=double(spm_read_vols(hdr)>0);
        m_int(:,sc_map) = trk_atlas(header, tracks_trk, map, atlas_maps{sc_map});
end 

[m,index]=max(m_int,[],2); % count with which atlas the .trk overlaped the most. 


index_help=(1:1:total_tracks_int);

% clean atlas bundles 
for atl=1:length(atlas_maps)
        atlas_index=(m~=0)&(index==atl) & (m>0.5); % voting system, more than 50% needs to be part of the atlas
        if sum(atlas_index)> 10
        tracks_atlas=tracks_trk(atlas_index);
        interpolated_trk=trk_interp(tracks_atlas,101,[],1); % interpolate tracts into 100 nodes, will make the number of nodes uneven because of the tie at center option 
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
        
        index_keep=(sum(weights<3)==101); % remove tracks that are more than 3 standard deviations away from tract mean 
        lengths=trk_length(tracks_atlas);%calculate the lengths of all tracks 
        
       int=index_help(atlas_index);
       tracks_atlas=track_mrtrix_tck;
       tracks_atlas.data=track_mrtrix_tck.data(int( (index_keep) & (lengths>prctile(lengths,50)) ));
       tracks_atlas.total_count=sum(int( (index_keep) & (lengths>prctile(lengths,50)) ));
       write_mrtrix_tracks(tracks_atlas,[p '/' atlas_maps{atl} '.tck']);
       
       % save indices for a .trk
       int2=zeros(1,total_tracks_int);
       int2(int((index_keep) & (lengths>prctile(lengths,50))))=1;
       atl_index_final(:,atl)=int2;
        elseif sum(atlas_index)~=0
       tracks_atlas=track_mrtrix_tck;
       tracks_atlas.data=track_mrtrix_tck.data(atlas_index);
       tracks_atlas.total_count=sum(atlas_index);
       write_mrtrix_tracks(tracks_atlas,[p '/' atlas_maps{atl} '.tck']);
       
       % save indices for a .trk
       atl_index_final(:,atl)=atlas_index; 
        end
        
end


track_density=zeros(dim(1),dim(2),dim(3),length(atlas_maps));

    for i=1:length(tracks_trk)
        if find(atl_index_final(i,:))
    tracks_trk(i).matrix(:,4)=find(atl_index_final(i,:))*ones(1,length(tracks_trk(i).matrix));
    vox=ceil(tracks_trk(i).matrix(:,1:3)./repmat(header.voxel_size, tracks_trk(i).nPoints,1));
    [~,ind]=unique(vox,'rows');
    track_density(sub2ind([dim(1),dim(2),dim(3),length(atlas_maps)],vox(ind,1),vox(ind,2),vox(ind,3),repmat(find(atl_index_final(i,:)),[length(ind),1])))=track_density(sub2ind([dim(1),dim(2),dim(3),length(atlas_maps)],vox(ind,1),vox(ind,2),vox(ind,3),repmat(find(atl_index_final(i,:)),[length(ind),1])))+1;
        else 
            delete_index(i)=1;
        end
    end
 
tracks_trk(delete_index>0)=[];
% save track density maps 
for atl=1:length(atlas_maps)
hdr.dt=[8 0];
hdr.fname=[p '/' atlas_maps{atl} '_TD.nii'];
spm_write_vol(hdr,track_density(:,:,:,atl));
end

% convert Trck Density maps to MNI space 
oldNormstring=cell(1,length(atlas_maps)+1);
oldNormstring{1}=nfa;
for par=1:length(atlas_maps)
DKI_par=[p '/' atlas_maps{par} '_TD.nii'];
oldNormstring(par+1)={DKI_par};
end
oldNormSub(oldNormstring,[ p '/wb' n x], 8, 8 );

% save language tracks in MNI space for surfice 
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
else
         fprintf('Skipping segmentation in MNI tracts: found %s\n',[p '/Language_tracks_MNI.trk']);
   
end

%% tractometry 

mkdir([p '/along_tract_metrics'])

total_tracks_new=zeros(length(combinations),1);
for j=1:length(combinations) % all combinations of ROIs; was originally paralelized but kept getting an error

fprintf('Calculating 4D connectomes: ROI %d to ROI %d\n',combinations(j,1),combinations(j,2));

%initialize    
begin_ROI1=[];
begin_ROI2=[];
end_ROI1=[];
end_ROI2=[];

tracks_trk=struct('matrix',[],'nPoints',[]);

    
ROI1=(ROI==combinations(j,1))& (gmwm>0);
ROI2=(ROI==combinations(j,2))& (gmwm>0); 

% read in correct tracts and convert to trk format, I did this because I am
% using John Colby .trk code to interpolate. We could probably change this
% part 
        for node=1:total_tracks(j+1)
        int=hdr.mat\[track_mrtrix_tck.data{node+sum(total_tracks(1:j))} ones(length(track_mrtrix_tck.data{node+sum(total_tracks(1:j))}),1)]'; % tracks in image space
        tracks_trk(node).matrix=int(1:3,:)'.*repmat(pixel_size,length(track_mrtrix_tck.data{node+sum(total_tracks(1:j))}),1); %  image space in mm 
        tracks_trk(node).nPoints=length(track_mrtrix_tck.data{node+sum(total_tracks(1:j))});
        end
        
   if ~any(structfun(@isempty, tracks_trk(1)) )
       % find indices for tracts that start in ROI1 and ROI2
                for iTrk=1:length(tracks_trk)              
                % Translate continuous vertex coordinates into discrete voxel coordinates
                vox = ceil(tracks_trk(iTrk).matrix(:,1:3) ./ repmat(pixel_size, tracks_trk(iTrk).nPoints,1));
              
                % Index into volume to extract scalar values
                inds                = sub2ind(hdr.dim, vox(:,1), vox(:,2), vox(:,3));
                scalars_int         = ROI1(inds);
                begin_ROI1(iTrk)    = scalars_int(1); % does it start in ROI1
                end_ROI1(iTrk)    = scalars_int(end); % does it end in ROI1
                scalars_int         = ROI2(inds);
                begin_ROI2(iTrk)    = scalars_int(1); % does it start in ROI2
                end_ROI2(iTrk)    = scalars_int(end); % does it end in ROI2
                end
                all=[begin_ROI1;end_ROI1;begin_ROI2;end_ROI2]';
                index= (sum(all==[1,0,0,1],2)==4)| (sum(all==[0,1,1,0],2)==4); % note that not all tracts will receive an index. This is because mrtrix interpolates ROI1 and ROI2 when tracking 
        
    %% interpolate tracks in nb_nodes 
    if sum(index)~=0
    interpolated_trk=trk_interp(tracks_trk(index),nb_nodes,[]); % interpolate tracts into 100 nodes, will make the number of nodes uneven because of the tie at center option 
    % flip tracts that start in ROI1 and end in ROI2 this might be more
    % complicated than it needs to be 
    interpolated_trk_flip=interpolated_trk;
    all_index=1:length(tracks_trk);
    all_index_subset=all_index(index); 
    int=all_index((sum(all==[1,0,0,1],2)==4)>0); 
    [~,B]=ismember(int,all_index_subset);
         
    interpolated_trk_flip(:,:,B)=interpolated_trk(fliplr(1:end),:,B); % flip tracks that start in ROI1 , as such when we interpolate all tracts start at the same location 
   
    % write tracts
    track_mrtrix_tsf=track_mrtrix_tck;
    track_mrtrix_tck2=track_mrtrix_tck;
    track_mrtrix_tck2.data=[];
    track_mrtrix_tsf.data=[];

     for nb=1:sum(index)
     int=(transform*[(interpolated_trk_flip(:,:,nb)./repmat(pixel_size,size(interpolated_trk_flip,1),1)) ones(size(interpolated_trk_flip,1),1)]'); % convert to .tck format (world coordinates in voxels) 
     track_mrtrix_tck2.data{nb}=int(1:3,:)';
     track_mrtrix_tsf.data{nb}=(1:size(interpolated_trk_flip,1))';
     end

% uncomment if you want to output all tracks seperately 
     
%     track_mrtrix_tck2.count=1;
%     track_mrtrix_tck2.total_count=sum(index);
%     write_mrtrix_tracks(track_mrtrix_tck2,[p '/along_tract_metrics/ROI' num2str(combinations(j,1)) '_ROI' num2str(combinations(j,2)) '.tck']);
%     track_mrtrix_tsf.count=1;
%     track_mrtrix_tsf.total_count=sum(index);
%     write_mrtrix_tsf(track_mrtrix_tsf,[p '/along_tract_metrics/ROI' num2str(combinations(j,1)) '_ROI' num2str(combinations(j,2)) '.tsf']);
%    
   track_mrtrix_all(j)=track_mrtrix_tck2;
   total_tracks_new(j)=length(track_mrtrix_tck2.data); % calculate how many tracts were found between the two ROIs  

    else
        
    end
    
    else
                    % might still put something here     
    end

    
end 

% write out tracts 

total_tracks_final_new=zeros(max(ROI(:)),max(ROI(:)));
  
%reindex because parfor didnt let me index it properly from the beginning
for i=1:length(combinations)
    total_tracks_final_new(combinations(i,1),combinations(i,2))=total_tracks_new(i);
end
save([p '/along_tract_metrics/total_tracks_cleaned.mat'],'total_tracks_final_new');


track_mrtrix_tck2=track_mrtrix_all(1);
track_mrtrix_tsf=track_mrtrix_all(1);
track_mrtrix_tsf2=track_mrtrix_all(1);
counter=0;
for i=1:length(combinations)
        for node=1:total_tracks_new(i)
        track_mrtrix_tck2.data{node+counter}=track_mrtrix_all(i).data{node};
        track_mrtrix_tsf.data{node+counter}=i*ones(size(track_mrtrix_all(i).data{node},1),1);
        track_mrtrix_tsf2.data{node+counter}=(1:100)';
        end
       counter=counter+total_tracks_new(i);
end
track_mrtrix_tck2.count=num2str(3);
track_mrtrix_tck2.total_count=num2str(counter);
write_mrtrix_tracks(track_mrtrix_tck2,[p '/along_tract_metrics/all.tck']);

track_mrtrix_tsf.count=num2str(1);
track_mrtrix_tsf.total_count=num2str(counter);
write_mrtrix_tsf(track_mrtrix_tsf,[p '/along_tract_metrics/all.tsf']);

track_mrtrix_tsf2.count=num2str(1);
track_mrtrix_tsf2.total_count=num2str(counter);
write_mrtrix_tsf(track_mrtrix_tsf2,[p '/along_tract_metrics/all_check_interp.tsf']);


%% calculate 4D Connectome 

connectome=zeros(max(ROI(:)),max(ROI(:)),length(scalar_maps),nb_nodes);
connectome_std=zeros(max(ROI(:)),max(ROI(:)),length(scalar_maps),nb_nodes);
total_tracks=zeros(1,length(combinations)); 

for cnt=1:length(combinations)
        total_tracks(cnt)=total_tracks_final_new(combinations(cnt,1),combinations(cnt,2));
end
total_tracks=[1,total_tracks];    % offset total_tracks with 1    

        for scal=1:length(scalar_maps)
          command=['tcksample ' p '/along_tract_metrics/all.tck ' p '/sDKIdu_' scalar_maps{scal} '_dki.nii ' p '/along_tract_metrics/' scalar_maps{scal} '.txt -force']  ;
          [~,~]=system(command);
          command=['tcksample ' p '/along_tract_metrics/all.tck ' p '/sDKIdu_' scalar_maps{scal} '_dki.nii ' p '/along_tract_metrics/' scalar_maps{scal} '.tsf -force']  ;
          [~,~]=system(command);

          sample=load([p '/along_tract_metrics/' scalar_maps{scal} '.txt' ]);
          if scal==2 || scal==3 || scal==4, sample=sample*1000; end % change units of diffusivity measures to um2/ms

            for j=1:length(combinations)
                if total_tracks(j+1) %skip if there where no tracts found between rois 
                connectome(combinations(j,1),combinations(j,2),scal,:)=mean(sample(sum(total_tracks(1:j)):sum(total_tracks(1:j+1))-1,:),1)';
                connectome_std(combinations(j,1),combinations(j,2),scal,:)=std(sample(sum(total_tracks(1:j)):sum(total_tracks(1:j+1))-1,:),1)';
                end
            end
            
        end
    save([p '/along_tract_metrics/Scalar_connectome.mat'],'connectome')
    save([p '/along_tract_metrics/Scalar_connectome_std.mat'],'connectome_std')

 
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

