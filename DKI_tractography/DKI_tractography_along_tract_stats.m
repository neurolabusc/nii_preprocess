%% calculate along tract measurements of MD/FA/MK between all pairs of a provided atlas (in diffusion space) 
function DKI_tractography_along_tract_stats(imgs,atlas,nb_nodes)
[p, n, x] = fileparts(imgs.T1);
if exist([p '/scalars_mean_' atlas '.mat'],'file')
    fprintf('Skipping DKI tractography with atlas: %s\n', atlas); 
    return
end
%create gm-wm interface for seeding purposes 
if ~exist([p '/5tt_' n x],'file')
    fprintf('Running mrtrix 5ttgen and 5tt2gmwmi to %s\n', n); 
command=['5ttgen fsl ' p '/wb' n x ' ' p '/w5tt_' n x ' -force -quiet'];
system(command)
command=['5tt2gmwmi ' p '/5tt_' n x ' ' p '/wgmwmi_' n x ' -force -quiet'];
system(command)
[p, n_DKI, x] = fileparts(imgs.DKI);
nfa=[p '/ns' n_DKI 'du_FAx' x];
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
    [p, n, x] = fileparts(imgs.DKI);
    hdr=spm_vol([ p '/' n  atlasext x]);
    ROI=spm_read_vols(hdr);  % read in atlas ROIs in native diffusion space
    index_ROI=[7 11 15 31 35 37 39 184 186 1 9 13 25 27 29 49 69 71 41 43]; % uncomment for only language specific and domain general ROIs  
   %index_ROI=[1:max(ROI(:))];
else
    atlasext = ['_roi_' atlas];
    [p, n, x] = fileparts(imgs.DKI);
    hdr=spm_vol([ p '/' n  atlasext x]);
    ROI=spm_read_vols(hdr);  % read in atlas ROIs in native diffusion space
    index_ROI=[1:max(ROI(:))]; % if JHU do only language specific and domain general ROIs  
    
end

%delete empty ROIs
count=1;
for i=1:length(index_ROI)
if sum(sum(sum((ROI==index_ROI(i)) & (gmwm>0))))
   index_ROI_no_empty(count)=index_ROI(i);
   count=count+1;
end
end


%%
scalar_maps={'du_FAx','du_MD','du_MK'}; % extension of parametric maps 

% initialize 
if mod(nb_nodes,2)==0, nb_nodes=nb_nodes+1;end % when using tie at center you need uneven number of points

meanInt=zeros(length(index_ROI_no_empty),length(index_ROI_no_empty),length(scalar_maps));
meanInt_final=zeros(max(ROI(:)),max(ROI(:)),length(scalar_maps));

stdInt=zeros(length(index_ROI_no_empty),length(index_ROI_no_empty),length(scalar_maps));
stdInt_final=zeros(max(ROI(:)),max(ROI(:)),length(scalar_maps));

total_tracks=zeros(length(index_ROI_no_empty),length(index_ROI_no_empty));
total_tracks_final=zeros(max(ROI(:)),max(ROI(:)));

scalar_mean=zeros(length(index_ROI_no_empty),length(index_ROI_no_empty),nb_nodes,length(scalar_maps)); 
scalar_mean_final=zeros(max(ROI(:)),max(ROI(:)),nb_nodes,length(scalar_maps)); 

scalar_sd=zeros(length(index_ROI_no_empty),length(index_ROI_no_empty),nb_nodes,length(scalar_maps));
scalar_sd_final=zeros(max(ROI(:)),max(ROI(:)),nb_nodes,length(scalar_maps));

track_mean_mrtrix=zeros(3,nb_nodes,length(index_ROI_no_empty),length(index_ROI_no_empty));
track_mean_mrtrix_final=zeros(3,nb_nodes,max(ROI(:)),max(ROI(:)));

track_mrtrix=cell(length(index_ROI_no_empty),10000,length(index_ROI_no_empty));
track_mrtrix_final=cell(max(ROI(:)),10000,max(ROI(:)));

pixel_size = abs(diag(hdr.mat))';
pixel_size = pixel_size(1:3); 
transform=hdr.mat; % transformation matrix used to calculate tranformation from world coordinate space to image space 

dim=hdr.dim;

mkdir([p '/temp']);



parfor i=1:length(index_ROI_no_empty)

    track_mean_mrtrix_temp=zeros(3,nb_nodes, length(index_ROI_no_empty));
    total_tracks_temp=zeros(1,length(index_ROI_no_empty));
    track_mrtrix_temp=cell(length(index_ROI_no_empty),10000); % 100000 is the max amount of tracks that can be found according to tracking settings 
        
    meanInt_temp=zeros(length(index_ROI_no_empty),length(scalar_maps));
    stdInt_temp=zeros(length(index_ROI_no_empty),length(scalar_maps));
    scalar_mean_temp=zeros(length(index_ROI_no_empty),nb_nodes,length(scalar_maps));
    scalar_std_temp=zeros(length(index_ROI_no_empty),nb_nodes,length(scalar_maps));
    
    % initialize for parfor .trk header 
    header=struct;
    header.voxel_size=pixel_size;
    header.dim=dim;
    header.n_scalars=0;

    for j=i+1:length(index_ROI_no_empty) % calculate only upper diagonal matrix 
        fprintf('%s: ROI %d to ROI %d\t',atlas,index_ROI_no_empty(i),index_ROI_no_empty(j));

        hdr_loop=spm_vol([ p '/' n  atlasext x]);            
        seed=(ROI==index_ROI_no_empty(i)) & (gmwm>0); hdr_loop.fname=[p '/temp/seed_' num2str(index_ROI_no_empty(i)) '_' num2str(index_ROI_no_empty(j)) '.nii']; spm_write_vol(hdr_loop,seed>0); % seed ROI 
        include=(ROI==index_ROI_no_empty(j))& (gmwm>0); hdr_loop.fname=[p '/temp/include_' num2str(index_ROI_no_empty(i)) '_' num2str(index_ROI_no_empty(j)) '.nii']; spm_write_vol(hdr_loop,include>0); % end ROI
        
%% Run MRtrix tractography - deterministic, dODF based, angular threshold 60, FA>0.1
% seed to end
command=['tckgen -algorithm SD_STREAM -seed_image ' p '/temp/seed_' num2str(index_ROI_no_empty(i)) '_' num2str(index_ROI_no_empty(j)) '.nii -cutoff 0.1 -include ' p '/temp/include_' num2str(index_ROI_no_empty(i)) '_' num2str(index_ROI_no_empty(j)) '.nii ' p '/SH_coeff.nii ' p '/temp/seed_include_' num2str(index_ROI_no_empty(i)) '_' num2str(index_ROI_no_empty(j)) '.tck -seeds 10000 -select 10000 -angle 60 -stop -force -seed_unidirectional -quiet'];
system(command);
% end to seed 
command=['tckgen -algorithm SD_STREAM -seed_image ' p '/temp/include_' num2str(index_ROI_no_empty(i)) '_' num2str(index_ROI_no_empty(j)) '.nii -cutoff 0.1 -include ' p '/temp/seed_' num2str(index_ROI_no_empty(i)) '_' num2str(index_ROI_no_empty(j)) '.nii ' p '/SH_coeff.nii ' p '/temp/include_seed_' num2str(index_ROI_no_empty(i)) '_' num2str(index_ROI_no_empty(j)) '.tck -seeds 10000 -select 10000 -angle 60 -stop -force -seed_unidirectional -quiet'];
system(command);
% combine both results into one .tck 
command=['tckedit ' p '/temp/seed_include_' num2str(index_ROI_no_empty(i)) '_' num2str(index_ROI_no_empty(j)) '.tck ' p '/temp/include_seed_' num2str(index_ROI_no_empty(i)) '_' num2str(index_ROI_no_empty(j)) '.tck ' p '/temp/combined_' num2str(index_ROI_no_empty(i)) '_' num2str(index_ROI_no_empty(j)) '.tck -force -quiet'];
system(command);
%% Calculate along tract metrics (cite John Colby work: https://www.sciencedirect.com/science/article/pii/S1053811911012833?via%3Dihub ) 

tracks_trk=struct;

tracks = read_mrtrix_tracks ([ p '/temp/combined_' num2str(index_ROI_no_empty(i)) '_' num2str(index_ROI_no_empty(j)) '.tck']); % tracks from mrtrix are in world coordinates (voxels) 
tracks_mm=tracks; % initialize the matrix needed to transform tracks from voxels to mm 
total_tracks_temp(j)=length(tracks.data); % calculate how many tracts were found between the two ROIs 
fprintf('Total Tracks: %d\n', total_tracks_temp(j));
 
if (total_tracks_temp(j)~=0) && (total_tracks_temp(j)>5) % do along tract measurements if there are more than 5 streamlines between ROIs 

        for nb_tracks=1:length(tracks.data)
          int=transform\[tracks.data{nb_tracks} ones(length(tracks.data{nb_tracks}),1)]'; % transform to image space
          tracks_mm.data{nb_tracks}=int(1:3,:)'.*repmat(pixel_size,length(tracks.data{nb_tracks}),1);  % convert to mm 
        end

        for nb_tracks=1:length(tracks.data) % change from .tck format to .trk format 
        tracks_trk(nb_tracks).matrix=tracks_mm.data{nb_tracks};
        tracks_trk(nb_tracks).nPoints=length(tracks_mm.data{nb_tracks});
        end
        header.n_count=length(tracks.data);

    interpolated_trk=trk_interp(tracks_trk,nb_nodes,[],1); % interpolate tracts into 100 nodes, will make the number of nodes uneven because of the tie at center option 
    tracks_interp     = trk_flip(header, interpolated_trk, [1 1 1]);

    % calculate weights for diffusion metrics based on distance to
    % geometric average (tracts further away from the mean contribute less
    % to the average metric)(borrowed strategy from AFQ/mrDiffusion: dtiFIberGroupPropertyWeightedAverage.m & AFQ_FiberTractGaussian.m) 
    weights=[]; 
    weights_param=[];
    track_mean=mean(tracks_interp,3);
        for node=1:nb_nodes
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
    weights=(weights); % parfor loops needs this??? 
    index_keep=(sum(weights<4)==nb_nodes); % remove tracks that are more than 4 standard deviations away from tract mean 
    weights_param=(weights_param); % parfor loops needs this??? 
    tracks_interp=tracks_interp(:,:,index_keep);
    
    weightsNormalized = weights_param(:,index_keep)./(repmat(sum(weights_param(:,index_keep), 2), [1 sum(index_keep)])); % normalize weights so the sum of all nodes is 1 (weighted average) 

    tracks_interp_str = trk_restruc(tracks_interp);
    header.n_count=length(tracks_interp_str); %recalculate header count because of removed tracks 

    % save tracks
    
    %mean track 
    track_mean = mean(tracks_interp, 3);
    int=(transform*[(track_mean./repmat(pixel_size,length(track_mean),1)) ones(length(track_mean),1)]'); % convert to .tck format (world coordinates in voxels) 
    track_mean_mrtrix_temp(:,:,j)=int(1:3,:);
    
    total_tracks_temp(j)=sum(index_keep); %update total tracks after filtering

    % save all bundles in .tck
        for nb=1:sum(index_keep)
        int=(transform*[(tracks_interp(:,:,nb)./repmat(pixel_size,length(track_mean),1)) ones(length(track_mean),1)]'); % convert to .tck format (world coordinates in voxels) 
        track_mrtrix_temp{j,nb}=int(1:3,:)';
        end

    for sc_map=1:length(scalar_maps) 
    % SCALAR MAPS 
    hdr=spm_vol([p '/s' n  scalar_maps{sc_map} '.nii']);
    map=spm_read_vols(hdr);
    if sc_map==2, map=map*1000;end % convert MD to more conventional um2/ms 

    [meanInt_temp(j,sc_map), stdInt_temp(j,sc_map), ~] = trk_stats(header, tracks_interp_str, map, 'nearest');
    [header_sc, tracks_sc] = trk_add_sc(header, tracks_interp_str, map, scalar_maps{sc_map});
    
    %initalize
    scalars = zeros(tracks_sc(1).nPoints, header_sc.n_count, header_sc.n_scalars);
    scalars_nw = zeros(tracks_sc(1).nPoints, header_sc.n_count, header_sc.n_scalars);

    fw = weightsNormalized(:);

        for sc=1:header_sc.n_scalars
            mat_long        = cat(1, tracks_sc.matrix);
            scalars(:,:,sc)  = reshape(mat_long(:,4).*fw, tracks_sc(1).nPoints, header_sc.n_count, header_sc.n_scalars);
            scalars_nw(:,:,sc)  = reshape(mat_long(:,4), tracks_sc(1).nPoints, header_sc.n_count, header_sc.n_scalars);
        end

    scalar_mean_temp(j,:,sc_map) = squeeze(sum(scalars, 2)); % sum because they are weighted scalars 
    scalar_std_temp(j,:,sc_map)  = squeeze(nanstd(scalars_nw, 0, 2));
    end
    
end
    end
    total_tracks(i,:)=total_tracks_temp;
    meanInt(i,:,:)=meanInt_temp;
    stdInt(i,:,:)=stdInt_temp;
    scalar_mean(i,:,:,:)=scalar_mean_temp;
    scalar_sd(i,:,:,:)=scalar_std_temp;

    track_mean_mrtrix(:,:,:,i)=track_mean_mrtrix_temp;
    track_mrtrix(:,:,i)=track_mrtrix_temp;
end

%reindex because parfor didnt let me index it properly from the beginning
for i=1:length(index_ROI_no_empty)
    for j=1:length(index_ROI_no_empty)
    total_tracks_final(index_ROI_no_empty(i),index_ROI_no_empty(j))=total_tracks(i,j);
    meanInt_final(index_ROI_no_empty(i),index_ROI_no_empty(j),:)=meanInt(i,j,:);
    stdInt_final(index_ROI_no_empty(i),index_ROI_no_empty(j),:)=stdInt(i,j,:);
    scalar_mean_final(index_ROI_no_empty(i),index_ROI_no_empty(j),:,:)=scalar_mean(i,j,:,:);
    scalar_sd_final(index_ROI_no_empty(i),index_ROI_no_empty(j),:,:)=scalar_sd(i,j,:,:);

    track_mean_mrtrix_final(:,:,index_ROI_no_empty(j),index_ROI_no_empty(i))=track_mean_mrtrix(:,:,j,i);
    track_mrtrix_final(index_ROI_no_empty(j),:,index_ROI_no_empty(i))=track_mrtrix(j,:,i);
    end
end

% save along tract stats as a structure, note that this is written
% specifically for three maps change this when adding metrics 
scalar_mean_struc=struct(scalar_maps{1},scalar_mean_final(:,:,:,1),scalar_maps{2},scalar_mean_final(:,:,:,2),scalar_maps{3},scalar_mean_final(:,:,:,3));
scalar_sd_struc=struct(scalar_maps{1},scalar_sd(:,:,:,1),scalar_maps{2},scalar_sd(:,:,:,2),scalar_maps{3},scalar_sd(:,:,:,3));

save([p '/scalars_mean_' atlas '.mat'],'scalar_mean_struc');
save([p '/scalars_sd_' atlas '.mat'],'scalar_sd_struc');
save([p '/mean_' atlas '.mat'],'meanInt_final') %mean tract values (FA,MD,MK)
save([p '/std_' atlas '.mat'],'stdInt_final') %std tract values (FA,MD,MK)
save([p '/total_tracks_' atlas '.mat'],'total_tracks_final');
save([p '/track_mrtrix_final_' atlas '.mat'],'track_mrtrix_final');
save([p '/track_mean_mrtrix_final_' atlas '.mat'],'track_mean_mrtrix_final');


tracks = read_mrtrix_tracks ([ p '/temp/combined_1_9.tck']); % tracks from mrtrix are in world coordinates (voxels) 
track_mean_mrtrix_tck=tracks;
for i=1:length(index_ROI_no_empty)
    for j=1:length(index_ROI_no_empty)
track_mean_mrtrix_tck.data{(i-1)*length(index_ROI_no_empty)+j}=track_mean_mrtrix_final(:,:,index_ROI_no_empty(i),index_ROI_no_empty(j))';
    end
end
track_mean_mrtrix_tck.count=length(index_ROI_no_empty)*length(index_ROI_no_empty);
track_mean_mrtrix_tck.total_count=length(index_ROI_no_empty)*length(index_ROI_no_empty);
write_mrtrix_tracks(track_mean_mrtrix_tck,[p '/mean_all_' atlas '.tck']);

tracks = read_mrtrix_tracks ([ p '/temp/combined_1_9.tck']); % tracks from mrtrix are in world coordinates (voxels) 
track_mrtrix_tck=tracks;
counter=0;
for i=1:length(index_ROI_no_empty)
    for j=1:length(index_ROI_no_empty)
        for node=1:total_tracks_final(index_ROI_no_empty(i),index_ROI_no_empty(j))
        track_mrtrix_tck.data{node+counter}=track_mrtrix_final{index_ROI_no_empty(j),node,index_ROI_no_empty(i)};
        end
       counter=counter+total_tracks_final(index_ROI_no_empty(i),index_ROI_no_empty(j));
    end
end
track_mrtrix_tck.count=counter;
track_mrtrix_tck.total_count=counter;
write_mrtrix_tracks(track_mrtrix_tck,[p '/all_' atlas '.tck']);

function oldNormSub(src, tar, smoref, reg, interp)
%coregister T2 to match T1 image, apply to lesion
if isempty(src) || isempty(tar), return; end;
if ~exist('smoref','var'), smoref = 0; end;
if ~exist('reg','var'), reg = 1; end;
if ~exist('interp','var'), interp = 1; end;
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

