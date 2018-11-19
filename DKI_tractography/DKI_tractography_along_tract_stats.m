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
command=['5ttgen fsl ' p '/wwb' n x ' ' p '/5tt_' n x ' -force -quiet']; % need segmentation of T1 in diffusion space, this step could probably be optimized
system(command)
command=['5tt2gmwmi ' p '/5tt_' n x ' ' p '/gmwmi_' n x ' -force -quiet'];
system(command)
spm_reslice({[p '/wwb' n x],[p '/gmwmi_' n x]}) % get in same space as wwbT1 (not sure why it's not in that space to begin with?) 

else
    fprintf('Skipping generating WM-GM interface: found %s\n', ['gmwmi_' n x]);
end

%%
% threshold gmwm interface at 50% chance 
gmwm=spm_read_vols(spm_vol([p '/rgmwmi_' n x]));
gmwm=gmwm>0.5; % threshold was picked arbitrary, could likely be optimized 

if strcmpi(atlas,'jhu')
    atlasext = '_roi';
else
   atlasext = ['_roi_' atlas];
end

[p, n, x] = fileparts(imgs.DKI);

hdr=spm_vol([ p '/' n  atlasext x]);
ROI=spm_read_vols(hdr);  % read in atlas ROIs in native diffusion space
%%
scalar_maps={'du_FAx','du_MD','du_MK'}; % extension of parametric maps 

% initialize 
meanInt=zeros(max(ROI(:)),max(ROI(:)));
stdInt=zeros(max(ROI(:)),max(ROI(:)));
total_tracks=zeros(max(ROI(:)),max(ROI(:)));
scalar_mean=zeros(max(ROI(:)),max(ROI(:)),nb_nodes+1,length(scalar_maps)); % when nb_nodes is uneven, remove +1 
scalar_sd=zeros(max(ROI(:)),max(ROI(:)),nb_nodes+1,length(scalar_maps)); % when nb_nodes is uneven, remove +1 

pixel_size = abs(diag(hdr.mat))';
pixel_size = pixel_size(1:3); 
transform=hdr.mat; % transformation matrix used to calculate tranformation from world coordinate space to image space 

% initialize .trk header 
header.voxel_size=pixel_size;
header.dim=hdr.dim;
header.n_scalars=0;
    
for i=1:max(ROI(:))
    for j=i:max(ROI(:)) % calculate only upper diagonal matrix 
        fprintf('%s: ROI %d to ROI %d\t',atlas,i,j);
        if i~=j % don't seed and end in the same ROI
            
seed=(ROI==i) & (gmwm>0); hdr.fname=[p '/seed.nii']; spm_write_vol(hdr,seed>0); % seed ROI 
include=(ROI==j)& (gmwm>0); hdr.fname=[p '/include.nii']; spm_write_vol(hdr,include>0); % end ROI

%% Run MRtrix tractography - deterministic, dODF based, angular threshold 60, FA>0.1
% seed to end
command=['tckgen -algorithm SD_STREAM -seed_image ' p '/seed.nii -cutoff 0.1 -include ' p '/include.nii ' p '/SH_coeff.nii ' p '/seed_include.tck -seeds 10000 -select 10000 -angle 60 -stop -force -seed_unidirectional -quiet'];
system(command);
% end to seed 
command=['tckgen -algorithm SD_STREAM -seed_image ' p '/include.nii -cutoff 0.1 -include ' p '/seed.nii ' p '/SH_coeff.nii ' p '/include_seed.tck -seeds 10000 -select 10000 -angle 60 -stop -force -seed_unidirectional -quiet'];
system(command);
% combine both results into one .tck 
command=['tckedit ' p '/seed_include.tck ' p '/include_seed.tck ' p '/combined.tck -force -quiet'];
system(command);
%% Calculate along tract metrics (cite John Colby work: https://www.sciencedirect.com/science/article/pii/S1053811911012833?via%3Dihub ) 
clear tracks_trk
tracks = read_mrtrix_tracks ([ p '/combined.tck']); % tracks from mrtrix are in world coordinates (voxels) 
tracks_mm=tracks; % initialize the matrix needed to transform tracks from voxels to mm 


total_tracks(i,j)=length(tracks.data); % calculate how many tracts were found between the two ROIs 
fprintf('Total Tracks: %d\n', total_tracks(i,j));
 
if (total_tracks(i,j)~=0) && (total_tracks(i,j)>5) % do along tract measurements if there are more than 5 streamlines between ROIs 

        for nb_tracks=1:length(tracks.data)
          int=inv(transform)*[tracks.data{nb_tracks} ones(length(tracks.data{nb_tracks}),1)]'; % transform to image space
          tracks_mm.data{nb_tracks}=int(1:3,:)'.*repmat(pixel_size,length(tracks.data{nb_tracks}),1);  % convert to mm 
        end

        for nb_tracks=1:length(tracks.data) % change from .tck format to .trk format 
        tracks_trk(nb_tracks).matrix=tracks_mm.data{nb_tracks};
        tracks_trk(nb_tracks).nPoints=length(tracks_mm.data{nb_tracks});
        end
        header.n_count=length(tracks.data);



    interpolated_trk=trk_interp(tracks_trk,nb_nodes,[],1); % interpolate tracts into 100 nodes, will make the number of nodes uneven because of the tie at center option 
    if mod(nb_nodes,2)==0, nb_nodes=nb_nodes+1;end % when using tie at center you need uneven number of points
    tracks_interp     = trk_flip(header, interpolated_trk, [1 1 1]);

    % calculate weights for diffusion metrics based on distance to
    % geometric average (tracts further away from the mean contribute less
    % to the average metric)(borrowed strategy from AFQ/mrDiffusion: dtiFIberGroupPropertyWeightedAverage.m & AFQ_FiberTractGaussian.m) 
    clear weights weights_param
    track_mean=mean(tracks_interp,3);
        for node=1:nb_nodes
        int=tril(cov(permute(tracks_interp(node,:,:),[3 2 1])));
        cov_matrix(:,node)=int([1:3 5:6 9]');
        sigma=[cov_matrix(1:3,node)'
               0 cov_matrix(4:5,node)'
               0 0 cov_matrix(6,node)'];
        sigma=sigma+sigma'-diag(diag(sigma));
        mu=track_mean(node,:);
        d=bsxfun(@minus,squeeze(tracks_interp(node,:,:))',mu);
        % mahalanobis distance of each point on each fiber from the tract core  
        weights(node,:)=sqrt(dot(d/(sigma),d,2))';
        weights_param(node,:)=mvnpdf(squeeze(tracks_interp(node,:,:))',mu,sigma)'; 
        end

    index_keep=(sum(weights<4)==nb_nodes); % remove tracks that are more than 4 standard deviations away from tract mean 
    tracks_interp=tracks_interp(:,:,index_keep);

    weightsNormalized = weights_param(:,index_keep)./(repmat(sum(weights_param(:,index_keep), 2), [1 sum(index_keep)])); % normalize weights so the sum of all nodes is 1 (weighted average) 
    tracks_interp_str = trk_restruc(tracks_interp);
    header.n_count=length(tracks_interp_str); %recalculate header count because of removed tracks 

    % save tracks
    
    %mean track 
    track_mean = mean(tracks_interp, 3);
    track_mean_mrtrix=tracks;
    int=(transform*[(track_mean./repmat(pixel_size,length(track_mean),1)) ones(length(track_mean),1)]'); % convert to .tck format (world coordinates in voxels) 
    track_mean_mrtrix.data={int(1:3,:)'};
    track_mean_mrtrix.count='1';
    track_mean_mrtrix.total_count='1';
    
    if ~exist([p '/mean_all.tck'],'file')
    write_mrtrix_tracks(track_mean_mrtrix,[p '/mean_all.tck']);
    else
    write_mrtrix_tracks(track_mean_mrtrix,[p '/mean.tck']); % save mean.tck containing the mean bundle between (i,j) 
    command=['tckedit ' p '/mean_all.tck ' p '/mean.tck '  p '/mean_all_int.tck -force -quiet']; % save all mean bundles in mean_all.tck 
    system(command);
    movefile([p '/mean_all_int.tck'],[p '/mean_all.tck']); % needed this workaround to concatenate .tck's 
    end
    
    % save all bundles in .tck
        track_mrtrix=tracks;
        for nb=1:sum(index_keep)
        int=(transform*[(tracks_interp(:,:,nb)./repmat(pixel_size,length(track_mean),1)) ones(length(track_mean),1)]'); % convert to .tck format (world coordinates in voxels) 
        track_mrtrix.data{nb}=int(1:3,:)';
        end

    track_mrtrix.count=num2str(sum(index_keep));
    track_mrtrix.total_count=num2str(sum(index_keep));
    track_mrtrix.data(sum(index_keep)+1:end)=[];

    if ~exist([p '/all_' atlas '.tck'],'file')
    write_mrtrix_tracks(track_mrtrix,[p '/all_' atlas '.tck']); % initialize all_atlas.tck
    else
    write_mrtrix_tracks(track_mrtrix,[p '/int.tck']); % all tracts between roi i and j 
    command=['tckedit ' p '/all.tck ' p '/int.tck '  p '/all_int.tck -force -quiet']; 
    system(command);
    movefile([p '/all_int.tck'],[p 'all_' atlas '.tck']); % concatenate tracts in all_atlas.tck 
    end

    % calculate scalars along the tract coordinates
    clear header_sc_trk tracks_sc_trk
    
    header.id_string='TRACK ';
    header.origin=[0 0 0];
    header.n_properties=0;
    header.property_name= repmat(blanks(20),[10,1]);
    header.vox_to_ras=transform; 
    header.reserved=blanks(444)'; 
    header.voxel_order='LAS ';
    header.pad2='LAS ';
    header.image_orientation_patient= [ 1 0 0 0 1 0];
    header.pad1='  ';
    header.invert_x=0;
    header.invert_y=0;
    header.invert_z=0;
    header.swap_xy=0;
    header.swap_yz=0;
    header.swap_zx=0;
    header.version=2;
    header.hdr_size=1000;
    header.scalar_name=repmat(blanks(20),[10,1]);
    
    for sc_map=1:length(scalar_maps) 
    % SCALAR MAPS 
    hdr=spm_vol([p '/s' n  scalar_maps{sc_map} '.nii']);
    map=spm_read_vols(hdr);
    if sc_map==2, map=map*1000;end % convert MD to more conventional um2/ms 

    [meanInt(i,j,sc_map), stdInt(i,j,sc_map), ~] = trk_stats(header, tracks_interp_str, map, 'nearest');
    [header_sc, tracks_sc] = trk_add_sc(header, tracks_interp_str, map, scalar_maps{sc_map});
        if sc_map==1, header_sc_trk=header_sc;tracks_sc_trk=tracks_sc; else [header_sc_trk, tracks_sc_trk] = trk_add_sc(header_sc_trk, tracks_sc_trk, map, scalar_maps{sc_map}); end

    %initalize
    scalars = zeros(tracks_sc(1).nPoints, header_sc.n_count, header_sc.n_scalars);
    scalars_nw = zeros(tracks_sc(1).nPoints, header_sc.n_count, header_sc.n_scalars);

    fw = weightsNormalized(:);

        for sc=1:header_sc.n_scalars
            mat_long        = cat(1, tracks_sc.matrix);
            scalars(:,:,sc)  = reshape(mat_long(:,4).*fw, tracks_sc(1).nPoints, header_sc.n_count, header_sc.n_scalars);
            scalars_nw(:,:,sc)  = reshape(mat_long(:,4), tracks_sc(1).nPoints, header_sc.n_count, header_sc.n_scalars);
        end

    scalar_mean(i,j,1:nb_nodes,sc_map) = squeeze(sum(scalars, 2)); % sum because they are weighted scalars 
    scalar_sd(i,j,1:nb_nodes,sc_map)  = squeeze(nanstd(scalars_nw, 0, 2));
    end
    
    %trk_write(header_sc_trk,tracks_sc_trk,[p '/all_' atlas '.trk']); can
    %be used in the future to output.trk files ( needs debugging)
    
end
        end

    end
end
% save along tract stats as a structure, note that this is written
% specifically for three maps change this when adding metrics 
scalar_mean_struc=struct(scalar_maps{1},scalar_mean(:,:,:,1),scalar_maps{2},scalar_mean(:,:,:,2),scalar_maps{3},scalar_mean(:,:,:,3));
scalar_sd_struc=struct(scalar_maps{1},scalar_sd(:,:,:,1),scalar_maps{2},scalar_sd(:,:,:,2),scalar_maps{3},scalar_sd(:,:,:,3));

save([p '/scalars_mean_' atlas '.mat'],'scalar_mean_struc')
save([p '/scalars_sd' atlas '.mat'],'scalar_sd_struc')
save([p '/mean_' atlas '.mat'],'meanInt') %mean tract values (FA,MD,MK)
save([p '/std' atlas '.mat'],'meanInt') %std tract values (FA,MD,MK)
save([p '/total_tracks' atlas '.mat'],'total_tracks')


