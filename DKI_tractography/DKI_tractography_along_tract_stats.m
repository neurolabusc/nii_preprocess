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
gmwm=gmwm>0.5;

if strcmpi(atlas,'jhu')
    atlasext = '_roi';
else
   atlasext = ['_roi_' atlas];
end

[p, n, x] = fileparts(imgs.DKI);

hdr=spm_vol([ p '/' n  atlasext x]);
ROI=spm_read_vols(hdr);
%%
scalar_maps={'du_FAx','du_MD','du_MK'};

meanInt=zeros(max(ROI(:)),max(ROI(:)));
stdInt=zeros(max(ROI(:)),max(ROI(:)));
total_tracks=zeros(max(ROI(:)),max(ROI(:)));
scalar_mean=zeros(max(ROI(:)),max(ROI(:)),101,length(scalar_maps));
scalar_sd=zeros(max(ROI(:)),max(ROI(:)),101,length(scalar_maps));

pixel_size = abs(diag(hdr.mat))';
pixel_size = pixel_size(1:3); 
transform=hdr.mat;

% initialize .trk header 
header.voxel_size=pixel_size;
header.dim=hdr.dim;
header.n_scalars=0;
    
for i=1:max(ROI(:))
    for j=i:max(ROI(:))
        fprintf('ROI %d to ROI %d\t', i,j);
        if i~=j
            
seed=(ROI==i) & (gmwm>0); hdr.fname=[p '/seed.nii']; spm_write_vol(hdr,seed>0);
include=(ROI==j)& (gmwm>0); hdr.fname=[p '/include.nii']; spm_write_vol(hdr,include>0);

%%

command=['tckgen -algorithm SD_STREAM -seed_image ' p '/seed.nii -cutoff 0.1 -include ' p '/include.nii ' p '/SH_coeff.nii ' p '/seed_include.tck -seeds 10000 -select 10000 -angle 60 -stop -force -seed_unidirectional -quiet'];
system(command);
command=['tckgen -algorithm SD_STREAM -seed_image ' p '/include.nii -cutoff 0.1 -include ' p '/seed.nii ' p '/SH_coeff.nii ' p '/include_seed.tck -seeds 10000 -select 10000 -angle 60 -stop -force -seed_unidirectional -quiet'];
system(command);
command=['tckedit ' p '/seed_include.tck ' p '/include_seed.tck ' p '/combined.tck -force -quiet'];
system(command);
%%
clear tracks_trk
tracks = read_mrtrix_tracks ([ p '/combined.tck']); % tracks from mrtrix are in world coordinate system in voxels
tracks_mm=tracks; 


total_tracks(i,j)=length(tracks.data);
fprintf('Total Tracks: %d\n', total_tracks(i,j));
 
if (total_tracks(i,j)~=0) && (total_tracks(i,j)>5) 

        for nb_tracks=1:length(tracks.data)
          int=inv(transform)*[tracks.data{nb_tracks} ones(length(tracks.data{nb_tracks}),1)]'; % transform to image space
          tracks_mm.data{nb_tracks}=int(1:3,:)'.*repmat(pixel_size,length(tracks.data{nb_tracks}),1);  % convert to mm 
        end

        for nb_tracks=1:length(tracks.data)
        tracks_trk(nb_tracks).matrix=tracks_mm.data{nb_tracks};
        tracks_trk(nb_tracks).nPoints=length(tracks_mm.data{nb_tracks});
        end


    header.n_count=length(tracks.data);
    interpolated_trk=trk_interp(tracks_trk,nb_nodes,[],1); % 
    if mod(nb_nodes,2)==0, nb_nodes=nb_nodes+1;end
    tracks_interp     = trk_flip(header, interpolated_trk, [1 1 1]);

    % 
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
        weights(node,:)=sqrt(dot(d/(sigma),d,2))';
        weights_param(node,:)=mvnpdf(squeeze(tracks_interp(node,:,:))',mu,sigma)'; % calculate weights based on distance to geometric average 
        end

    index_keep=(sum(weights<4)==nb_nodes); % remove tracks that are more than 4 standard deviations away from mean 
    tracks_interp=tracks_interp(:,:,index_keep);

    weightsNormalized = weights_param(:,index_keep)./(repmat(sum(weights_param(:,index_keep), 2), [1 sum(index_keep)]));
    tracks_interp_str = trk_restruc(tracks_interp);
    header.n_count=length(tracks_interp_str); %recalculate header count because of removed tracks 

    % save tracks
    
    %mean track 
    track_mean = mean(tracks_interp, 3);
    track_mean_mrtrix=tracks;
    int=(transform*[(track_mean./repmat(pixel_size,length(track_mean),1)) ones(length(track_mean),1)]');
    track_mean_mrtrix.data={int(1:3,:)'};
    track_mean_mrtrix.count='1';
    track_mean_mrtrix.total_count='1';
    
    if ~exist([p '/mean_all.tck'],'file')
    write_mrtrix_tracks(track_mean_mrtrix,[p '/mean_all.tck']);
    else
    write_mrtrix_tracks(track_mean_mrtrix,[p '/mean.tck']);
    command=['tckedit ' p '/mean_all.tck ' p '/mean.tck '  p '/mean_all_int.tck -force -quiet'];
    system(command);
    movefile([p '/mean_all_int.tck'],[p '/mean_all.tck']);
    end
    
    % save whole brain 
        track_mrtrix=tracks;
        for nb=1:sum(index_keep)
        int=(transform*[(tracks_interp(:,:,nb)./repmat(pixel_size,length(track_mean),1)) ones(length(track_mean),1)]');
        track_mrtrix.data{nb}=int(1:3,:)';
        end

    track_mrtrix.count=num2str(sum(index_keep));
    track_mrtrix.total_count=num2str(sum(index_keep));
    track_mrtrix.data(sum(index_keep)+1:end)=[];

    if ~exist([p '/all_' atlas '.tck'],'file')
    write_mrtrix_tracks(track_mrtrix,[p '/all_' atlas '.tck']);
    else
    write_mrtrix_tracks(track_mrtrix,[p '/int.tck']);
    command=['tckedit ' p '/all.tck ' p '/int.tck '  p '/all_int.tck -force -quiet'];
    system(command);
    movefile([p '/all_int.tck'],[p 'all_' atlas '.tck']);
    end


    for sc_map=1:length(scalar_maps)
    % SCALAR MAPS 
    hdr=spm_vol([p '/s' n  scalar_maps{sc_map} '.nii']);
    map=spm_read_vols(hdr);

    [meanInt(i,j), stdInt(i,j), ~] = trk_stats(header, tracks_trk, map, 'nearest');
    [header_sc, tracks_sc] = trk_add_sc(header, tracks_interp_str, map, scalar_maps{sc_map});

    scalars = zeros(tracks_sc(1).nPoints, header_sc.n_count, header_sc.n_scalars);
    scalars_nw = zeros(tracks_sc(1).nPoints, header_sc.n_count, header_sc.n_scalars);

    fw = weightsNormalized(:);

        for sc=1:header_sc.n_scalars
            mat_long        = cat(1, tracks_sc.matrix);
            scalars(:,:,sc)  = reshape(mat_long(:,4).*fw, tracks_sc(1).nPoints, header_sc.n_count, header_sc.n_scalars);
            scalars_nw(:,:,sc)  = reshape(mat_long(:,4), tracks_sc(1).nPoints, header_sc.n_count, header_sc.n_scalars);
        end

    scalar_mean(i,j,1:nb_nodes,sc_map) = squeeze(sum(scalars, 2));
    scalar_sd(i,j,1:nb_nodes,sc_map)  = squeeze(nanstd(scalars_nw, 0, 2));
    end
end
        end

    end
end

command=['tcksample ' p '/mean_all.tck ' p '/s' n 'du_FAx.nii mean_color.tsf -force'];
system(command)

command=['tcksample ' p '/all.tck ' p '/s' n 'du_FAx.nii all_color.tsf -force'];
system(command)

scalar_mean_struc=struct(scalar_maps{1},scalar_mean(:,:,:,1),scalar_maps{2},scalar_mean(:,:,:,2),scalar_maps{3},scalar_mean(:,:,:,3));
scalar_sd_struc=struct(scalar_maps{1},scalar_sd(:,:,:,1),scalar_maps{2},scalar_sd(:,:,:,2),scalar_maps{3},scalar_sd(:,:,:,3));

save([p '/scalars_mean_' atlas '.mat'],'scalar_mean_struc')
save([p '/scalars_sd' atlas '.mat'],'scalar_sd_struc')
save([p '/total_tracks' atlas '.mat'],'total_tracks')


