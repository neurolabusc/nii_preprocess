
function DKI_tractometry(imgs,atlas,nb_nodes,scalar_maps)
%% initialize
index_ROI=[1:2:93,96,97,99,113,115,167,182:2:188]; % all left gray matter jhu
combinations=nchoosek(index_ROI,2);
[p, n, ~] = fileparts(imgs.T1);

mkdir([p '/along_tract_metrics'])
hdr=spm_vol([p '/DKI_roi.nii']);
ROI=spm_read_vols(hdr); 
pixel_size = abs(diag(hdr.mat))';
pixel_size = pixel_size(1:3); 
dim=hdr.dim;
transform=hdr.mat; % transformation matrix used to calculate tranformation from world coordinate space to image space 

if exist([p '/all_' atlas '_excl.tck'],'file')
track_mrtrix = read_mrtrix_tracks([p '/all_' atlas '_excl.tck']); % check if lesion was excluded  
load([p '/total_tracks_' atlas '_excl.mat']);
else 
track_mrtrix = read_mrtrix_tracks([p '/all_' atlas '.tck']); % assume lesion was not excluded
load([p '/total_tracks_' atlas '.mat']);
end


gmwm=spm_read_vols(spm_vol([p '/gmwmi_' n '.nii']));   
total_tracks=zeros(1,length(combinations)+1); 

 % offset total_tracks with 0 
for cnt=1:length(combinations)
    total_tracks(cnt+1)=total_tracks_final(combinations(cnt,1),combinations(cnt,2));
end
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
        int=hdr.mat\[track_mrtrix.data{node+sum(total_tracks(1:j))} ones(length(track_mrtrix.data{node+sum(total_tracks(1:j))}),1)]'; % tracks in image space
        tracks_trk(node).matrix=int(1:3,:)'.*repmat(pixel_size,length(track_mrtrix.data{node+sum(total_tracks(1:j))}),1); %  image space in mm 
        tracks_trk(node).nPoints=length(track_mrtrix.data{node+sum(total_tracks(1:j))});
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
         
    interpolated_trk_flip(:,:,B)=interpolated_trk(fliplr(1:end),:,B); % flip tracks that start in ROI1 
   
    % write tracts
    track_mrtrix_tck=track_mrtrix;
    track_mrtrix_tsf=track_mrtrix;
    track_mrtrix_tck.data=[];
    track_mrtrix_tsf.data=[];

     for nb=1:sum(index)
     int=(transform*[(interpolated_trk_flip(:,:,nb)./repmat(pixel_size,size(interpolated_trk_flip,1),1)) ones(size(interpolated_trk_flip,1),1)]'); % convert to .tck format (world coordinates in voxels) 
     track_mrtrix_tck.data{nb}=int(1:3,:)';
     track_mrtrix_tsf.data{nb}=(1:size(interpolated_trk_flip,1))';
     end

% uncomment if you want to output all tracks seperately 
     
%     track_mrtrix_tck.count=1;
%     track_mrtrix_tck.total_count=sum(index);
%     write_mrtrix_tracks(track_mrtrix_tck,[subdir '/BL/along_tract_metrics/ROI' num2str(combinations(j,1)) '_ROI' num2str(combinations(j,2)) '.tck']);
%     track_mrtrix_tsf.count=1;
%     track_mrtrix_tsf.total_count=sum(index);
%     write_mrtrix_tsf(track_mrtrix_tsf,[subdir '/BL/along_tract_metrics/ROI' num2str(combinations(j,1)) '_ROI' num2str(combinations(j,2)) '.tsf']);
%    
   track_mrtrix_all(j)=track_mrtrix_tck;
   total_tracks_new(j)=length(track_mrtrix_tck.data); % calculate how many tracts were found between the two ROIs  

   
     % uncomment if you want to calculate the mean geometry and output all
     % of them
     
%     mean_trk=mean(interpolated_trk_flip,3)   ;
%    
%     track_mrtrix_tck=track_mrtrix;
%     track_mrtrix_tsf=track_mrtrix;
%     track_mrtrix_tck.data=[];
%     track_mrtrix_tsf.data=[];
% 
%         nb=1;
%         int=(transform*[(mean_trk./repmat(pixel_size,size(interpolated_trk_flip,1),1)) ones(size(interpolated_trk_flip,1),1)]'); % convert to .tck format (world coordinates in voxels) 
%         track_mrtrix_tck.data{nb}=int(1:3,:)';
%         track_mrtrix_tsf.data{nb}=(1:size(interpolated_trk_flip,1))';
%            
  
%     track_mrtrix_tck.count=1;
%     track_mrtrix_tck.total_count=1;
%     write_mrtrix_tracks(track_mrtrix_tck,[subdir '/BL/along_tract_metrics/mean_ROI' num2str(combinations(j,1)) '_ROI' num2str(combinations(j,2)) '.tck']);
%     track_mrtrix_tsf.count=1;
%     track_mrtrix_tsf.total_count=1;
%     write_mrtrix_tsf(track_mrtrix_tsf,[subdir '/BL/along_tract_metrics/mean_ROI' num2str(combinations(j,1)) '_ROI' num2str(combinations(j,2)) '.tsf']); 
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


track_mrtrix_tck=track_mrtrix_all(1);
track_mrtrix_tsf=track_mrtrix_all(1);
track_mrtrix_tsf2=track_mrtrix_all(1);
counter=0;
for i=1:length(combinations)
        for node=1:total_tracks_new(i)
        track_mrtrix_tck.data{node+counter}=track_mrtrix_all(i).data{node};
        track_mrtrix_tsf.data{node+counter}=i*ones(size(track_mrtrix_all(i).data{node},1),1);
        track_mrtrix_tsf2.data{node+counter}=(1:100)';
        end
       counter=counter+total_tracks_new(i);
end
track_mrtrix_tck.count=num2str(3);
track_mrtrix_tck.total_count=num2str(counter);
write_mrtrix_tracks(track_mrtrix_tck,[p '/along_tract_metrics/all.tck']);

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
                connectome(combinations(j,1),combinations(j,2),scal,:)=mean(sample(sum(total_tracks(1:j)):sum(total_tracks(1:j+1))-1,:),1)';
                connectome_std(combinations(j,1),combinations(j,2),scal,:)=std(sample(sum(total_tracks(1:j)):sum(total_tracks(1:j+1))-1,:),1)';
            end
            
        end
    save([p '/along_tract_metrics/Scalar_connectome.mat'],'connectome')
    save([p '/along_tract_metrics/Scalar_connectome_std.mat'],'connectome_std')
