%% initialize 
basedir='/Volumes/Data2/Aphasia_Project';
folders=dir(basedir);
index_ROI=[1:2:93,96,97,99,113,115,167,182:2:188]; % left gray matter
combinations=nchoosek(index_ROI,2);
scalar_maps={'fa','md','dax','drad','mk','kax','krad','kfa'}; % calculate along tract stats of these metrics 
nb_nodes=100;

%% loop over all subjects 
for subj=1:length(folders)
subdir=[basedir '/' folders(subj).name '/output_redo/output_redo/'];
if exist([subdir 'BL/all_jhu_excl.tck'],'file') % check if BL folder exists 

mkdir([subdir 'BL/along_tract_metrics'])
hdr=spm_vol([subdir 'BL/DKI_roi.nii']);
pixel_size = abs(diag(hdr.mat))';
pixel_size = pixel_size(1:3); 
dim=hdr.dim;

track_mrtrix = read_mrtrix_tracks([subdir 'BL/all_jhu_excl.tck']); % tracks from mrtrix are in world coordinates (voxels) 
load([subdir 'BL/total_tracks_jhu_excl.mat']);
transform=hdr.mat; % transformation matrix used to calculate tranformation from world coordinate space to image space 

ROI=spm_read_vols(hdr); 
gmwm=spm_read_vols(spm_vol([subdir 'BL/gmwmi_T1_P' folders(subj).name '.nii']));
   
total_tracks=zeros(1,length(combinations)+1); 
 % offset total_tracks with 0 
for cnt=1:length(combinations)
    total_tracks(cnt+1)=total_tracks_final(combinations(cnt,1),combinations(cnt,2));
end
total_tracks_new=zeros(length(combinations),1);

parfor j=1:length(combinations) % all combinations of ROIs 

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
save([subdir '/BL/along_tract_metrics/total_tracks_cleaned.mat'],'total_tracks_final_new');


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
write_mrtrix_tracks(track_mrtrix_tck,[subdir '/BL/along_tract_metrics/all.tck']);

track_mrtrix_tsf.count=num2str(1);
track_mrtrix_tsf.total_count=num2str(counter);
write_mrtrix_tsf(track_mrtrix_tsf,[subdir '/BL/along_tract_metrics/all.tsf']);

track_mrtrix_tsf2.count=num2str(1);
track_mrtrix_tsf2.total_count=num2str(counter);
write_mrtrix_tsf(track_mrtrix_tsf2,[subdir '/BL/along_tract_metrics/all_check_interp.tsf']);


%% calculate 4D Connectome 

connectome=zeros(189,189,length(scalar_maps),nb_nodes);
connectome_std=zeros(189,189,length(scalar_maps),nb_nodes);

load([subdir '/BL/along_tract_metrics/total_tracks_cleaned.mat']);
total_tracks=zeros(1,length(combinations)); 

for cnt=1:length(combinations)
        total_tracks(cnt)=total_tracks_final_new(combinations(cnt,1),combinations(cnt,2));
end
total_tracks=[1,total_tracks];    % offset total_tracks with 1    

        for scal=1:length(scalar_maps)
          command=['tcksample ' subdir '/BL/along_tract_metrics/all.tck ' subdir '/BL/sDKIdu_' scalar_maps{scal} '_dki.nii ' subdir '/BL/along_tract_metrics/' scalar_maps{scal} '.txt -force']  ;
          [~,~]=system(command);
          command=['tcksample ' subdir '/BL/along_tract_metrics/all.tck ' subdir '/BL/sDKIdu_' scalar_maps{scal} '_dki.nii ' subdir '/BL/along_tract_metrics/' scalar_maps{scal} '.tsf -force']  ;
          [s,t]=system(command);

          sample=load([subdir 'BL/along_tract_metrics/' scalar_maps{scal} '.txt' ]);
          if scal==2 || scal==3 || scal==4, sample=sample*1000; end % change units of diffusivity measures to um2/ms

            for j=1:length(combinations)
                connectome(combinations(j,1),combinations(j,2),scal,:)=mean(sample(sum(total_tracks(1:j)):sum(total_tracks(1:j+1))-1,:),1)';
                connectome_std(combinations(j,1),combinations(j,2),scal,:)=std(sample(sum(total_tracks(1:j)):sum(total_tracks(1:j+1))-1,:),1)';
            end
            
        end
    save([subdir 'BL/along_tract_metrics/Scalar_connectome.mat'],'connectome')
    save([subdir 'BL/along_tract_metrics/Scalar_connectome_std.mat'],'connectome_std')
end
end

%% 1W

%% loop over all subjects 
for subj=1:length(folders)
subdir=[basedir '/' folders(subj).name '/output_redo/output_redo/'];
if exist([subdir '1W/all_jhu_excl.tck'],'file') % check if BL folder exists 

mkdir([subdir '1W/along_tract_metrics'])
hdr=spm_vol([subdir '1W/DKI_roi.nii']);
pixel_size = abs(diag(hdr.mat))';
pixel_size = pixel_size(1:3); 
dim=hdr.dim;

track_mrtrix = read_mrtrix_tracks([subdir '1W/all_jhu_excl.tck']); % tracks from mrtrix are in world coordinates (voxels) 
load([subdir '1W/total_tracks_jhu_excl.mat']);
transform=hdr.mat; % transformation matrix used to calculate tranformation from world coordinate space to image space 

ROI=spm_read_vols(hdr); 
gmwm=spm_read_vols(spm_vol([subdir '1W/gmwmi_T1_P' folders(subj).name '.nii']));
   
total_tracks=zeros(1,length(combinations)+1); 
 % offset total_tracks with 0 
for cnt=1:length(combinations)
    total_tracks(cnt+1)=total_tracks_final(combinations(cnt,1),combinations(cnt,2));
end
total_tracks_new=zeros(length(combinations),1);

parfor j=1:length(combinations) % all combinations of ROIs 

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
%     write_mrtrix_tracks(track_mrtrix_tck,[subdir '/1W/along_tract_metrics/ROI' num2str(combinations(j,1)) '_ROI' num2str(combinations(j,2)) '.tck']);
%     track_mrtrix_tsf.count=1;
%     track_mrtrix_tsf.total_count=sum(index);
%     write_mrtrix_tsf(track_mrtrix_tsf,[subdir '/1W/along_tract_metrics/ROI' num2str(combinations(j,1)) '_ROI' num2str(combinations(j,2)) '.tsf']);
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
%     write_mrtrix_tracks(track_mrtrix_tck,[subdir '/1W/along_tract_metrics/mean_ROI' num2str(combinations(j,1)) '_ROI' num2str(combinations(j,2)) '.tck']);
%     track_mrtrix_tsf.count=1;
%     track_mrtrix_tsf.total_count=1;
%     write_mrtrix_tsf(track_mrtrix_tsf,[subdir '/1W/along_tract_metrics/mean_ROI' num2str(combinations(j,1)) '_ROI' num2str(combinations(j,2)) '.tsf']); 
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
save([subdir '/1W/along_tract_metrics/total_tracks_cleaned.mat'],'total_tracks_final_new');


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
write_mrtrix_tracks(track_mrtrix_tck,[subdir '/1W/along_tract_metrics/all.tck']);

track_mrtrix_tsf.count=num2str(1);
track_mrtrix_tsf.total_count=num2str(counter);
write_mrtrix_tsf(track_mrtrix_tsf,[subdir '/1W/along_tract_metrics/all.tsf']);

track_mrtrix_tsf2.count=num2str(1);
track_mrtrix_tsf2.total_count=num2str(counter);
write_mrtrix_tsf(track_mrtrix_tsf2,[subdir '/1W/along_tract_metrics/all_check_interp.tsf']);


%% calculate 4D Connectome 

connectome=zeros(189,189,length(scalar_maps),nb_nodes);
connectome_std=zeros(189,189,length(scalar_maps),nb_nodes);

load([subdir '/1W/along_tract_metrics/total_tracks_cleaned.mat']);
total_tracks=zeros(1,length(combinations)); 

for cnt=1:length(combinations)
        total_tracks(cnt)=total_tracks_final_new(combinations(cnt,1),combinations(cnt,2));
end
total_tracks=[1,total_tracks];    % offset total_tracks with 1    

        for scal=1:length(scalar_maps)
          command=['tcksample ' subdir '/1W/along_tract_metrics/all.tck ' subdir '/1W/sDKIdu_' scalar_maps{scal} '_dki.nii ' subdir '/1W/along_tract_metrics/' scalar_maps{scal} '.txt -force']  ;
          [~,~]=system(command);
          command=['tcksample ' subdir '/1W/along_tract_metrics/all.tck ' subdir '/1W/sDKIdu_' scalar_maps{scal} '_dki.nii ' subdir '/1W/along_tract_metrics/' scalar_maps{scal} '.tsf -force']  ;
          [s,t]=system(command);

          sample=load([subdir '1W/along_tract_metrics/' scalar_maps{scal} '.txt' ]);
          if scal==2 || scal==3 || scal==4, sample=sample*1000; end % change units of diffusivity measures to um2/ms

            for j=1:length(combinations)
                connectome(combinations(j,1),combinations(j,2),scal,:)=mean(sample(sum(total_tracks(1:j)):sum(total_tracks(1:j+1))-1,:),1)';
                connectome_std(combinations(j,1),combinations(j,2),scal,:)=std(sample(sum(total_tracks(1:j)):sum(total_tracks(1:j+1))-1,:),1)';
            end
            
        end
    save([subdir '1W/along_tract_metrics/Scalar_connectome.mat'],'connectome')
    save([subdir '1W/along_tract_metrics/Scalar_connectome_std.mat'],'connectome_std')
end
end

% for  j=1:length(combinations)
%     connectome(combinations(j,1),combinations(j,2),:,:)=connectome_int(j,:,:);
% end
       
%         %% see if there is any clusters in bundle ( messes up interpolation)
%         test_index=[1:length(tracks_trk)];
%         test_index_subset=test_index(index); 
%         
%         tracks_trk_matrix=trk_restruc(tracks_trk);
%         tracks_trk_matrix_flip=tracks_trk_matrix;
%         tracks_trk_matrix_flip(:,:,(sum(all==[1,0,0,1],2)==4)>0)=tracks_trk_matrix(fliplr(1:end),:,(sum(all==[1,0,0,1],2)==4)>0); % flip tracks that start in ROI1 
%         tracks_trk_test=trk_restruc(tracks_trk_matrix_flip);
% 
%         int_mean=zeros(sum(index),3);
%         for i=1:sum(index)
%         int=tracks_trk_test(test_index_subset).matrix;
%         int(isnan(int(:,1)),:)=[];
%         int_mean(i,:)=mean(int); % average geometry for each track
%         tracks_trk_test(test_index_subset).matrix=int;
%         end
%         [idx,c]=kmeans(int_mean,2);
%         dist=sqrt((sum(c(1,:)-c(2,:)).^2)); % distance between two centroids 
%         
%        if dist>5 % if centroids are more than 5mm apart split up 
%          Cluster1=(idx==1) ; 
%          interpolated_trk=trk_interp(tracks_trk_test(test_index_subset(Cluster1)),nb_nodes,[],1); % interpolate tracts into 100 nodes, will make the number of nodes uneven because of the tie at center option 
%          
% 
%          
%          Cluster2=(idx==2);  
% 
%            
%        end
%        
       
%     tracks_interp_str = trk_restruc(interpolated_trk_flip);
%     
%        % calculate weights for diffusion metrics based on distance to
%     % geometric average (tracts further away from the mean contribute less
%     % to the average metric)(borrowed strategy from AFQ/mrDiffusion: dtiFIberGroupPropertyWeightedAverage.m & AFQ_FiberTractGaussian.m) 
%     
%     weights=[]; 
%     weights_param=[];
%     track_mean=mean(tracks_interp,3);
%     
%         for node=1:nb_nodes
%         int=tril(cov(permute(tracks_interp(node,:,:),[3 2 1])));
%         cov_matrix=int([1:3 5:6 9]');
%         sigma=[cov_matrix(1:3)'
%                0 cov_matrix(4:5)'
%                0 0 cov_matrix(6)'];
%         sigma=sigma+sigma'-diag(diag(sigma));
%         mu=track_mean(node,:);
%         d=bsxfun(@minus,squeeze(tracks_interp(node,:,:))',mu);
%         % mahalanobis distance of each point on each fiber from the tract core  
%         weights(node,:)=sqrt(dot(d/(sigma),d,2))';
%         weights_param(node,:)=mvnpdf(squeeze(tracks_interp(node,:,:))',mu,sigma)';     
%         end
%         
%     index_keep=(sum(weights<4)==nb_nodes); % remove tracks that are more than 4 standard deviations away from tract mean 
%     weights_param=(weights_param); % parfor loops needs this??? 
%     tracks_interp=tracks_interp(:,:,index_keep);
%     
    
       
%     % calculate weights for diffusion metrics based on distance to
%     % geometric average (tracts further away from the mean contribute less
%     % to the average metric)(borrowed strategy from AFQ/mrDiffusion: dtiFIberGroupPropertyWeightedAverage.m & AFQ_FiberTractGaussian.m) 
%     
%     weights=[]; 
%     weights_param=[];
%     track_mean=mean(tracks_interp,3);
%     
%         for node=1:nb_nodes
%         int=tril(cov(permute(tracks_interp(node,:,:),[3 2 1])));
%         cov_matrix=int([1:3 5:6 9]');
%         sigma=[cov_matrix(1:3)'
%                0 cov_matrix(4:5)'
%                0 0 cov_matrix(6)'];
%         sigma=sigma+sigma'-diag(diag(sigma));
%         mu=track_mean(node,:);
%         d=bsxfun(@minus,squeeze(tracks_interp(node,:,:))',mu);
%         % mahalanobis distance of each point on each fiber from the tract core  
%         weights(node,:)=sqrt(dot(d/(sigma),d,2))';
%         weights_param(node,:)=mvnpdf(squeeze(tracks_interp(node,:,:))',mu,sigma)';     
%         end
%         
%     index_keep=(sum(weights<4)==nb_nodes); % remove tracks that are more than 4 standard deviations away from tract mean 
%     weights_param=(weights_param); % parfor loops needs this??? 
%     tracks_interp=tracks_interp(:,:,index_keep);
%     
%     weightsNormalized = weights_param(:,index_keep)./(repmat(sum(weights_param(:,index_keep), 2), [1 sum(index_keep)])); % normalize weights so the sum of all nodes is 1 (weighted average) 
% 
%     tracks_interp_str = trk_restruc(tracks_interp);
%     header.n_count=length(tracks_interp_str); %recalculate header count because of removed tracks 
%      
%     %mean track 
%     track_mean = mean(tracks_interp, 3);
%     int=(transform*[(track_mean./repmat(pixel_size,length(track_mean),1)) ones(length(track_mean),1)]'); % convert to .tck format (world coordinates in voxels) 
%     track_mean_mrtrix(:,:,i)=int(1:3,:);
%  
%     total_tracks(i)=sum(index_keep); %update total tracks after filtering
% 
%     % save all bundles in .tck
%     
%         for nb=1:sum(index_keep)
%         int=(transform*[(tracks_interp(:,:,nb)./repmat(pixel_size,length(track_mean),1)) ones(length(track_mean),1)]'); % convert to .tck format (world coordinates in voxels) 
%         track_mrtrix_temp{nb}=int(1:3,:)';
%         end
%         
%     track_mrtrix(i,:)=track_mrtrix_temp';
%     
%         % SCALAR MAPS 
%         for sc_map=1:length(scalar_maps) 
%         hdr=spm_vol([p '/' dwi_name '_' scalar_maps{sc_map} '_dki.nii']);
%         map=spm_read_vols(hdr);
%     
%         if sc_map==2, map=map*1000;end % convert MD to more conventional um2/ms 
% 
%         [meanInt_temp(sc_map), stdInt_temp(sc_map), ~] = trk_stats(header, tracks_interp_str, map, 'nearest');
%         [header_sc, tracks_sc] = trk_add_sc(header, tracks_interp_str, map, scalar_maps{sc_map});
%     
%         %initalize
%         scalars = zeros(tracks_sc(1).nPoints, header_sc.n_count, header_sc.n_scalars);
%         scalars_nw = zeros(tracks_sc(1).nPoints, header_sc.n_count, header_sc.n_scalars);
% 
%         fw = weightsNormalized(:);
% 
%             for sc=1:header_sc.n_scalars
%                 mat_long        = cat(1, tracks_sc.matrix);
%                 scalars(:,:,sc)  = reshape(mat_long(:,4).*fw, tracks_sc(1).nPoints, header_sc.n_count, header_sc.n_scalars);
%                 scalars_nw(:,:,sc)  = reshape(mat_long(:,4), tracks_sc(1).nPoints, header_sc.n_count, header_sc.n_scalars);
%             end
% 
%         scalar_mean_temp(:,sc_map) = squeeze(sum(scalars, 2)); % sum because they are weighted scalars 
%         scalar_sd_temp(:,sc_map)  = squeeze(nanstd(scalars_nw, 0, 2));
%         end
%         
%     meanInt(i,:)=meanInt_temp';
%     stdInt(i,:)=stdInt_temp';
%     scalar_mean(i,:,:)=scalar_mean_temp;
%     scalar_sd(i,:,:)=scalar_sd_temp;
