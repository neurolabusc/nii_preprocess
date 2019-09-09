function kODF_nii_preprocess(FT_parameters,dki_dt,dki_kt,dki_fa)
%dke_ft(FT_PARAMETERS) optimizes the kurtosis dODF and performs
%DKI-based white matter fiber tractography
%
%See DKE FT MODULE USER'S GUIDE for more info
%
%TECHNICAL PUBLICATIONS:
%     1.	Jensen et al. Leading non-Gaussian corrections for diffusion orientation diffusion orientation distribution function. NMR Biomed. 2014;27:202-11. 
%     2.	Glenn et al. Optimization of white matter fiber tractography with diffusional kurtosis imaging. NMR Biomed. [In Press]. 
%
%Author: Russell Glenn
%Medical University of South Carolina
%July 2015

%TEMP-------------------
%Come back and add to params file
kfa_colormap = 1; 
%-----------------------
FTVersion = GetFTVersion;

%Get input parameters-------------------------------------------------------
input_struc = readvariables(FT_parameters);

%Create output structure---------------------------------------------------
fn_out_struc.fib='dki';
fn_out_struc.gfa='gfa';
fn_out_struc.nfd='nfd';
fn_out_struc.odf_k_min='odf_k_min';
fn_out_struc.odf_k='odf_k';
fn_out_struc.odf_k_max='odf_k_max';
fn_out_struc.odf_d='odf_d';
fn_out_struc.dki_odf_coeff='odf_coeff';
fn_out_struc.gfa_rgb='gfa_rgb';
fn_out_struc.fa_rgb='fa_rgb';
fn_out_struc.kfa_rgb = 'kfa_rgb'; 
fn_out_struc.FT_struct='FT_struct';
fn_out_struc.FT_dki='DKI';
fn_out_struc.FT_dti='DTI';
fn_out_struc.SH='SH_coeff';

fnames = fieldnames(fn_out_struc);

% %Check for parallel computing
% try if matlabpool('size')==0; matlabpool open; fprintf('\n'); end; 
% catch ME
%     ME.message;
% end 

warning('off','all')

%--------------------------------------------------------------------------
%GET DATA
%--------------------------------------------------------------------------

%Initialize data that is constant for all subjects
    if input_struc.odf_optimization&&input_struc.quasiNewton; 
        BFGS_options = optimset('largescale','off','GradObj','off','Display','off','TolX',1E-9);    %BFGS
    end
    
    %Get necessary permutations / inversions to convert image volume and gradient table to LPS
    input_struc.image_orientation = lower(input_struc.image_orientation); 
    input_struc.odf_orientation = lower(input_struc.odf_orientation);

    permute_img = [[strfind(input_struc.image_orientation,'l') strfind(input_struc.image_orientation,'r')],[strfind(input_struc.image_orientation,'a') strfind(input_struc.image_orientation,'p')],[strfind(input_struc.image_orientation,'s') strfind(input_struc.image_orientation,'i')]]; 
    input_struc.image_orientation = input_struc.image_orientation(permute_img); 
    invert_img = 2*[strcmpi(input_struc.image_orientation(1),'l'), strcmpi(input_struc.image_orientation(2),'p'), strcmpi(input_struc.image_orientation(3),'s')]-1; 

    permute_odf = [[strfind(input_struc.odf_orientation,'l') strfind(input_struc.odf_orientation,'r')],[strfind(input_struc.odf_orientation,'a') strfind(input_struc.odf_orientation,'p')],[strfind(input_struc.odf_orientation,'s') strfind(input_struc.odf_orientation,'i')]]; 
    input_struc.image_orientation = input_struc.image_orientation(permute_odf); 
    invert_odf = 2*[strcmpi(input_struc.odf_orientation(1),'l'), strcmpi(input_struc.odf_orientation(2),'p'), strcmpi(input_struc.odf_orientation(3),'s')]-1; 
    
    if input_struc.odf_optimization
        %Get Sampling distribution and other info
        [S, IDX, idx8, AREA, ~, separation_angle] = feval(str2func(sprintf( ['sphericalgrid' num2str(input_struc.sd)])));
        AREA = AREA./sum(AREA);    %For GFA calculation    

        %Go ahead and get kurtosis dODF functions: These don't need to keep being defined
        [~,fK,fK2] = getODF_FCN(zeros(6,1),zeros(15,1),4); 
        R = [sin(S(:,1)).*cos(S(:,2)), sin(S(:,1)).*sin(S(:,2)),cos(S(:,1))]; %Sampling distribution in Cartesian coordinates

        if input_struc.make_fib_file
            if strcmpi(input_struc.odf_res,'low')||strcmpi(input_struc.odf_res,'high')&&input_struc.sd==3
                idx_fib = idx8;
                [SX IDXX , ~, ~, odf_faces] = sphericalgrid3;

            elseif strcmpi(input_struc.odf_res,'high')&&input_struc.sd>=4
                idx_fib = 1:1281; 
                [SX IDXX , ~, ~, odf_faces] = sphericalgrid4;
            else
                error('Could not find a correct odf_res value. Please input ''low'' or ''high'' for the odf_res input'); 
            end

            V = [sin(SX(IDXX(:,1),1)).*cos(SX(IDXX(:,1),2)), sin(SX(IDXX(:,1),1)).*sin(SX(IDXX(:,1),2)), cos(SX(IDXX(:,1),1))]';
            V = bsxfun(@times,V,invert_odf'); 
            odf_vertices = [V -V];
            clear V SX IDXX
        end
        
        if kfa_colormap
            %Define some functions to use
            
            %Directional dependence of the kurtosis tensor (n = m x 3 Cartesian coordinates)
            W_fcn = @(n,kt)n(:,1).^4.*kt(1)+n(:,2).^4.*kt(2)+n(:,3).^4.*kt(3)+4.*n(:,1).^3.*n(:,2).*kt(4)+4.*n(:,1).^3.*n(:,3).*kt(5)+4.*n(:,1).*n(:,2).^3.*kt(6)+...
                4.*n(:,1).*n(:,3).^3.*kt(7)+4.*n(:,2).^3.*n(:,3).*kt(8)+4.*n(:,2).*n(:,3).^3.*kt(9)+6.*n(:,1).^2.*n(:,2).^2.*kt(10)+6.*n(:,1).^2.*n(:,3).^2.*kt(11)+...
                6.*n(:,2).^2.*n(:,3).^2.*kt(12)+12.*n(:,1).^2.*n(:,2).*n(:,3).*kt(13)+12.*n(:,1).*n(:,2).^2.*n(:,3).*kt(14)+12.*n(:,1).*n(:,2).*n(:,3).^2.*kt(15); 
            
            %KFA function (equivalent to DKE Calculation ~ takes no time to redo)
            kfa_fcn  = @(kt)sqrt(((kt(1,:)-(kt(1,:)+kt(2,:)+kt(3,:)+2.*kt(10,:)+2.*kt(11,:)+2.*kt(12,:))./5).^2+(kt(2,:)-(kt(1,:)+kt(2,:)+kt(3,:)+2.*kt(10,:)+...
                2.*kt(11,:)+2.*kt(12,:))./5).^2+(kt(3,:)-(kt(1,:)+kt(2,:)+kt(3,:)+2.*kt(10,:)+2.*kt(11,:)+2.*kt(12,:))./5).^2+...
                6.*(kt(10,:)-(kt(1,:)+kt(2,:)+kt(3,:)+2.*kt(10,:)+2.*kt(11,:)+2.*kt(12,:))./5./3).^2+6.*(kt(11,:)-(kt(1,:)+kt(2,:)+kt(3,:)+2.*kt(10,:)+...
                2.*kt(11,:)+2.*kt(12,:))./5./3).^2+6.*(kt(12,:)-(kt(1,:)+kt(2,:)+kt(3,:)+2.*kt(10,:)+2.*kt(11,:)+2.*kt(12,:))./5./3).^2+...
                4.*kt(4,:).^2+4.*kt(5,:).^2+4.*kt(6,:).^2+4.*kt(7,:).^2+4.*kt(8,:).^2+4.*kt(9,:).^2+12.*kt(13,:).^2+12.*kt(14,:).^2+12.*kt(15,:).^2)./(kt(1,:).^2+...
                kt(2,:).^2+kt(3,:).^2+6.*kt(10,:).^2+6.*kt(11,:).^2+6.*kt(12,:).^2+4.*kt(4,:).^2+4.*kt(5,:).^2+4.*kt(6,:).^2+4.*kt(7,:).^2+4.*kt(8,:).^2+...
                4.*kt(9,:).^2+12.*kt(13,:).^2+12.*kt(14,:).^2+12.*kt(15,:).^2));
                
        end
    end
    
    
%Iterate through each subject
for isubject = 1:length(input_struc.subject_list)
    
    subj_dir = fullfile(input_struc.studydir, input_struc.subject_list{isubject});
    cd(subj_dir);
    
    subj_name = ''; 
    if ~isempty(input_struc.subject_list{isubject}); subj_name = [input_struc.subject_list{isubject} '_']; end

    diary off
    fn_diary = fullfile(subj_dir, 'FT.log');
    fid = fopen(fn_diary, 'w');
    if fid < 0
        error('Cannot open output file %s! Output directory does not exist or is write-protected.', fn_diary);
    end
    fclose(fid);
    
    diary(fn_diary)

    fprintf('Start date and time: %s\n', datestr(now, 'mmmm dd, yyyy HH:MM:SS'))
    fprintf('%s\n',FTVersion)% EM
    
    %Get output filenames
    fn_out_struc_subject = cell2struct(cellfun(@(x)[subj_name input_struc.pre_name fn_out_struc.(x) input_struc.post_name],fnames,'uniformoutput',0),fnames);
    
    fpb = fprintf('Optimizing kurtosis dODF...');   
    hdr=spm_vol(dki_fa); % header template   
        
    fa = spm_read_vols(hdr); 
    dimension = hdr.dim; 
    voxel_size = sqrt(sum(hdr.mat(1:3,1:3).^2)); 
    nvox = prod(dimension);                %number of voxels
    
    %GET Tensors and mask for ODF optimization
    if input_struc.odf_optimization
        DT=spm_read_vols(spm_vol(dki_dt));    %Diffusion tensor (DT)
        KT=spm_read_vols(spm_vol(dki_kt));     %Kurtosis tensor (KT)
        KT=reshape(KT,[prod(hdr.dim),15])';
        DT=reshape(DT,[prod(hdr.dim),6])';
 
        %Check dimensions: This can be off, for example, if the output fa
        %image was interpolated
        if nvox~=size(DT,2); 
           error(sprintf('The number of voxels in the FA image and the tensors must match! \n\nMake sure map_interpolation_method.flag = 0 for tensor fitting or reconstruct the fa image from the DT.mat using the Nifti header from the diffusion weighted images')); 
        end
        
        %Reorder mask_idx_input to go through mask index in LPS orientation (for
        %odfn variables in .fib file)
        mask_idx_lps = permute(reshape(1:prod(dimension),dimension).*double(fa~=0),permute_img); 
        for i = find(invert_img==-1); mask_idx_lps = flipdim(mask_idx_lps,i); end
        mask_idx_lps = mask_idx_lps(mask_idx_lps>0);
        
        if kfa_colormap
            kfa = kfa_fcn(KT); kfa(isnan(kfa))=0; 
        end

    end

    %Make Tractography Input structure; 
    if input_struc.tractography_flg
        
        %Vox_to_RAS transormation: TRK vox is in LPS
        vox_to_ras = hdr.mat([permute_img 4],:)*[diag(invert_img) double(invert_img==-1)'.*hdr.dim(permute_img)'-0.5*[1;1;1];0 0 0 1];
        
        FT_struc = struct('fa_threshold',input_struc.fa_threshold,'angle_threshold',input_struc.angle_threshold,...
            'trk_length',input_struc.trk_length,'step_size',input_struc.step_size,'trk_mask',input_struc.trk_mask,'seed_mask',input_struc.seed_mask,...
            'shift',input_struc.shift,'name','','seednum',input_struc.seednum,'permute_odf',permute_odf,'invert_odf',invert_odf,'permute_img',permute_img,...
            'invert_img',invert_img,'hdr',hdr,'vox_to_ras',vox_to_ras,'SEED',[],'outdir',subj_dir,'reset_memory',input_struc.release_memory,...
            'output_DTI_trks',input_struc.output_DTI_trks,'pre_name',[subj_name input_struc.pre_name],'post_name',input_struc.post_name); 
        
        FT_struc.name = [FT_struc.pre_name 'FT_' fn_out_struc.FT_dki FT_struc.post_name '.trk'];
        
        FT_struc.SEED = random_seed_FT(FT_struc); 
                
        if FT_struc.step_size==0; FT_struc.step_size = mean(voxel_size)/2; end
    end
    
%OPTIMIZE ODF--------------------------------------------------------------
if input_struc.odf_optimization
     
        %Initialize .fib file
        if input_struc.make_fib_file == 1         
            save([fn_out_struc_subject.fib '.fib'], 'dimension', 'voxel_size', 'odf_vertices', 'odf_faces', '-v4')     
        end

    %Initialize output parameters
        odf_k = cell(1,nvox);
        odf_k_max = cell(1,nvox); 
        odf_k_min = zeros(1,nvox); 
        odf_d = cell(1,nvox);

        gfa = zeros(1,nvox);
        nfd = zeros(1,nvox);
        odf_coeff = zeros(29,nvox);
        gfa_rgb = zeros(3,nvox);
        fa_rgb = zeros(3,nvox);
        
        Harm_id={1,5:9,17:25,37:49,65:81}; % max degree = 8 
        B_dki=getSH(input_struc.degree,[ atan2(R(:,2),R(:,1)) acos(R(:,3))],'real');%EM
        B_dki=B_dki(:,cell2mat(Harm_id(1:input_struc.degree/2+1)));%EM
       % B_dki(:,[3,5,8,10,12,14,17,19,21,23,25,27])=B_dki(:,[3,5,8,10,12,14,17,19,21,23,25,27])*-1;
       % still trying to figure this one out
        SH_coeff=zeros(length(cell2mat(Harm_id(1:input_struc.degree/2+1))),nvox); %EM
        
        if kfa_colormap
            kfa_rgb = zeros(3,nvox); 
        end

    %--------------------------------------------------------------------------
    %PROCESS KURTOSIS DODF
    %--------------------------------------------------------------------------
    %Note: data is handled in discrete volumes based on odf_size to build the 
    %.fib file since each odfn variable is saved as a block
    tic
    for n = 0:ceil(numel(mask_idx_lps)/input_struc.odf_size)-1 %get data for odfn
    % for n = 2
        range = n*input_struc.odf_size+1:min((n+1)*input_struc.odf_size,numel(mask_idx_lps));

        idxn = mask_idx_lps(range); %indices in the actual 3D image volume (needed to store things in right place)
        DTN = DT(:,idxn);     %reslice for parloop ~ needed for parfor loop per matlab restrictions
        KTN = KT(:,idxn);     %reslice for parloop ~ needed for parfor loop per matlab restrictions

        if input_struc.make_fib_file
            odfs = zeros(length(idx_fib),numel(idxn)); %initialize odf values to save
        end

        %Re-Initialize variables for this iteration per  matlab parfor loop
        %restrictions on indexing
        
        odf_kn = cell(1,numel(idxn));
        odf_k_maxn = cell(1,numel(idxn)); 
        odf_k_minn = zeros(1,numel(idxn)); 
        odf_dn = cell(1,numel(idxn));
        gfa_n = zeros(1,numel(idxn));
        fa_n = fa(idxn); 
        nfd_n = zeros(1,numel(idxn));

        odf_coeff_n = zeros(29,numel(idxn));
        gfa_rgb_n = zeros(3,numel(idxn));
        fa_rgb_n = zeros(3,numel(idxn));
        
        A2l=zeros(length(cell2mat(Harm_id(1:input_struc.degree/2+1))),numel(idxn));
        
        if kfa_colormap
            kfa_n = kfa(idxn); 
            kfa_rgb_n = zeros(3,numel(idxn)); 
        end

        fprintf(repmat('\b',1,fpb));
        fpb = fprintf('Optimizing kurtosis dODF ~ %s odf%d @ %0.2f min...',input_struc.subject_list{isubject},n,toc/60);

        %Iterate through all values in this volumes
        %par
        for i = 1:length(idxn)
            try
                dt = DTN(:,i); kt = KTN(:,i);           %Get Tensors 
                if fa_n(i)>0.90; %Very high fa dODFs wit particularly low 
                    %eigenvalues may 'explode' due to U = DavgD^-1;
                    %This could potentially be fixed during tensor estimation
                    [dt kt] = regularize_tensors(dt,kt,0.9); 
                end
                
                %EVALUATE KURTOSIS dODF------------------------------------
                A = getODF_FCN(dt,kt,input_struc.radial_weight);       %Get dODF coefficients              
                DKI_P = -fK2(R,A);  %Evaluate dODF over spherical grid (invert bc fK is inverted for minimization)
                %----------------------------------------------------------


                %Find local maxima of dODF over spherical grid (brute force)
                max_idx = IDX(DKI_P(IDX(:,1))==max(DKI_P(IDX),[],2),1); 
 
                if input_struc.quasiNewton==1 %Refine peak estimates with non-linear optimization (quasi-Newton Method)
                    DKI_V = S(max_idx,:); %Switch back to spherical coordinates for non-linear
                    DKI_max = DKI_P(max_idx); 

                     for j = 1:size(DKI_V,1);  
%                         [p,fval,ef] = fminunc(fK,DKI_V(j,:),BFGS_options);
                        [p,fval,ef] = fminunc(@(x)fK(x,A),DKI_V(j,:),BFGS_options);
                        if ef==1||ef==2||ef==3||ef==5
                            DKI_max(j) = -fval;         %Invert because we found minimum (cf getODF_FCN) 
                            DKI_V(j,:) = p ; %Refined seed point
                        end                    
                     end
                    
                    %Convert to Cartesian
                    DKI_V = [sin(DKI_V(:,1)).*cos(DKI_V(:,2)), sin(DKI_V(:,1)).*sin(DKI_V(:,2)), cos(DKI_V(:,1))]'; 
                
                else
                    DKI_V = R(max_idx,:)'; 
                    DKI_max = DKI_P(max_idx);                
                end

                %Sort by magnitude
                [DKI_max idx] = sort(DKI_max,'descend');
                DKI_V = DKI_V(:,idx);

                 % get SH coefficients for mrtrix % EM 09/12/2018 
        
                A2l(:,i)=(B_dki'*B_dki)^-1*B_dki'*DKI_P/DKI_max(1) *fa_n(i); % scale the odfs based on fa so I can control stopping criteria in mrtrix 
                
                %CHECK FOR DUPLICATES: This effects the number of peaks detected.
                chk_i = 1; 
                while chk_i < size(DKI_V,2)

                    %Get vertices within a certain angle threshold
                    deg = real(acosd(DKI_V(:,chk_i)'*DKI_V)); 
                    degi = deg<separation_angle*2; 

                    %Average directions that fall within that threshold
                    DKI_V(:,chk_i) = mean(DKI_V(:,degi),2); 
                    DKI_V(:,chk_i) = DKI_V(:,chk_i)./sqrt(sum(DKI_V(:,chk_i).^2)); 
                    DKI_max(chk_i) = mean(DKI_max(degi));

                    %Remove any duplicates
                    degi(chk_i)=0; 
                    DKI_V = DKI_V(:,~degi);       
                    DKI_max = DKI_max(~degi); 
                    chk_i = chk_i+1; 
                end      


                %Get DTI directions
                [DTI_V L] = eig([dt(1) dt(4) dt(5); dt(4) dt(2) dt(6); dt(5) dt(6) dt(3)]); 
                [L, Li] = sort(diag(L),'descend');

                gfa_i = sqrt(1-sum(AREA.*DKI_P).^2/sum(AREA.*DKI_P.^2));


                %SAVE THINGS
                odf_kn{i} = DKI_V; 
                odf_k_maxn{i} = DKI_max; 
                odf_k_minn(i) = min(DKI_P);
                odf_dn{i} = DTI_V(:,Li(1)); 

                gfa_n(i) = gfa_i; 
                nfd_n(i) = size(DKI_V,2);

                odf_coeff_n(:,i) = A; 
               
                %Permutation used so L-R is red, A-P is green, and I-S is blue
                gfa_rgb_n(:,i) = gfa_i*abs(DKI_V(permute_odf,1));
                fa_rgb_n(:,i) = fa_n(i)*abs(DTI_V(permute_odf,Li(1)));
                
                if kfa_colormap
                    W = W_fcn(R(IDX(:,1),:),kt); 
                    W_V = R(find(W==max(W),1),:)'; 
                    kfa_rgb_n(:,i) = kfa_n(i).*abs(W_V(permute_odf));
                end
                    

                if input_struc.make_fib_file && input_struc.save_odfs
                    odfs(:,i) = DKI_P(idx_fib); 
                end
              
            catch  ME
                ME.message;
            end
        end
                
        %There is a potential memory leak on matlabs parallel computing
        %toolbox. This quick fix releases memory accumulated during the parfor
        %loop. 
        if input_struc.release_memory==2; try matlabpool close; matlabpool open; catch ;end; end %#ok

        if input_struc.make_fib_file
            %NOTE: In case images volumes were co-registered, some voxels  (particularly from the
            %first or last slice) may cause errors with the dODF calculation. For
            %visualization with DSI Studio these need to be removed so they can
            %accuratly match up the fa values (stored as a 1xn array) with the
            %proper odfs (stored as a 12818xm array). This occurs in only a
            %small number of voxels in boundary regions of the image volume
            %following co-registration.
            %
            %This is double and tripple checking compatability between odfn
            %variables and fa, which DSI studio indexes differently.
            %ODF-derived parameters are also included as these voxels have some
            %funny stuff going on. This does not re-write or modify the fa.nii
            %file, but may remove some fa values stored in the .fib file. 

            x = fa_n==0|sum(odfs.^2)'==0; %funny stuff going on in brain mask        
            fa_n(x)=0; gfa_n(x) = 0; nfd_n(x)=0; odf_k_minn(x) = 0;
            [odf_kn{x}] = deal([]); [odf_k_maxn{x}] = deal([]); [odf_dn{x}] = deal([]);
            odfs = odfs(:,~x);

            odfs = input_struc.scale_odf.*bsxfun(@rdivide,odfs,max(odfs)); 

            eval(sprintf('odf%d=odfs;',n)) %Update odfn to new values
            save([fn_out_struc_subject.fib '.fib'],sprintf('odf%d',n),'-v4','-append');  %add odfn to .fib file
            eval(sprintf('clear odf%d odfs;',n)) %clear up some memory for next time        
        end


        %Update output parameters based on values calculated for this block
        [odf_k{idxn}] = deal(odf_kn{:});
        [odf_k_max{idxn}] = deal(odf_k_maxn{:});
        odf_k_min(idxn) = odf_k_minn; 
        [odf_d{idxn}] = deal(odf_dn{:});

        fa(idxn) = fa_n; 
        gfa(idxn) = gfa_n;
        nfd(idxn) = nfd_n; 

        odf_coeff(:,idxn) = deal(odf_coeff_n(:,1:end));
        gfa_rgb(:,idxn) = deal(gfa_rgb_n(:,1:end)); 
        fa_rgb(:,idxn) = deal(fa_rgb_n(:,1:end)); 
        
        SH_coeff(:,idxn)=A2l; %EM
        
        if kfa_colormap
            kfa_rgb(:,idxn) = deal(kfa_rgb_n(:,1:end)); 
        end

    end

    fprintf(repmat('\b',1,fpb));
    fpb = fprintf('\nOptimization of kurtosis dODF complete @ %0.2f min.\n Saving data...',toc/60);

    %--------------------------------------------------------------------------
    %WRITE OUTPUTS FROM ODF OPTIMIZATION 
    %--------------------------------------------------------------------------

    if input_struc.make_fib_file ==1 

        fprintf(repmat('\b',1,fpb));
        fpb = fprintf('Finalizing .fib file...%0.2f min',toc/60);

        %CALCULATE DIR AND FAN PARAMETERS FOR NEW .FIB FILE
        %NOTE: QA values used by DSI studio are replaced by FA values
        %calculated for DKI
        fa_fib = permute(fa,permute_img); 
        C = permute(reshape(odf_k,dimension),permute_img); 
        inv_dim = find(invert_img==-1); 
        for i = inv_dim; fa_fib = flipdim(fa_fib,i); C = flipdim(C,i); end
        idx_fa_fib = find(fa_fib>0); 
        fa_fib = fa_fib(idx_fa_fib); 
        A = cell(size(idx_fa_fib)); 
        [A{:}] = deal(C{idx_fa_fib}); 
        A = cellfun(@(x)bsxfun(@times,x(permute_odf,:),invert_odf'),A,'UniformOutput',0); 

        fax = zeros(5,numel(idx_fa_fib));
        dirx = zeros(3,5,numel(idx_fa_fib));

        parfor i = 1:length(idx_fa_fib);  
            v = A{i};     
            diri = zeros(3,5);    
            xidx = 1:(min(5,size(v,2))); 
            diri(:,xidx)=v(:,xidx);
            dirx(:,:,i) = diri; 
            z = [0 0 0 0 0]'; 
            z(xidx)=1; 
            fax(:,i)=fa_fib(i).*z; 
        end


        fa0 = zeros(1,prod(dimension)); fa1 = fa0; fa2 = fa0; fa3 = fa0; fa4 = fa0; 
        dir0 = zeros(3,prod(dimension)); dir1 = dir0; dir2 = dir0; dir3 = dir0; dir4 = dir0; 

        fa0(idx_fa_fib) = fax(1,:);
        fa1(idx_fa_fib) = fax(2,:);
        fa2(idx_fa_fib) = fax(3,:);
        fa3(idx_fa_fib) = fax(4,:);
        fa4(idx_fa_fib) = fax(5,:);

        dir0(:,idx_fa_fib) = permute(dirx(:,1,:),[1 3 2]); 
        dir1(:,idx_fa_fib) = permute(dirx(:,2,:),[1 3 2]);
        dir2(:,idx_fa_fib) = permute(dirx(:,3,:),[1 3 2]);
        dir3(:,idx_fa_fib) = permute(dirx(:,4,:),[1 3 2]);
        dir4(:,idx_fa_fib) = permute(dirx(:,5,:),[1 3 2]);

        save( [fn_out_struc_subject.fib '.fib'], 'dir0', 'dir1', 'dir2', 'dir3', 'dir4', 'fa0', 'fa1', 'fa2' ,'fa3', 'fa4', '-v4', '-append')
        clear C A fa_fib fa0 fa1 fa2 fa3 fa4 dir0 dir1 dir2 dir3 dir4 fax dirx

    end

    if input_struc.wrt_flg    
        %Build image volumes from vectors
        odf_k = reshape(odf_k,dimension); 
        odf_k_max = reshape(odf_k_max,dimension);
        odf_k_min = reshape(odf_k_min,dimension); 

        odf_d = reshape(odf_d,dimension);

        gfa = reshape(gfa,dimension);
        nfd = reshape(nfd,dimension);

        SH_coeff = reshape(SH_coeff',[dimension length(cell2mat(Harm_id(1:input_struc.degree/2+1)))]);
        
        gfa_rgb = cat(3,reshape(gfa_rgb(1,:),[dimension(1:2) 1 dimension(3)]),reshape(gfa_rgb(2,:),[dimension(1:2) 1 dimension(3)]),reshape(gfa_rgb(3,:),[dimension(1:2) 1 dimension(3)]));
        fa_rgb = cat(3,reshape(fa_rgb(1,:),[dimension(1:2) 1 dimension(3)]),reshape(fa_rgb(2,:),[dimension(1:2) 1 dimension(3)]),reshape(fa_rgb(3,:),[dimension(1:2) 1 dimension(3)]));

        if kfa_colormap
            kfa_rgb = cat(3,reshape(kfa_rgb(1,:),[dimension(1:2) 1 dimension(3)]),reshape(kfa_rgb(2,:),[dimension(1:2) 1 dimension(3)]),reshape(kfa_rgb(3,:),[dimension(1:2) 1 dimension(3)]));
        end
        
        %Write output files
        hdr.dt = [16 0];
        hdr.fname = [fn_out_struc_subject.gfa '.nii']; spm_write_vol(hdr,gfa);    
        hdr.fname = [fn_out_struc_subject.nfd '.nii']; spm_write_vol(hdr,nfd);
        hdr.fname = [fn_out_struc_subject.odf_k_min '.nii']; spm_write_vol(hdr,odf_k_min); 
       
        hdr.fname = [fn_out_struc_subject.SH '.nii'];hdr_4D=repmat(hdr,[1,length(cell2mat(Harm_id(1:input_struc.degree/2+1)))]);
        
        % DO NOT USE 
%         count=0;
%         for l=0:2:input_struc.degree
%             
%             for m=2:2:2*l+1
%              SH_coeff(:,:,:,m+count) = SH_coeff(:,:,:,m+(2*(l-2)))*-1; % mrtrix does not have (-1)^m in formulas! 
%             end
%          
%             count=count+2*l+1;
%             
%         end
        
        for ii=1:(length(cell2mat(Harm_id(1:input_struc.degree/2+1))))
         hdr_4D(ii).n=[ii 1];
         spm_write_vol(hdr_4D(ii),SH_coeff(:,:,:,ii));
        end
        save([fn_out_struc_subject.odf_k '.mat'],'odf_k');
        save([fn_out_struc_subject.odf_k_max '.mat'],'odf_k_max');
        save([fn_out_struc_subject.odf_k_min '.mat'],'odf_k_min');
        save([fn_out_struc_subject.odf_d '.mat'],'odf_d');
        save([fn_out_struc_subject.dki_odf_coeff '.mat'],'odf_coeff');
        save([fn_out_struc_subject.gfa_rgb '.mat'],'gfa_rgb');
        save([fn_out_struc_subject.fa_rgb '.mat'],'fa_rgb');
        
        if kfa_colormap
            save([fn_out_struc_subject.kfa_rgb '.mat'],'kfa_rgb');
        end
    end

end %FINISH PROCESSING ODF-------------------------------------------------

%Perform Tractography
if input_struc.tractography_flg
    if ~input_struc.odf_optimization
        load(fullfile(subj_dir,[fn_out_struc_subject.odf_k '.mat']))
        if FT_struc.output_DTI_trks
            load(fullfile(subj_dir,[fn_out_struc_subject.odf_d '.mat']))
        end
    end
   
    save(fn_out_struc_subject.FT_struct,'-struct','FT_struc');
    fprintf(repmat('\b',1,fpb));
    
    odf_k = reshape(odf_k,dimension);       %Make sure odf's are the right size (in case write flag is off)
   
    EULER_DKE(FT_struc,odf_k, fn_out_struc.FT_dki);
    
    if FT_struc.output_DTI_trks
        odf_d = reshape(odf_d,dimension);   %Make sure odf's are the right size (in case write flag is off)
        EULER_DKE(FT_struc,odf_d,fn_out_struc.FT_dti);
    end

    fpb = 0; 
end


fprintf(repmat('\b',1,fpb))
if input_struc.release_memory; try matlabpool close; matlabpool open; fprintf('\n'); catch ;end; end %#ok
fprintf('\nOptimization Complete...\n\n')
diary off

%Clear things for this subject
clear DT DTN KT KTN subj_dir subj_name fn_diary fid fn_out_struc_subject ...
    dimension voxel_size nvox mask_idx_lps inv_dim FT_struc range idxn odfs ...
    odf_k odf_k_max odf_k_min odf_d gfa nfd odf_coeff fa_rgb fa ...
    odf_kn odf_k_maxn odf_k_minn odf_dn fa_n nfd_n odf_coeff_n gfa_rbg_n ...
    A DKI_P max_idx DKI_V DKI_max idx chi_i degi DTI_V L Li x

end

