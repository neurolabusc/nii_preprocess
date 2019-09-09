function [X,param] = estimate_parameters(DT,KT,param)
%DT is diffusion tensor (6 x nvox) and KT is kurtosis tensor (15 x nvox);
%Parameters to estimate are supplied by param. Currently supported
%parameters are: {'md','dax','drad','mk','kax','krad','fa','kfa','fa_rgb'};
%
%Each calculation has been vectorized for speed.
%
%Author: Russell Glenn
% Medical University of South Carolina
%Author has provided permission for distribution with nii_preprocess
if ~exist('rotate_tensors', 'file')
    p = fileparts(which(mfilename));
    addpath(fullfile(p, 'DKI_tractography'));
    if ~exist('rotate_tensors', 'file')
        error('Unable to find rotate_tensors');
    end
end


%Initialize output
X = zeros(length(param),size(DT,2)); if ismember('fa_rgb',param); X = [X;zeros(2,size(DT,2))]; end

%compress DT and KT
midx = find(sqrt(sum([DT;KT].^2)));
DT = DT(:,midx);
KT = KT(:,midx);

%Check Names---------------------------------------------------------------
    param = cellfun(@(x)lower(x),param,'UniformOutput',0);
    keep = ismember(param,{'md','dax','drad','mk','kax','krad','fa','kfa'}); %Supported Parameters
    keep_2 = ismember(param,{'fa_rgb'}); %bonus :)
    keep_name = cell(1,sum([keep,keep_2]));
    [keep_name{1:sum(keep)}] = deal(param{keep});
    if sum(keep_2); keep_name{end} = 'fa_rgb'; end
    param = keep_name; clear keep_name keep keep_2

%GET ADDITIONAL INFORMATION THAT IS NECESSARY FOR PARTICULAR PARAMETERS----------------------------------
    if sum(ismember({'fa_rgb'},param)); DTX = DT; end %Copy in case rotated
    if sum(ismember({'dax','drad','kax','krad'},param)); parfor i = 1:size(DT,2); [DT(:,i),KT(:,i)] = rotate_tensors(DT(:,i),KT(:,i)); end; end %Rotate Tensors
    if sum(ismember({'mk','fa_rgb'},param));
        %Sampling Distribution
        [S IDX,~,AREA] = sphericalgrid4; AREA = AREA(IDX(:,1))./sum(AREA);
        R = [sin(S(IDX(:,1),1)).*cos(S(IDX(:,1),2)),sin(S(IDX(:,1),1)).*sin(S(IDX(:,1),2)),cos(S(IDX(:,1),1))];
        %Directional Diffusivity over R
        DDR = [R(:,1).^2,R(:,2).^2,R(:,3).^2,2.*[R(:,1).*R(:,2),R(:,1).*R(:,3),R(:,2).*R(:,3)]];
        if ismember({'mk'},param)
        %Directional kurtosis over R
        DKR = [R(:,1).^4,R(:,2).^4,R(:,3).^4,...
        4.*[R(:,1).^3.*R(:,2),R(:,1).^3.*R(:,3),R(:,1).*R(:,2).^3,...
        R(:,1).*R(:,3).^3,R(:,2).^3.*R(:,3),R(:,2).*R(:,3).^3],...
        6.*[R(:,1).^2.*R(:,2).^2,R(:,1).^2.*R(:,3).^2,R(:,2).^2.*R(:,3).^2],...
        12.*[R(:,1).^2.*R(:,2).*R(:,3),R(:,1).*R(:,2).^2.*R(:,3),R(:,1).*R(:,2).*R(:,3).^2]];
        end;
    end
    if ismember({'kax'},param);
        %Directional Diffusivity and Kurtosis in principal direction
        DDA = [1,0,0,0,0,0];
        DKA = [1,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
    end
    if ismember({'krad'},param)
        nrad = 360;
        rad = linspace(0,2*pi-2*pi/nrad,nrad);
        RAD = [zeros(length(rad),1),cos(rad)',sin(rad)'];
        %Directional Diffusivity and Kurtosis in direction orthogonal to principal direction
        DDRAD = [RAD(:,1).^2,RAD(:,2).^2,RAD(:,3).^2,2.*[RAD(:,1).*RAD(:,2),RAD(:,1).*RAD(:,3),RAD(:,2).*RAD(:,3)]]; %Directional Diffusivity
        DKRAD = [RAD(:,1).^4,RAD(:,2).^4,RAD(:,3).^4,...
        4.*[RAD(:,1).^3.*RAD(:,2),RAD(:,1).^3.*RAD(:,3),RAD(:,1).*RAD(:,2).^3,...
        RAD(:,1).*RAD(:,3).^3,RAD(:,2).^3.*RAD(:,3),RAD(:,2).*RAD(:,3).^3],...
        6.*[RAD(:,1).^2.*RAD(:,2).^2,RAD(:,1).^2.*RAD(:,3).^2,RAD(:,2).^2.*RAD(:,3).^2],...
        12.*[RAD(:,1).^2.*RAD(:,2).*RAD(:,3),RAD(:,1).*RAD(:,2).^2.*RAD(:,3),RAD(:,1).*RAD(:,2).*RAD(:,3).^2]];
    end
%-------------------------------------------------------------------------------------------------------

%GET PARAMETERS:
%NOTE: There are a few places where this requires a good bit of memory. In
%those instances it will try to solve everything simultaneoulsy but in case
%memory is not avaliable, it will do voxelwise
for nidx = 1:length(param)
    ni = param{nidx};
    fprintf('Calculating: %s\n',ni)
    if strcmpi(ni,'md'); X(nidx,midx) = mean(DT(1:3,:));
    elseif strcmpi(ni,'dax'); X(nidx,midx) = DT(1,:);
    elseif strcmpi(ni,'drad'); X(nidx,midx) = mean(DT(2:3,:));
    elseif strcmpi(ni,'fa'); X(nidx,midx) = sqrt(3.*sum([DT(1:3,:)-repmat(sum(DT(1:3,:))./3,3,1);repmat(DT(4:6,:),2,1)].^2)./sum([DT(1:3,:);repmat(DT(4:6,:),2,1)].^2)./2);
    elseif strcmpi(ni,'kfa')
        X(nidx,midx) = sqrt(((KT(1,:)-(KT(1,:)+KT(2,:)+KT(3,:)+2.*KT(10,:)+2.*KT(11,:)+2.*KT(12,:))./5).^2+(KT(2,:)-(KT(1,:)+KT(2,:)+KT(3,:)+2.*KT(10,:)+...
        2.*KT(11,:)+2.*KT(12,:))./5).^2+(KT(3,:)-(KT(1,:)+KT(2,:)+KT(3,:)+2.*KT(10,:)+2.*KT(11,:)+2.*KT(12,:))./5).^2+...
        6.*(KT(10,:)-(KT(1,:)+KT(2,:)+KT(3,:)+2.*KT(10,:)+2.*KT(11,:)+2.*KT(12,:))./15).^2+6.*(KT(11,:)-(KT(1,:)+KT(2,:)+KT(3,:)+2.*KT(10,:)+...
        2.*KT(11,:)+2.*KT(12,:))./15).^2+6.*(KT(12,:)-(KT(1,:)+KT(2,:)+KT(3,:)+2.*KT(10,:)+2.*KT(11,:)+2.*KT(12,:))./15).^2+...
        4.*KT(4,:).^2+4.*KT(5,:).^2+4.*KT(6,:).^2+4.*KT(7,:).^2+4.*KT(8,:).^2+4.*KT(9,:).^2+12.*KT(13,:).^2+12.*KT(14,:).^2+12.*KT(15,:).^2)./(KT(1,:).^2+...
        KT(2,:).^2+KT(3,:).^2+6.*KT(10,:).^2+6.*KT(11,:).^2+6.*KT(12,:).^2+4.*KT(4,:).^2+4.*KT(5,:).^2+4.*KT(6,:).^2+4.*KT(7,:).^2+4.*KT(8,:).^2+...
        4.*KT(9,:).^2+12.*KT(13,:).^2+12.*KT(14,:).^2+12.*KT(15,:).^2));
    elseif strcmpi(ni,'mk');
        try %all at once
            mk = sum(repmat(AREA,1,length(midx)).*(repmat(mean(DT(1:3,:)),size(DKR,1),1).^2).*(DKR*KT)./((DDR*DT).^2)); mk(mk<0)=0; mk(mk>3)=3;
            X(nidx,midx) = mk;
        catch %voxelwise
            parfor i = 1:size(X,2); X(nidx,i) = sum(AREA.*repmat(mean(DT(1:3,i)).^2,size(DKR,1),1).*(DKR*KT(:,i))./((DDR*DT(:,i)).^2)); end
        end
    elseif strcmpi(ni,'kax'); X(nidx,midx) = repmat(mean(DT(1:3,:)).^2,size(DKA,1),1).*(DKA*KT)./((DDA*DT.^2));
    elseif strcmpi(ni,'krad');
        try %all at once
            X(nidx,midx) = mean(repmat(mean(DT(1:3,:)).^2,size(DKRAD,1),1).*(DKRAD*KT)./((DDRAD*DT.^2)));
        catch %voxelwise
            parfor i = 1:size(X,2); X(nidx,i) = mean(repmat(mean(DT(1:3,i)).^2,size(DKRAD,1),1).*(DKRAD*KT(:,i))./((DDRAD*DT(:,i)).^2)); end
        end
    elseif strcmpi(ni,'fa_rgb')
        try %all at once
            dd_n = DDR*DTX;
            fa = sqrt(3.*sum([DTX(1:3,:)-repmat(sum(DTX(1:3,:))./3,3,1);repmat(DTX(4:6,:),2,1)].^2)./sum([DTX(1:3,:);repmat(DTX(4:6,:),2,1)].^2)./2);
            [v1,~] = ind2sub(size(dd_n),find(repmat(max(dd_n),size(dd_n,1),1)==dd_n));
            X(nidx:nidx+2,midx) = (repmat(fa',1,3).*abs(R(v1,:)))';
        catch %voxel-wise
            x1 = zeros(1,size(DTX,2)); x2 = zeros(1,size(DTX,2)); x3 = zeros(1,size(DTX,2));
            parfor i = 1:size(DTX,2);
                fa = sqrt(3.*sum([DTX(1:3,i)-repmat(sum(DTX(1:3,i))./3,3,1);repmat(DTX(4:6,i),2,1)].^2)./sum([DTX(1:3,i);repmat(DTX(4:6,i),2,1)].^2)./2);
                [L V] = eig([DTX(1,i),DTX(4,i),DTX(5,i);DTX(4,i),DTX(2,i),DTX(6,i);DTX(5,i),DTX(6,i),DTX(3,i)]);
                [~,vi] = sort(diag(V),'descend');
                x1(i) = fa*L(1,vi(1)); x2(i) = fa*L(2,vi(1)); x3(i) = fa*L(3,vi(1));
            end
            X(nidx:nidx+2,midx) = abs([x1;x2;x3]);
        end
    end
end
