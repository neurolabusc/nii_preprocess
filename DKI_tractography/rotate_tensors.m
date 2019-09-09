function [dt_rot kt_rot varargout] = rotate_tensors(dt,kt)
%Rotate diffusion and kurtosis tensors to diagonize D so that D = [l1 0 0;
%0 l2 0; 0 0 l3] with l1 >= l2 >= l3 with l1,l2,l3 being the eigenvalues of
%D. 
%
%INPUTS: 
%    dt = [D11; D22; D33; D12; D13; D23];
%    kt = [W1111; W2222; W3333; W1112; W1113; W1222; W1333; W2223; W2333;
%    W1122; W1133; W2233; W1123; W1223; W1233]; 
%
%OUTPUTS: 
%    dt_rot = [l1; l2; l3]; 
%    kt_rot = [WT1111; WT2222; WT3333; WT1122; WT1133; WT2233]; with WT corresponding to the rotated kurtosis tensor
%    V; eigenvectors corresponding to l1; l2; l3 from the original diffusion tensor; ie D = V*diag(dt_rot)*V'; 
%
%Author: Russell Glenn
%Medical University of South Carolina

%ROTATE DIFFUSION AND KURTOSIS TENSORS
D = [dt(1) dt(4) dt(5);dt(4) dt(2) dt(6); dt(5) dt(6) dt(3)]; 
[V,L] = eig(D);
[L idx] = sort(diag(L),'descend'); 
% [L idx] = sort(diag(L),'ascend'); 
V = V(:,idx); 

dt_rot = [L;0;0;0];

W = zeros(3,3,3,3);
W(1,1,1,1) = kt(1);
W(2,2,2,2) = kt(2);
W(3,3,3,3) = kt(3);
W(1,1,1,2) = kt(4);  W(1,1,2,1) = W(1,1,1,2); W(1,2,1,1) = W(1,1,1,2); W(2,1,1,1) = W(1,1,1,2);
W(1,1,1,3) = kt(5);  W(1,1,3,1) = W(1,1,1,3); W(1,3,1,1) = W(1,1,1,3); W(3,1,1,1) = W(1,1,1,3);
W(1,2,2,2) = kt(6);  W(2,2,1,2) = W(1,2,2,2); W(2,1,2,2) = W(1,2,2,2); W(2,2,2,1) = W(1,2,2,2);
W(1,3,3,3) = kt(7);  W(3,3,1,3) = W(1,3,3,3); W(3,1,3,3) = W(1,3,3,3); W(3,3,3,1) = W(1,3,3,3);
W(2,2,2,3) = kt(8);  W(2,2,3,2) = W(2,2,2,3); W(2,3,2,2) = W(2,2,2,3); W(3,2,2,2) = W(2,2,2,3);
W(2,3,3,3) = kt(9);  W(3,2,3,3) = W(2,3,3,3); W(3,3,2,3) = W(2,3,3,3); W(3,3,3,2) = W(2,3,3,3);
W(1,1,2,2) = kt(10); W(1,2,1,2) = W(1,1,2,2); W(1,2,2,1) = W(1,1,2,2); W(2,1,2,1) = W(1,1,2,2); W(2,2,1,1) = W(1,1,2,2); W(2,1,1,2) = W(1,1,2,2); 
W(1,1,3,3) = kt(11); W(1,3,1,3) = W(1,1,3,3); W(1,3,3,1) = W(1,1,3,3); W(3,1,3,1) = W(1,1,3,3); W(3,3,1,1) = W(1,1,3,3); W(3,1,1,3) = W(1,1,3,3); 
W(2,2,3,3) = kt(12); W(2,3,2,3) = W(2,2,3,3); W(2,3,3,2) = W(2,2,3,3); W(3,2,3,2) = W(2,2,3,3); W(3,3,2,2) = W(2,2,3,3); W(3,2,2,3) = W(2,2,3,3); 
W(1,1,2,3) = kt(13); W(1,1,3,2) = W(1,1,2,3); W(1,2,1,3) = W(1,1,2,3); W(1,2,3,1) = W(1,1,2,3); W(1,3,1,2) = W(1,1,2,3); W(1,3,2,1) = W(1,1,2,3); W(2,1,1,3) = W(1,1,2,3); W(2,1,3,1) = W(1,1,2,3); W(2,3,1,1) = W(1,1,2,3); W(3,1,1,2) = W(1,1,2,3); W(3,1,2,1) = W(1,1,2,3); W(3,2,1,1) = W(1,1,2,3);
W(1,2,2,3) = kt(14); W(2,2,1,3) = W(1,2,2,3); W(2,2,3,1) = W(1,2,2,3); W(2,1,2,3) = W(1,2,2,3); W(2,1,3,2) = W(1,2,2,3); W(2,3,2,1) = W(1,2,2,3); W(2,3,1,2) = W(1,2,2,3); W(1,2,3,2) = W(1,2,2,3); W(1,3,2,2) = W(1,2,2,3); W(3,2,2,1) = W(1,2,2,3); W(3,2,1,2) = W(1,2,2,3); W(3,1,2,2) = W(1,2,2,3);
W(1,2,3,3) = kt(15); W(3,3,2,1) = W(1,2,3,3); W(3,3,1,2) = W(1,2,3,3); W(3,2,3,1) = W(1,2,3,3); W(3,2,1,3) = W(1,2,3,3); W(3,1,3,2) = W(1,2,3,3); W(3,1,2,3) = W(1,2,3,3); W(2,3,3,1) = W(1,2,3,3); W(2,3,1,3) = W(1,2,3,3); W(2,1,3,3) = W(1,2,3,3); W(1,3,3,2) = W(1,2,3,3); W(1,3,2,3) = W(1,2,3,3); 

kt_rot = zeros(15,1);
for i = 1:3; for j = 1:3; for k = 1:3; for l = 1:3
    kt_rot(1) = kt_rot(1) + V(i, 1) * V(j, 1) * V(k, 1) * V(l, 1) * W(i, j, k, l);  % Wt1111
    kt_rot(2) = kt_rot(2) + V(i, 2) * V(j, 2) * V(k, 2) * V(l, 2) * W(i, j, k, l);  % Wt2222
    kt_rot(3) = kt_rot(3) + V(i, 3) * V(j, 3) * V(k, 3) * V(l, 3) * W(i, j, k, l);  % Wt3333
    kt_rot(4) = kt_rot(4) + V(i, 1) * V(j, 1) * V(k, 1) * V(l, 2) * W(i, j, k, l);  % Wt1112
    kt_rot(5) = kt_rot(5) + V(i, 1) * V(j, 1) * V(k, 1) * V(l, 3) * W(i, j, k, l);  % Wt1113
    kt_rot(6) = kt_rot(6) + V(i, 1) * V(j, 2) * V(k, 2) * V(l, 2) * W(i, j, k, l);  % Wt1222
    kt_rot(7) = kt_rot(7) + V(i, 1) * V(j, 3) * V(k, 3) * V(l, 3) * W(i, j, k, l);  % Wt1333
    kt_rot(8) = kt_rot(8) + V(i, 2) * V(j, 2) * V(k, 2) * V(l, 3) * W(i, j, k, l);  % Wt2223
    kt_rot(9) = kt_rot(9) + V(i, 2) * V(j, 3) * V(k, 3) * V(l, 3) * W(i, j, k, l);  % Wt2333
    kt_rot(10) = kt_rot(10) + V(i, 1) * V(j, 1) * V(k, 2) * V(l, 2) * W(i, j, k, l);  % Wt1122
    kt_rot(11) = kt_rot(11) + V(i, 1) * V(j, 1) * V(k, 3) * V(l, 3) * W(i, j, k, l);  % Wt1133                
    kt_rot(12) = kt_rot(12) + V(i, 2) * V(j, 2) * V(k, 3) * V(l, 3) * W(i, j, k, l);  % Wt2233   
    kt_rot(13) = kt_rot(13) + V(i, 1) * V(j, 1) * V(k, 2) * V(l, 3) * W(i, j, k, l);  % Wt1123
    kt_rot(14) = kt_rot(14) + V(i, 1) * V(j, 2) * V(k, 2) * V(l, 3) * W(i, j, k, l);  % Wt1223
    kt_rot(15) = kt_rot(15) + V(i, 1) * V(j, 2) * V(k, 3) * V(l, 3) * W(i, j, k, l);  % Wt1233
end; end; end; end

if nargout>2
    varargout{1} = V; 
end