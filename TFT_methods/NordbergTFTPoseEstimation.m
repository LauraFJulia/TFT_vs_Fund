function [R_t_2,R_t_3,Reconst,T,iter]=NordbergTFTPoseEstimation(Corresp,CalM)
% Pose estimation of 3 views from corresponding triplets of points using
% the minimal TriFocal Tensor parameterization by K. Nordberg.
%
% An initial trifocal tensor is computed linearly from the trilinearities
% using the triplets of correspondences. Then the error is minimized using
% Gauss-Helmert model to impose the minimal constraints of the K. Nordberg
% parameterization. After the optimization the essential matrices are
% computed from the tensor and the orientations are extracted by SVD.
%
% Input arguments:
%  Corresp  - 6xN matrix containing in each column, the 3 projections of
%             the same space point onto the 3 images.
%  CalM     - 9x3 matrix containing the M calibration 3x3 matrices for 
%             each camera concatenated.
%
% Output arguments: 
%  R_t_2    - 3x4 matrix containing the rotation matrix and translation 
%             vector [R2,t2] for the second camera.
%  R_t_3    - 3x4 matrix containing the rotation matrix and translation 
%             vector [R3,t3] for the third camera.
%  Reconst  - 3xN matrix containing the 3D reconstruction of the
%             correspondences.
%  T        - 3x3x3 array containing the trifocal tensor associated to 
%             this triplet of cameras.
% iter      - number of iterations needed in GH algorithm to reach minimum 
%
% Copyright (c) 2017 Laura F. Julia      

% Normalization of the data
[x1,Normal1]=Normalize2Ddata(Corresp(1:2,:));
[x2,Normal2]=Normalize2Ddata(Corresp(3:4,:));
[x3,Normal3]=Normalize2Ddata(Corresp(5:6,:));

% Model to estimate T: linear equations
[T,P1,P2,P3]=linearTFT(x1,x2,x3);

% apply projective transformation so matrix B ( P3=[B|b] ) has full rank
H=eye(4);
if rank(P3(:,1:3))<3 % should always be the case
    H(4,1:3)=null(P3(:,1:3))';
elseif rank(P2(:,1:3))<3
    H(4,1:3)=null(P2(:,1:3))';
end
P1=P1*H; P2=P2*H; P3=P3*H;
    
% Compute K. Nordberg param
A=P2(:,1:3); a=P2(:,4); r=A\a;
B=P3(:,1:3); b=P3(:,4); s=B\b;

U=[r, crossM(r)^2*s, crossM(r)*s];      U=U*(U'*U)^(-1/2); U=sign(det(U))*U;
V=[a, crossM(a)*A*s, crossM(a)^2*A*s];  V=V*(V'*V)^(-1/2); V=sign(det(V))*V;
W=[b, crossM(b)*B*r, crossM(b)^2*B*r];  W=W*(W'*W)^(-1/2); W=sign(det(W))*W;

% good representation of U V W
[~,~,v]=svd(U-eye(3));  vec_u=v(:,3);
o_u=atan2(vec_u'*[U(3,2)-U(2,3);U(1,3)-U(3,1);U(2,1)-U(1,2)]/2, (trace(U)-1)/2);
[~,~,v]=svd(V-eye(3));  vec_v=v(:,3);
o_v=atan2(vec_v'*[V(3,2)-V(2,3);V(1,3)-V(3,1);V(2,1)-V(1,2)]/2, (trace(V)-1)/2);
[~,~,v]=svd(W-eye(3));  vec_w=v(:,3);
o_w=atan2(vec_w'*[W(3,2)-W(2,3);W(1,3)-W(3,1);W(2,1)-W(1,2)]/2, (trace(W)-1)/2);

% param0=[vec_u*o_u; vec_v*o_v; vec_w*o_w];
% func=@(x)ls_orthog_matrices(x,T);
% param = lsqnonlin(func,param0);

Ts=transform_TFT(T,U,V,W,1);
paramT=Ts([1,7,10,12,16,19:22,25])';
[~,iT]=max(abs(paramT));
paramT=paramT/paramT(iT);

param0=[vec_u*o_u; vec_v*o_v; vec_w*o_w;...
    paramT([1:iT-1,iT+1:10])];
func=@(x)ls_trilin(x,iT,[x1;x2;x3]);
options = optimoptions(@lsqnonlin,'Display','off');
[param,~,~,~,output] = lsqnonlin(func,param0,[],[],options);
iter=output.iterations;

% orthogonal matrices
o_u=norm(param(1:3)); vec_u=param(1:3)/o_u;
o_v=norm(param(4:6)); vec_v=param(4:6)/o_v;
o_w=norm(param(7:9)); vec_w=param(7:9)/o_w;
U=eye(3)+sin(o_u)*crossM(vec_u)+(1-cos(o_u))*crossM(vec_u)^2;
V=eye(3)+sin(o_v)*crossM(vec_v)+(1-cos(o_v))*crossM(vec_v)^2;
W=eye(3)+sin(o_w)*crossM(vec_w)+(1-cos(o_w))*crossM(vec_w)^2;

% sparse tensor
paramT([1:iT-1,iT+1:10])=param(10:18);
Ts=zeros(3,3,3);
Ts([1,7,10,12,16,19:22,25])=paramT;
T=transform_TFT(Ts,U',V',W',1);

% denormalization
T= transform_TFT(T,Normal1,Normal2,Normal3,1);

% Find orientation using calibration and TFT
[R_t_2,R_t_3]=R_t_from_TFT(T,CalM,Corresp);

% Find 3D points by triangulation
Reconst=triangulation3D({CalM(1:3,:)*eye(3,4),CalM(4:6,:)*R_t_2,CalM(7:9,:)*R_t_3},Corresp);
Reconst=Reconst(1:3,:)./repmat(Reconst(4,:),3,1);

end

function f=ls_orthog_matrices(x,T0)

o_u=norm(x(1:3)); vec_u=x(1:3)/o_u;
o_v=norm(x(4:6)); vec_v=x(4:6)/o_v;
o_w=norm(x(7:9)); vec_w=x(7:9)/o_w;
U=eye(3)+sin(o_u)*crossM(vec_u)+(1-cos(o_u))*crossM(vec_u)^2;
V=eye(3)+sin(o_v)*crossM(vec_v)+(1-cos(o_v))*crossM(vec_v)^2;
W=eye(3)+sin(o_w)*crossM(vec_w)+(1-cos(o_w))*crossM(vec_w)^2;

Ts=transform_TFT(T0,U,V,W,1);

f=[ Ts(1,2,1); Ts(2,:,1)'; Ts(3,:,1)';...
    Ts(1,2,2); Ts(2,:,2)'; Ts(3,2:3,2)';...
    Ts(2:3,2,3); Ts(2:3,3,3)];
end
    
function f=ls_trilin(x,iT,obs)

% orthogonal matrices
o_u=norm(x(1:3)); vec_u=x(1:3)/o_u;
o_v=norm(x(4:6)); vec_v=x(4:6)/o_v;
o_w=norm(x(7:9)); vec_w=x(7:9)/o_w;
U=eye(3)+sin(o_u)*crossM(vec_u)+(1-cos(o_u))*crossM(vec_u)^2;
V=eye(3)+sin(o_v)*crossM(vec_v)+(1-cos(o_v))*crossM(vec_v)^2;
W=eye(3)+sin(o_w)*crossM(vec_w)+(1-cos(o_w))*crossM(vec_w)^2;

% sparse tensor
paramT=ones(10,1);
paramT([1:iT-1,iT+1:10])=x(10:18);
Ts=zeros(3,3,3);
Ts([1,7,10,12,16,19:22,25])=paramT;

% original tensor
T=transform_TFT(Ts,U',V',W',1);


N=size(obs,2);
f=zeros(4*N,1);
for i=1:N
    % points in the three images for correspondance i
    x1=obs(1:2,i); x2=obs(3:4,i); x3=obs(5:6,i);
    
    % 4 trilinearities
    ind2=4*(i-1);
    S2=[0 -1; -1 0; x2(2) x2(1)];
    S3=[0 -1; -1 0; x3(2) x3(1)];
    f(ind2+1:ind2+4)=reshape( S2.'*(x1(1)*T(:,:,1)+x1(2)*T(:,:,2)+T(:,:,3))*S3,4,1);
end

end
% 
% function f=ls_repr(x,iT,obs,CalM)
% 
% % orthogonal matrices
% o_u=norm(x(1:3)); vec_u=x(1:3)/o_u;
% o_v=norm(x(4:6)); vec_v=x(4:6)/o_v;
% o_w=norm(x(7:9)); vec_w=x(7:9)/o_w;
% U=eye(3)+sin(o_u)*crossM(vec_u)+(1-cos(o_u))*crossM(vec_u)^2;
% V=eye(3)+sin(o_v)*crossM(vec_v)+(1-cos(o_v))*crossM(vec_v)^2;
% W=eye(3)+sin(o_w)*crossM(vec_w)+(1-cos(o_w))*crossM(vec_w)^2;
% 
% % sparse tensor
% paramT=ones(10,1);
% paramT([1:iT-1,iT+1:10])=x(10:18);
% Ts=zeros(3,3,3);
% Ts([1,7,10,12,16,19:22,25])=paramT;
% 
% % original tensor
% T=transform_TFT(Ts,U',V',W',1);
% [R_t_2,R_t_3]=R_t_from_TFT(T,CalM,obs);
% P1=CalM(1:3,:)*eye(3,4);
% P2=CalM(1:3,:)*R_t_2;
% P3=CalM(7:9,:)*R_t_3;
% Pcam={P1,P2,P3};
% X=triangulation3D({P1,P2,P3},obs);
% 
% N=size(obs,2);
% f=zeros(2*3*N,1);
% for i=1:N
%     P=X(:,i);
%     for j=1:3
%         p=obs(2*(j-1)+1,i);
%         p_est=Pcam{j}*P; p_est=p_est(1:2)/p_est(3);
%         f(6*(i-1)+2*(j-1)+(1:2))=p-p_est;
%     end
% end
% 
% end




