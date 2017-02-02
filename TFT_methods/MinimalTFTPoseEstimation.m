function [R_t_2,R_t_3,Reconst,T,iter]=MinimalTFTPoseEstimation(Corresp,CalM)
% Pose estimation of 3 views from corresponding triplets of points using
% the minimal TriFocal Tensor.
%
% An initial trifocal tensor is computed linearly from the trilinearities
% using the triplets of correspondences. Then the error is minimized using
% Gauss-Helmert model to impose the minimal constraints of the C. Ressl TFT
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

% Normalization of the data
[x1,Normal1]=Normalize2Ddata(Corresp(1:2,:));
[x2,Normal2]=Normalize2Ddata(Corresp(3:4,:));
[x3,Normal3]=Normalize2Ddata(Corresp(5:6,:));

% Model to estimate T: linear equations
[T,P1,P2,P3]=linearTFT(x1,x2,x3);

% Compute Ressl param
e21=P2(:,4);
[~,Ind]=max(abs(e21)); e21=e21/e21(Ind);
e31=P3(:,4);           e31=e31/norm(e31);

S=[T(Ind,:,1).',T(Ind,:,2).',T(Ind,:,3).'];
aux=norm(S(:));     S=S./aux;
T=T/aux;

Ind2=1:3; Ind2=Ind2(Ind2~=Ind);
mn=[e31.'*(T(:,:,1).'-S(:,1)*e21.');...
    e31.'*(T(:,:,2).'-S(:,2)*e21.');...
    e31.'*(T(:,:,3).'-S(:,3)*e21.')];
mn=mn(:,Ind2);

% compute 3d estimated points to have initial estimated reprojected image
% points
points3D=triangulation3D({P1,P2,P3},[x1;x2;x3]);
p1_est=P1*points3D; p1_est=p1_est(1:2,:)./repmat(p1_est(3,:),2,1);
p2_est=P2*points3D; p2_est=p2_est(1:2,:)./repmat(p2_est(3,:),2,1);
p3_est=P3*points3D; p3_est=p3_est(1:2,:)./repmat(p3_est(3,:),2,1);

% minimize reprojection error with Gauss-Helmert
N=size(x1,2);
p=[S(:);e21(Ind2);mn(:);e31];
x=reshape([x1(1:2,:);x2(1:2,:);x3(1:2,:)],6*N,1);
x_est=reshape([p1_est;p2_est;p3_est],6*N,1);
y=zeros(0,1);
func=@(x1,x2,x3)constraintsGH_MTFT(x1,x2,x3,Ind);
[~,p_opt,~,iter]=Gauss_Helmert(func,x_est,p,y,x,eye(6*N));

% recover parameters
S=reshape(p_opt(1:9),3,3);
e21=ones(3,1); e21(Ind2)=p_opt(10:11);
mn=zeros(3);   mn(:,Ind2)=reshape(p_opt(12:17),3,2);
e31=p_opt(18:20);
% estimated tensor
T(:,:,1)=(S(:,1)*e21.'+e31*mn(1,:)).';
T(:,:,2)=(S(:,2)*e21.'+e31*mn(2,:)).';
T(:,:,3)=(S(:,3)*e21.'+e31*mn(3,:)).';
% denormalization
T= transform_TFT(T,Normal1,Normal2,Normal3,1);

% Find orientation using calibration and TFT
[R_t_2,R_t_3]=R_t_from_T(T,CalM,Corresp);

% Find 3D points by triangulation
Reconst=triangulation3D({CalM(1:3,:)*eye(3,4),CalM(4:6,:)*R_t_2,CalM(7:9,:)*R_t_3},Corresp);

end




