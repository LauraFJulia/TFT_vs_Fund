function [R_t_2,R_t_3,Reconst,T,iter]=PiColPoseEstimation(Corresp,CalM)
% Pose estimation of 3 views from corresponding triplets of points using
% the Pi matrices from Ponce & Hebert parameterizing three views in the
% special case of collinear camera centers.
%
% An initial trifocal tensor is computed linearly from the trilinearities
% using the triplets of correspondences. From it, initial Pi matrices are
% computed. Then the error is minimized using the Gauss-Helmert model to
% impose the minimal constraints of the Pose & Hebert Pi matrices 
% parameterization. After the optimization, a final TFT is computed and the
% essential matrices are extracted. Finally, the orientations are retrieved
% by SVD.
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

% Number of correspondences
N=size(Corresp,2);

% Normalization of the data
[x1,Normal1]=Normalize2Ddata(Corresp(1:2,:));
[x2,Normal2]=Normalize2Ddata(Corresp(3:4,:));
[x3,Normal3]=Normalize2Ddata(Corresp(5:6,:));

% First approximation of T: linear equations
[~,P1,P2,P3]=linearTFT(x1,x2,x3);

%%% Epipolar-Trinocular error with Gauss-Helmert
% find homography H sending camera centers to fundamental points
M=[null(P1), null(P2)];
coeff=M\null(P3);
M=[coeff(1)*M(:,1),coeff(2)*M(:,2), null(M.')];
% transform camera matrices by M
P1=P1*M; P2=P2*M; P3=P3*M;

% find Pi matrices (vectors)
Pi1=inv(P1(:,2:4)); Pi2=inv(P2(:,[1 3 4])); Pi3=inv(P3(:,[2 3 4]));
Pi1=[0 0 0; Pi1];  
Pi2=[Pi2(1,:); 0 0 0;Pi2(2:3,:)];
Pi3=[0 0 0; Pi3];
% minimal parameterization
Pi1=Pi1./(norm(Pi1(4,:))); Pi2=Pi2./(norm(Pi2(4,:))); Pi3=Pi3./(norm(Pi3(4,:)));
Q1=eye(4);
u1=Pi1(3,:).';  v1=Pi1(4,:).';
u2=Pi2(3,:).';  v2=Pi2(4,:).';
tol=1e-10;
% Make pi31 perpendicular to pi41 and 
% pi32 perpendicular to pi42
A=(v1.'*v1)*(u2.'*v2)-(u1.'*v1)*(v2.'*v2);
B=(v1.'*v1)*(u2.'*u2)-(u1.'*u1)*(v2.'*v2);
C=(u1.'*v1)*(u2.'*u2)-(u1.'*u1)*(u2.'*v2);
if abs(A)>tol && (B^2-4*A*C)>=0 && abs(C)>tol
    Q1(3,4)=(-B+sqrt(B^2-4*A*C))/(2*A);
    Q1(4,3)=(-B+sqrt(B^2-4*A*C))/(2*C);
else
    error('The minimal param could not be found\n');
end
A=(u1*v1.'-v1*u1.'); B=(u2*v2.'-v2*u2.');
Q1(2,4)=(Pi1(2,:)*A*u1)/(u1.'*A*v1);
Q1(2,3)=(Pi1(2,:)*A.'*v1)/(u1.'*A*v1);
Q1(1,4)=(Pi2(1,:)*B*u2)/(u1.'*B*v2);
Q1(1,3)=(Pi2(1,:)*B.'*v2)/(u1.'*B*v2);

Pi1=Q1*Pi1;  Pi2=Q1*Pi2;  Pi3=Q1*Pi3;
Pi1=Pi1/norm(Pi1(2,:));
Pi2=Pi2/norm(Pi2(1,:));
Pi3=Pi3/norm(Pi3(2,:)-Pi3(1,:));
Q2=eye(4); 
Q2(3,3)=1/norm(Pi3(3,:));
Q2(4,4)=1/norm(Pi3(4,:));
Pi1=Q2*Pi1;  Pi2=Q2*Pi2;  Pi3=Q2*Pi3;
Pi3(1:2,:)=Pi3(1:2,:)-Pi3([1 1],:);

P1=P1*inv(Q2*Q1); P2=P2*inv(Q2*Q1); P3=P3*inv(Q2*Q1); 
points3D=triangulation3D({P1,P2,P3},[x1;x2;x3]);
p1_est=P1*points3D; p1_est=p1_est(1:2,:)./repmat(p1_est(3,:),2,1);
p2_est=P2*points3D; p2_est=p2_est(1:2,:)./repmat(p2_est(3,:),2,1);
p3_est=P3*points3D; p3_est=p3_est(1:2,:)./repmat(p3_est(3,:),2,1);

% minimize reproj error using Gauss-Helmert
pi=[reshape(Pi1(2:4,:).',9,1);reshape(Pi2([1 3 4],:).',9,1);reshape(Pi3(2:4,:).',9,1)];
x=reshape([x1(1:2,:);x2(1:2,:);x3(1:2,:)],6*N,1);
x_est=reshape([p1_est;p2_est;p3_est],6*N,1);
y=zeros(0,1);
P=eye(6*N);
[~,pi_opt,~,iter]=Gauss_Helmert(@constraintsGH_PiCM,x_est,pi,y,x,P);

% retrieve geometry from parameters 
Pi1=(reshape(pi_opt(1:9),3,3)).';
Pi2=(reshape(pi_opt(10:18),3,3)).';
Pi3=(reshape(pi_opt(19:27),3,3)).';
P1=zeros(3,4); P2=zeros(3,4); P3=zeros(3,4);
P1(:,2:4)=inv(Pi1); 
P2(:,[1 3 4])=inv(Pi2);
P3(:,2:4)=inv(Pi3);P3(:,1)=-P3(:,2);
T=TFT_from_P(P1,P2,P3);

% denormalization
T= transform_TFT(T,Normal1,Normal2,Normal3,1);

% Find orientation using calibration and TFT
[R_t_2,R_t_3]=R_t_from_T(T,CalM,Corresp);

% Find 3D points by triangulation
Reconst=triangulation3D({CalM(1:3,:)*eye(3,4),CalM(4:6,:)*R_t_2,CalM(7:9,:)*R_t_3},Corresp);
end







