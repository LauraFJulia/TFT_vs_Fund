function [R_t_2,R_t_3,Reconst,T,iter]=PiPoseEstimation(Corresp,CalM)
% Pose estimation of 3 views from corresponding triplets of points using
% the Pi matrices from Ponce & Hebert parameterizing three views.
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
M=[null(P1), null(P2), null(P3)];
M=[M, null(M.')];
% transform camera matrices by M
P1=P1*M; P2=P2*M; P3=P3*M;

% find Pi matrices (vectors)
Pi1=inv(P1(:,2:4)); Pi2=inv(P2(:,[1 3 4])); Pi3=inv(P3(:,[1 2 4]));
Pi1=[0 0 0; Pi1];  
Pi2=[Pi2(1,:); 0 0 0; Pi2(2:3,:)];
Pi3=[Pi3(1:2,:); 0 0 0; Pi3(3,:)];
% minimal parameterization
Pi1=Pi1./(norm(Pi1(4,:))); Pi2=Pi2./(norm(Pi2(4,:))); Pi3=Pi3./(norm(Pi3(4,:)));
Q=eye(4);
Q(1,1)=1./norm(Pi3(1,:)-dot(Pi3(1,:),Pi3(4,:))*Pi3(4,:));   Q(1,4)=-Q(1,1)*dot(Pi3(1,:),Pi3(4,:));
Q(2,2)=1./norm(Pi1(2,:)-dot(Pi1(2,:),Pi1(4,:))*Pi1(4,:));   Q(2,4)=-Q(2,2)*dot(Pi1(2,:),Pi1(4,:));
Q(3,3)=1./norm(Pi2(3,:)-dot(Pi2(3,:),Pi2(4,:))*Pi2(4,:));   Q(3,4)=-Q(3,3)*dot(Pi2(3,:),Pi2(4,:));
Pi1=Q*Pi1;  Pi2=Q*Pi2;  Pi3=Q*Pi3;

P1=P1*inv(Q); P2=P2*inv(Q); P3=P3*inv(Q); 
points3D=triangulation3D({P1,P2,P3},[x1;x2;x3]);
p1_est=P1*points3D; p1_est=p1_est(1:2,:)./repmat(p1_est(3,:),2,1);
p2_est=P2*points3D; p2_est=p2_est(1:2,:)./repmat(p2_est(3,:),2,1);
p3_est=P3*points3D; p3_est=p3_est(1:2,:)./repmat(p3_est(3,:),2,1);

% minimize reproj error using Gauss-Helmert
pi=[reshape(Pi1(2:4,:).',9,1);reshape(Pi2([1 3 4],:).',9,1);reshape(Pi3([1 2 4],:).',9,1)];
x=reshape([x1;x2;x3],6*N,1);
x_est=reshape([p1_est;p2_est;p3_est],6*N,1);
y=zeros(0,1);
P=eye(6*N);
[~,pi_opt,~,iter]=Gauss_Helmert(@constraintsGH_PiM,x_est,pi,y,x,P);

% retrieve geometry from parameters 
Pi1=(reshape(pi_opt(1:9),3,3)).';
Pi2=(reshape(pi_opt(10:18),3,3)).';
Pi3=(reshape(pi_opt(19:27),3,3)).';
% P1=zeros(3,4); P2=zeros(3,4); P3=zeros(3,4);
% P1(:,2:4)=inv(Pi1); 
% P2(:,[1 3 4])=inv(Pi2);
% P3(:,[1 2 4])=inv(Pi3);
% T=TFT_from_P(P1,P2,P3);
% 
% % denormalization
% T=transform_TFT(T,Normal1,Normal2,Normal3,1);
% 
% % Find orientation using calibration and TFT
% [R_t_2,R_t_3]=R_t_from_T(T,CalM,Corresp);

F21=Normal2.'*(Pi2(2,:).'*Pi1(3,:)-Pi2(3,:).'*Pi1(2,:))*Normal1;
F31=Normal3.'*(Pi3(2,:).'*Pi1(3,:)-Pi3(3,:).'*Pi1(1,:))*Normal1;
[R2,t2]=recover_R_t(CalM(1:3,:),CalM(4:6,:),F21,Corresp(1:2,:),Corresp(3:4,:));
[R3,t3]=recover_R_t(CalM(1:3,:),CalM(7:9,:),F31,Corresp(1:2,:),Corresp(5:6,:));

% Find the norm of t3 using the image points and reconstruction from
% images 1 and 2
u3=CalM(7:9,:)*t3;
X=triangulation3D({CalM(1:3,:)*eye(3,4),CalM(4:6,:)*[R2,t2]},Corresp(1:4,:));
X=X(1:3,:)./repmat(X(4,:),3,1);
X3=CalM(7:9,:)*R3*X;
lam=-sum(dot(cross([Corresp(5:6,:);ones(1,N)],X3,1),cross([Corresp(5:6,:);ones(1,N)],repmat(u3,1,N)),1))/...
    sum(sum(cross([Corresp(5:6,:);ones(1,N)],repmat(u3,1,N)).^2));
t3=lam*t3;

R_t_2=[R2,t2]; R_t_3=[R3,t3];

% tensor
T=TFT_from_P(CalM(1:3,:)*eye(3,4),CalM(4:6,:)*R_t_2,CalM(7:9,:)*R_t_3);
% Find 3D points by triangulation
Reconst=triangulation3D({CalM(1:3,:)*eye(3,4),CalM(4:6,:)*R_t_2,CalM(7:9,:)*R_t_3},Corresp);

end



function [R_f,t_f]=recover_R_t(K1,K2,F21,x1,x2)

E21=K2.'*F21*K1;
W=[0 -1 0; 1 0 0; 0 0 1];
[U,~,V]=svd(E21);
R=U*W*V.';  Rp=U*W.'*V.';
R=R*sign(det(R)); Rp=Rp*sign(det(Rp));
t=U(:,3);

%from the 4 possible solutions find the correct one using the image points
num_points_seen=0;
for k=1:4
    if k==2 || k==4
        t=-t;
    elseif k==3
        R=Rp;
    end
    X1=triangulation3D({[K1 [0;0;0]],K2*[R,t]},[x1;x2]); X1=X1./repmat(X1(4,:),4,1);
    X2=[R t]*X1;
    if sum(sign(X1(3,:))+sign(X2(3,:)))>=num_points_seen
        R_f=R; t_f=t;
        num_points_seen=sum(sign(X1(3,:))+sign(X2(3,:)));
    end
end

end





