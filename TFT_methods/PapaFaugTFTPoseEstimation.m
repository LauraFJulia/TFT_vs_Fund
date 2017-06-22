function [R_t_2,R_t_3,Reconst,T,iter]=PapaFaugTFTPoseEstimation(Corresp,CalM)
% Pose estimation of 3 views from corresponding triplets of points using
% the TriFocal Tensor with minimal constraints by T. Papadopoulo and O.
% Faugeras.
%
% An initial trifocal tensor is computed linearly from the trilinearities
% using the triplets of correspondences. Then a minimal parameterization is
% computed using the constraints and minimisation model presented in "A non
% linear method for estimatin the projective geometry of three views" by T.
% Papadopoulo and O. Faugeras. After the optimization the essential 
% matrices are computed from the tensor and the orientations are extracted
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
%
% Copyright (c) 2017 Laura F. Julia      

% Normalization of the data
[x1,Normal1]=Normalize2Ddata(Corresp(1:2,:));
[x2,Normal2]=Normalize2Ddata(Corresp(3:4,:));
[x3,Normal3]=Normalize2Ddata(Corresp(5:6,:));

% Model to estimate T: linear equations
[T,~,~,~]=linearTFT(x1,x2,x3);

param0=T(:);
lam=10^5;
func=@(x)ls_trilin_constr(x,[x1;x2;x3],lam);
options = optimoptions(@lsqnonlin,'Display','off');
param = lsqnonlin(func,param0,[],[],options);
T=reshape(param,3,3,3);

% denormalization
T= transform_TFT(T,Normal1,Normal2,Normal3,1);

% Find orientation using calibration and TFT
[R_t_2,R_t_3]=R_t_from_TFT(T,CalM,Corresp);

% Find 3D points by triangulation
Reconst=triangulation3D({CalM(1:3,:)*eye(3,4),CalM(4:6,:)*R_t_2,CalM(7:9,:)*R_t_3},Corresp);
Reconst=Reconst(1:3,:)./repmat(Reconst(4,:),3,1);

end





function f=ls_trilin_constr(x,obs,lam)


% sparse 
T=reshape(x,3,3,3);

N=size(obs,2);
f=zeros(4*N + 12,1);
for i=1:N
    % points in the three images for correspondance i
    x1=obs(1:2,i); x2=obs(3:4,i); x3=obs(5:6,i);
    
    % 4 trilinearities
    ind2=4*(i-1);
    S2=[0 -1; -1 0; x2(2) x2(1)];
    S3=[0 -1; -1 0; x3(2) x3(1)];
    f(ind2+1:ind2+4)=reshape( S2.'*(x1(1)*T(:,:,1)+x1(2)*T(:,:,2)+T(:,:,3))*S3,4,1);
end

% constraints with weight lambda
for i=1:3
    f(4*N +i,:)=lam*det(T(:,:,i));
end

i=0;
for k2=1:2
    for k3=1:2
        for l2=k2+1:3
            for l3=k3+1:3
                i=i+1;
                f(4*N +3 +i,1)= lam*(...
                    det(reshape([T(k2,k3,:),T(k2,l3,:),T(l2,l3,:)],3,3))*...
                    det(reshape([T(k2,k3,:),T(l2,k3,:),T(l2,l3,:)],3,3))-...
                    det(reshape([T(l2,k3,:),T(k2,l3,:),T(l2,l3,:)],3,3))*...
                    det(reshape([T(k2,k3,:),T(l2,k3,:),T(k2,l3,:)],3,3)));
            end
        end
    end
end


end


