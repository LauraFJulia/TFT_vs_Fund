function Corresp=project3Dpoints(Points3D,CalM,R_t)

M=size(CalM,1)/3;% number of cameras
N=size(Points3D,2); %number of points to project

Corresp=zeros(2*M,N);
for m=1:M
    x=CalM(3*(m-1)+(1:3),:)*R_t(3*(m-1)+(1:3),:)*[Points3D;ones(1,N)];
    Corresp(2*(m-1)+(1:2),:)=bsxfun(@rdivide,x(1:2,:),x(3,:));
end
end