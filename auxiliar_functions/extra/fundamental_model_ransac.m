function [model,vector_errors]=fundamental_model_ransac(Data,sample)

p1=Data{1}; p2=Data{2};
N=size(p1,2);
F=linearF(p1(:,sample),p2(:,sample));

if size(p1,1)==2
    p1=[p1;ones(1,N)];
    p2=[p2;ones(1,N)];
end
    
vector_errors=zeros(N,1);
for i=1:N
    line=F*p1(:,i);
    dist2=abs(line.'*p2(:,i))/norm(line(1:2));
    line=F.'*p2(:,i);
    dist1=abs(line.'*p1(:,i))/norm(line(1:2));
    vector_errors(i)=max(dist1,dist2);
end

model=F;