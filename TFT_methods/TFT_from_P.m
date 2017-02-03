function T=TFT_from_P(P1,P2,P3)
% Trifocal tensor from the projection matrices
%
% General formula to compute the trifocal tensor from any three projection
% matrices.
%
% Copyright (c) 2017 Laura F. Julia

T=zeros(3,3,3);
for i=1:3
    for j=1:3
        for k=1:3
            T(j,k,i)=(-1)^(i+1)*det([P1([1:(i-1) (i+1):3],:);P2(j,:);P3(k,:)]);
        end
    end
end
T=T/norm(T(:));
end
