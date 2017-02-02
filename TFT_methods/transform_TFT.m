function T_new = transform_TFT(T_old,M1,M2,M3,inverse)
% Tranformed TFT
%
% short function to transform the TFT when an algebraic transformation
% has been aplied to the image points.
%
% if inverse==0 :
%   from a TFT T_old assossiated to P1_old, P2_old, P_old, find the new TFT
%   T_new associated to P1_new=M1*P1_old, P2_new=M2*P2_old, P3_new=M3*P3_old.
% if inverse==1 :
%   from a TFT T_old assossiated to P1_old, P2_old, P_old, find the new TFT
%   T_new associated to M1*P1_new=P1_old, M2*P2_new=P2_old, M3*P3_new=P3_old.

if inverse==0
    M1i=inv(M1); T_new=zeros(3,3,3);
    T_new(:,:,1)=M2*(M1i(1,1)*T_old(:,:,1) + M1i(2,1)*T_old(:,:,2) + M1i(3,1)*T_old(:,:,3) )*M3.';
    T_new(:,:,2)=M2*(M1i(1,2)*T_old(:,:,1) + M1i(2,2)*T_old(:,:,2) + M1i(3,2)*T_old(:,:,3) )*M3.';
    T_new(:,:,3)=M2*(M1i(1,3)*T_old(:,:,1) + M1i(2,3)*T_old(:,:,2) + M1i(3,3)*T_old(:,:,3) )*M3.';

elseif inverse==1
    M2i=inv(M2); M3i=inv(M3); T_new=zeros(3,3,3);
    T_new(:,:,1)=M2i*(M1(1,1)*T_old(:,:,1) + M1(2,1)*T_old(:,:,2) + M1(3,1)*T_old(:,:,3) )*M3i.';
    T_new(:,:,2)=M2i*(M1(1,2)*T_old(:,:,1) + M1(2,2)*T_old(:,:,2) + M1(3,2)*T_old(:,:,3) )*M3i.';
    T_new(:,:,3)=M2i*(M1(1,3)*T_old(:,:,1) + M1(2,3)*T_old(:,:,2) + M1(3,3)*T_old(:,:,3) )*M3i.';
end

T_new=T_new/norm(T_new(:));

end
