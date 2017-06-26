syms t1 t2 t3 t4 t5 t6 t7 t8 t9 t10 t11 t12 t13 t14 t15 t16 t17 t18...
    t19 t20 t21 t22 t23 t24 t25 t26 t27
p= [t1 t2 t3 t4  t5 t6 t7 t8 t9 t10 t11 t12 t13 t14 t15 t16 t17 t18...
    t19 t20 t21 t22 t23 t24 t25 t26 t27].';
T=reshape(p,3,3,3);

C2=zeros(9,1)*t1;
i=0;
for k2=1:2
    for k3=1:2
        for l2=k2+1:3
            for l3=k3+1:3
                i=i+1;
                C2(i,1)= det(reshape([T(k2,k3,:),T(k2,l3,:),T(l2,l3,:)],3,3))*...
                    det(reshape([T(k2,k3,:),T(l2,k3,:),T(l2,l3,:)],3,3))-...
                    det(reshape([T(l2,k3,:),T(k2,l3,:),T(l2,l3,:)],3,3))*...
                    det(reshape([T(k2,k3,:),T(l2,k3,:),T(k2,l3,:)],3,3));
            end
        end
    end
end

J=zeros(9,27)*t1;
for i=1:27
    J(:,i)=diff(C2,p(i));
end
simplify(J)
