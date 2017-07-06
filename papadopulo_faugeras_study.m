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

%%
clear;
syms u1 u2 u3 v1 v2 v3 w1 w2 w3
u=[u1;u2;u3]; v=[v1;v2;v3]; w=[w1;w2;w3];
o_u=sqrt(sum(u.^2)); vec_u=u/o_u;
o_v=sqrt(sum(v.^2)); vec_v=v/o_v;
o_w=sqrt(sum(w.^2)); vec_w=w/o_w;
%U=eye(3)+sin(o_u)*crossM(vec_u)+(1-cos(o_u))*crossM(vec_u)^2;
U=eye(3)*cos(o_u)+sin(o_u)*crossM(vec_u)+(1-cos(o_u))*(vec_u*vec_u.');
V=eye(3)*cos(o_v)+sin(o_v)*crossM(vec_v)+(1-cos(o_v))*(vec_v*vec_v.');
W=eye(3)*cos(o_w)+sin(o_w)*crossM(vec_w)+(1-cos(o_w))*(vec_w*vec_w.');

U=simplify(U); V=simplify(V); W=simplify(W);

syms p1 p2 p3 p4 p5 p6 p7 p8 p9 p10
t=[p1 p2 p3 p4 p5 p6 p7 p8 p9 p10].';
Ts=zeros(3,3,3)*p1;
param_ind=[1,7,10,12,16,19:22,25];
Ts(param_ind)=t;

% original tensor
T=transf_t(Ts,U.',V.',W.');
T=simplify(T);

J_true=zeros(27,19)*p1;
var=[u;v;w;t];
for i=1:19
    J_true(:,i)=diff(T(:),var(i));
end


J_T=zeros(27,10)*p1;
for i=1:10
    J_T(:,i)=diff(T(:),t(i));
end



J=zeros(27,19)*p1;
% derivatives of the parameterization w.r.t. the sparse tensor
for i=1:10
    e=zeros(3,3,3)*p1; e(param_ind(i))=1;
    J(:,i+9)=reshape(transf_t(e,U.',V.',W.'),27,1);
end
% derivatives of the parameterization w.r.t. the orthogonal matrices
dU=zeros(3,3,3)*p1; dV=zeros(3,3,3)*p1; dW=zeros(3,3,3)*p1;
dU_true=zeros(3,3,3)*p1; dV_true=zeros(3,3,3)*p1; dW_true=zeros(3,3,3)*p1;
e=eye(3);
for i=1:3
    dU(:,:,i)=-vec_u(i)*sin(o_u)*eye(3)+ vec_u(i)*cos(o_u)*crossM(vec_u)+...
        sin(o_u)*(1/o_u)*(crossM(e(:,i))- vec_u(i)*crossM(vec_u))+...
        vec_u(i)*sin(o_u)*(vec_u*vec_u.')+...
        (1-cos(o_u))*(1/o_u)*(vec_u*e(i,:)+e(:,i)*vec_u.'-2*vec_u(i)*(vec_u*vec_u.'));
    dU_true(:,:,i)=diff(U,u(i));
    
    dV(:,:,i)=-vec_v(i)*sin(o_v)*eye(3)+ vec_v(i)*cos(o_v)*crossM(vec_v)+...
        sin(o_v)*(1/o_v)*(crossM(e(:,i))- vec_v(i)*crossM(vec_v))+...
        vec_v(i)*sin(o_v)*(vec_v*vec_v.')+...
        (1-cos(o_v))*(1/o_v)*(vec_v*e(i,:)+e(:,i)*vec_v.'-2*vec_v(i)*(vec_v*vec_v.'));
    dV_true(:,:,i)=diff(V,v(i));
    
    dW(:,:,i)=-vec_w(i)*sin(o_w)*eye(3)+ vec_w(i)*cos(o_w)*crossM(vec_w)+...
        sin(o_w)*(1/o_w)*(crossM(e(:,i))- vec_w(i)*crossM(vec_w))+...
        vec_w(i)*sin(o_w)*(vec_w*vec_w.')+...
        (1-cos(o_w))*(1/o_w)*(vec_w*e(i,:)+e(:,i)*vec_w.'-2*vec_w(i)*(vec_w*vec_w.'));
    dW_true(:,:,i)=diff(W,w(i));
end


    
for i=1:3
    dTu=transf_t(Ts,dU(:,:,i).',V.',W.');
    dTv=transf_t(Ts,U.',dV(:,:,i).',W.');
    dTw=transf_t(Ts,U.',V.',dW(:,:,i).');
    J(:,i+[0,3,6])=[dTu(:),dTv(:),dTw(:)];
end

simplify(J_true-J)

j=J_true-J; j=j(:);

%%

T0=rand(3,3,3);

u=rand(3,1); v=rand(3,1); w=rand(3,1);
o_u=norm(u); vec_u=u/o_u;
o_v=norm(v); vec_v=v/o_v;
o_w=norm(w); vec_w=w/o_w;
U=eye(3)+sin(o_u)*crossM(vec_u)+(1-cos(o_u))*crossM(vec_u)^2;
V=eye(3)+sin(o_v)*crossM(vec_v)+(1-cos(o_v))*crossM(vec_v)^2;
W=eye(3)+sin(o_w)*crossM(vec_w)+(1-cos(o_w))*crossM(vec_w)^2;

T=zeros(3,3,3);
for i=1:3
    T(:,:,i)=V' *(U(1,i)*T0(:,:,1) +U(2,i)*T0(:,:,2) +U(3,i)*T0(:,:,3))* W;
end

t=reshape(T,27,1,1);




