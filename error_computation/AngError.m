function [rot_err,t_err]=AngError(R_t_true,R_t_est)

R_true=R_t_true(:,1:3);     t_true=R_t_true(:,4);
R_est=R_t_est(:,1:3);       t_est=R_t_est(:,4);

% determine angle difference between rotations
rot_err=abs(180*acos((trace(R_true.'*R_est)-1)/2)/pi);

% determine angle difference between translations
t_err=abs(180*acos(dot( t_est/norm(t_est), t_true/norm(t_true) ))/pi);
end
