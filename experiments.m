% Script to recreate experiments on synthetic data for the ICIP paper 
%  submission "A Critical Review of the Trifocal Tensor Estimation"

clear; close all;

%% Initial parameters
N=12;       % number of 3D points
noise=1;    % sigma for the added Gaussian noise in pixels
f=50;       % focal length in mm
p_coll=0;   % no collinearity of camera centers
n_sim=10;   % number of simulations of data

%% Change the interval to reproduce different experiments

interval=0:0.25:3;      % for varying noise
% interval=20:20:300;     % for varying focal length
% interval=0.94:0.01:1;   % for making camera centers collinear
% interval=8:25;          % for varying number of initial points

%% Test the methods

method={...
    @LinearTFTPoseEstimation,...    % 1 - Linear TFT
    @MinimalTFTPoseEstimation,...   % 2 - Minimal TFT (Ressl)
    @NordbergTFTPoseEstimation,...  % 3 - Minimal TFT (Nordberg)
    @PapaFaugTFTPoseEstimation,...  % 4 - Papadopoulo Faugeras TFT
    @PiPoseEstimation,...           % 5 - Pi matrices (Ponce&Hebert)
    @PiColPoseEstimation,...        % 6 - Pi matrices - collinear (Ponce&Hebert)
    @LinearFPoseEstimation,...      % 7 - Linear Fundamental matrices
    @OptimFPoseEstimation};         % 8 - Fundamental matrices

methods_to_test=1:4;

% error vectors
repr_err=zeros(length(interval),length(method),2);
rot_err=zeros(length(interval),length(method),2);
t_err=zeros(length(interval),length(method),2);
iterBA=zeros(length(interval),length(method));

for i=1:length(interval)
    % change variable to the one to be varyied
    noise=interval(i);
    fprintf('Noise= %f\n', noise);
    
    for it=1:n_sim
        % Generate random data for a triplet of images
        [CalM,R_t0,Corresp]=generateSyntheticScene(N+100,noise,it,f,p_coll);
        rng(it);
        Corresp=Corresp(:,randsample(N+100,N));
        
        for m=methods_to_test
            
            % pose estimation by method m
            [R_t_2,R_t_3,Reconst]=method{m}(Corresp,CalM);
            % reprojection error
            repr_err(i,m,1)= repr_err(i,m,1)+...
                ReprError({CalM(1:3,:)*eye(3,4),...
                CalM(4:6,:)*R_t_2,CalM(7:9,:)*R_t_3},Corresp,Reconst)/n_sim;
            % angular errors
            [rot2_err,t2_err]=AngError(R_t0{1},R_t_2);
            [rot3_err,t3_err]=AngError(R_t0{2},R_t_3);
            
            rot_err(i,m,1)=rot_err(i,m,1)+(rot2_err+rot3_err)/(2*n_sim);
            t_err(i,m,1)=t_err(i,m,1)+(t2_err+t3_err)/(2*n_sim);
            
            % Apply Bundle Adjustment
            [R_t_ref,~,iter,repr_errBA]=BundleAdjustment(CalM,...
                [eye(3,4);R_t_2;R_t_3],Corresp,Reconst);
            % reprojection error
            repr_err(i,m,2)=repr_err(i,m,2)+repr_errBA/n_sim;
            % angular errors
            [rot2_err,t2_err]=AngError(R_t0{1},R_t_ref(4:6,:));
            [rot3_err,t3_err]=AngError(R_t0{2},R_t_ref(7:9,:));
            rot_err(i,m,2)=rot_err(i,m,2)+(rot2_err+rot3_err)/(2*n_sim);
            t_err(i,m,2)=t_err(i,m,2)+(t2_err+t3_err)/(2*n_sim);
            % iterations
            iterBA(i,m)=iterBA(i,m)+iter/n_sim;
        end
    end
end
            


%% Plot results
% methods_to_plot=1:6;        % All methods
% methods_to_plot=[1:3,5:6];  % no collinear method
methods_to_plot=methods_to_test;

method_names={'Linear TFT','Ressl TFT','Nordberg','PapadFaug','Ponce&Hebert',...
    'Ponce&Hebert-Col', 'Linear F', 'Optim F', 'Bundle Adj.'};

% reprojection error plot
figure;
plot(interval,repr_err(:,methods_to_plot,1))
title('Reprojection error')
legend(method_names(methods_to_plot),'Location','Best')

% rotation error plot
figure;
plot(interval,rot_err(:,methods_to_plot,1))
title('Angular error in rotations')
legend(method_names(methods_to_plot),'Location','Best')

% translation error plot
figure;
plot(interval,t_err(:,methods_to_plot,1))
title('Angular error in translations')
legend(method_names(methods_to_plot),'Location','Best')

% plots for Bundle Adjustment
% reprojection error plot
figure;
plot(interval,repr_err(:,methods_to_plot,2))
title('Reprojection error-BA')
legend(method_names(methods_to_plot),'Location','Best')

% rotation error plot
figure;
plot(interval,rot_err(:,methods_to_plot,2))
title('Angular error in rotations-BA')
legend(method_names(methods_to_plot),'Location','Best')

% translation error plot
figure;
plot(interval,t_err(:,methods_to_plot,2))
title('Angular error in translations-BA')
legend(method_names(methods_to_plot),'Location','Best')

% iterations in bundle adjustment plot
figure;
plot(interval,iterBA(:,methods_to_plot))
title('Iterations in bundle adjustment')
legend(method_names(methods_to_plot),'Location','Best')





