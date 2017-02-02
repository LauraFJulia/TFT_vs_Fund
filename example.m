% Example on how to use the code

clear;

%% Generate some random data for a triplet of images
N=100;      % number of 3D points
noise=1;    % sigma for the added Gaussian noise in pixels
seed=1;     % seed for random generation
f=50;       % focal length in mm
p_coll=0;   % no collinearity of camera centers
[CalM,R_t0,Corresp,points3D]=generateSyntheticScene(N,noise,seed,f,p_coll);

%% Test a method

method={...
    @LinearTFTPoseEstimation,...    % 1 - Linear TFT
    @MinimalTFTPoseEstimation,...   % 2 - Minimal TFT (Ressl)
    @PiPoseEstimation,...           % 3 - Pi matrices (Ponce&Hebert)
    @PiPoseEstimation,...           % 4 - Pi matrices - collinear (Ponce&Hebert)
    @LinearFPoseEstimation,...      % 5 - Linear Fundamental matrices
    @OptimFPoseEstimation};         % 6 - Fundamental matrices

[R_t_2,R_t_3,Reconst]=method{1}(Corresp,CalM);

%% Compute the errors

% reprojection error
repr_err=ReprError({CalM(1:3,:)*eye(3,4),CalM(4:6,:)*R_t_2,CalM(7:9,:)*R_t_3},...
    Corresp,Reconst);
fprintf('Reprojection error is %f\n',repr_err);

% angular errors
[rot2_err,t2_err]=AngError(R_t0{1},R_t_2);
[rot3_err,t3_err]=AngError(R_t0{2},R_t_3);

fprintf('Angular errors in rotations are %f and %f, for translations %f and %f\n',...
    rot2_err,rot3_err,t2_err,t3_err);

%% Apply Bundle Adjustment

[R_t_ref,Reconst_ref,iter,repr_err]=BundleAdjustment(CalM,[eye(3,4);R_t_2;R_t_3],Corresp,Reconst);
fprintf('Reprojection error is %f after Bundle Adjustment\n',repr_err);

% angular errors
[rot2_err,t2_err]=AngError(R_t0{1},R_t_ref(4:6,:));
[rot3_err,t3_err]=AngError(R_t0{2},R_t_ref(7:9,:));
fprintf('Angular errors in rotations after BA are %f and %f, for translations %f and %f\n',...
    rot2_err,rot3_err,t2_err,t3_err);




