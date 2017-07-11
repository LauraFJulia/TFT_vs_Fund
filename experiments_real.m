% Script to recreate experiments on real data for the PSIVT paper
%  submission "A Critical Review of the Trifocal Tensor Estimation"

clear; close all;


%% Here uncomment the dataset to use.

% dataset='fountain-P11';  
dataset='Herz-Jesu-P8';

%% Some parameters

path_to_data=strcat('Data/',dataset,'/');

switch dataset
    case 'fountain-P11'
        triplets_to_test=1:70;
    case 'Herz-Jesu-P8'
        triplets_to_test=1:50;
end

initial_sample_size=100;
bundle_adj_size=50;
repr_err_th=1;


%% Recover correspondances

corresp_file = matfile(strcat(path_to_data,'Corresp_triplets','.mat'));
indexes_sorted = corresp_file.indexes_sorted;
corresp_by_triplet = corresp_file.Corresp;
im_names = corresp_file.im_names;
clear corresp_file;

%% methods to evaluate

methods={...
    @LinearTFTPoseEstimation,...    % 1 - TFT - Linear estimation
    @ResslTFTPoseEstimation,...     % 2 - TFT - Ressl
    @NordbergTFTPoseEstimation,...  % 3 - TFT - Nordberg
    @PapaFaugTFTPoseEstimation,...  % 4 - TFT - Papadopoulo&Faugeras 
    @PiPoseEstimation,...           % 5 - Pi matrices - Ponce&Hebert
    @PiColPoseEstimation,...        % 6 - Pi matrices - Ponce&Hebert for collinear cameras
    @LinearFPoseEstimation,...      % 7 - Fundamental matrices - Linear estimation
    @OptimFPoseEstimation};         % 8 - Fundamental matrices - Optimized

methods_to_test=[1:5,7:8]; % no method for collinear cameras

%% error vectors
repr_err = zeros(length(triplets_to_test),length(methods),2);
rot_err  = zeros(length(triplets_to_test),length(methods),2);
t_err    = zeros(length(triplets_to_test),length(methods),2);
iter     = zeros(length(triplets_to_test),length(methods),2);
time     = zeros(length(triplets_to_test),length(methods),2);

%% evaluation

for it=1:length(triplets_to_test)
    
    % Triplet information and correspondances
    triplet=indexes_sorted(triplets_to_test(it) ,1:3); 
    im1=triplet(1); im2=triplet(2);  im3=triplet(3);
    Corresp=corresp_by_triplet{im1,im2,im3}';
    N=size(Corresp,2);
    fprintf('Triplet %d/%d (%d,%d,%d) with %d matching points.\n',...
        it,length(triplets_to_test),im1,im2,im3,N);
    
    % Ground truth poses and calibration
    [K1,R1_true,t1_true,im_size]=readCalibrationOrientation_Strecha(path_to_data,im_names{im1});
    [K2,R2_true,t2_true]=readCalibrationOrientation_Strecha(path_to_data,im_names{im2});
    [K3,R3_true,t3_true]=readCalibrationOrientation_Strecha(path_to_data,im_names{im3});
    CalM=[K1;K2;K3];
    R_t0={[R2_true*R1_true.', t2_true-R2_true*R1_true.'*t1_true],...
          [R3_true*R1_true.', t3_true-R3_true*R1_true.'*t1_true]};
    
    % Discart correspondances with repr_err > 1 pix
    Reconst0=triangulation3D({K1*eye(3,4),K2*R_t0{1},K3*R_t0{2}},Corresp);
    Reconst0=bsxfun(@rdivide,Reconst0(1:3,:),Reconst0(4,:));
    Corresp_new=project3Dpoints(Reconst0,CalM,[eye(3,4);R_t0{1};R_t0{2}]);
    residuals=Corresp_new-Corresp; % reprojection error
    Corresp_inliers=Corresp(:, sum(abs(residuals)>repr_err_th,1)==0 );
    N=size(Corresp_inliers,2);
    REr=ReprError({K1*eye(3,4),K2*R_t0{1},K3*R_t0{2}},Corresp_inliers);
    fprintf('%d valid correspondances with reprojection error %f.\n',N,REr);

    % Samples for initial estimation and bundle adjustment
    rng(it);
    init_sample=randsample(1:N,min(initial_sample_size,N));
    rng(it);
    ref_sample=randsample(init_sample,min(bundle_adj_size,length(init_sample)));
    Corresp_init=Corresp(:,init_sample);
    Corresp_ref=Corresp(:,ref_sample);
    
    % estimate the pose using all methods
    fprintf('method... ');
    for m=methods_to_test
        fprintf('%d ',m);
        
        % if there are not enough matches for initial estimation, inf error
        if (m>6 && N<8) || N<7 
            repr_err(it,m,:)=inf;    rot_err(it,m,:)=inf;
            t_err(it,m,:)=inf;       iter(it,m,:)=inf;
            time(it,m,:)=inf;
            continue;
        end
        
        % % Pose estimation by method m, measuring time            
        t0=cputime;
        [R_t_2,R_t_3,~,~,nit]=methods{m}(Corresp_init,CalM);
        t=cputime-t0;
        
        % reprojection error with all inliers
        repr_err(it,m,1)= ReprError({CalM(1:3,:)*eye(3,4),...
            CalM(4:6,:)*R_t_2,CalM(7:9,:)*R_t_3},Corresp_inliers);
        % angular errors
        [rot2_err,t2_err]=AngError(R_t0{1},R_t_2);
        [rot3_err,t3_err]=AngError(R_t0{2},R_t_3);
        rot_err(it,m,1)=(rot2_err+rot3_err)/2;
        t_err(it,m,1)=(t2_err+t3_err)/2;
        % iterations & time
        iter(it,m,1)=nit; time(it,m,1)=t;
        
        
        % % Apply Bundle Adjustment
        fprintf('(ref)... ');
        t0=cputime;
        [R_t_ref,~,nit,repr_errBA]=BundleAdjustment(CalM,...
            [eye(3,4);R_t_2;R_t_3],Corresp_ref);
        t=cputime-t0;

        % reprojection error with all inliers
        repr_err(it,m,2)= ReprError({CalM(1:3,:)*R_t_ref(1:3,:),...
            CalM(4:6,:)*R_t_ref(4:6,:),...
            CalM(7:9,:)*R_t_ref(7:9,:)},Corresp_inliers);
        % angular errors
        [rot2_err,t2_err]=AngError(R_t0{1},R_t_ref(4:6,:));
        [rot3_err,t3_err]=AngError(R_t0{2},R_t_ref(7:9,:));
        rot_err(it,m,2)=(rot2_err+rot3_err)/2;
        t_err(it,m,2)=(t2_err+t3_err)/2;
        % iterations & time
        iter(it,m,2)=nit; time(it,m,2)=t;
        
    end
    fprintf('\n');
    
end

%%

save Results/errors_1_to_70_triplets_fountain_100i_50r_1_reprinliers.mat repr_err rot_err t_err iter time...
   methods_to_test methods dataset triplets_to_test initial_sample_size...
   bundle_adj_size repr_err_th

% save Results/errors_1_to_50_triplets_herzjesu_100i_50r_1_reprinliers.mat repr_err rot_err t_err iter time...
%     methods_to_test method dataset triplets_to_evaluate initial_sample_size...
%     bundle_adj_size repr_err_th

%% Means for all triplets

means_all=zeros(8,5,2);
for m=methods_to_test
    means_all(m,:,1)=[mean(repr_err(:,m,1)), mean(rot_err(:,m,1)),...
        mean(t_err(:,m,1)), mean(iter(:,m,1)), mean(time(:,m,1))];
    means_all(m,:,2)=[mean(repr_err(:,m,2)), mean(rot_err(:,m,2)),...
        mean(t_err(:,m,2)), mean(iter(:,m,2)), mean(time(:,m,2))];
end



%% latex table

mod=[1:5,7:8];

input.data=zeros(8,5);
input.data=[[means_all(mod,[1:3,5],1),means_all(mod,4,2)];...
    [means_all(1,1:3,2),NaN,NaN]];

input.tableColLabels = {'repr. err. (px)', '\textsf{R} err.','\textsf{t} err.','init. time (s)','iter. BA'};
input.tableRowLabels = {'\textsf{TFT-L}','\textsf{TFT-R}','\textsf{TFT-N}',...
    '\textsf{TFT-PF}','\textsf{TFT-PH}', '\textsf{F-L}', '\textsf{F-O}', '\textsf{BA}'};
input.dataFormat = {'%.3f',4,'%.2f',1}; 
input.dataNanString = ' ';


latex = latexTable(input);




%% plot scene

switch dataset
    case 'fountain-P11'
        num_cam=11;
    case 'Herz-Jesu-P8'
        num_cam=8;
end
    
cam_centers=zeros(3,num_cam);
for im=1:num_cam
    [K,R,t]=readCalibrationOrientation_Strecha(path_to_data,im_names{im});
    cam_centers(:,im)=-R'*t;
end
figure;
scatter3(cam_centers(1,:),cam_centers(2,:),cam_centers(3,:));

%% compute angles
angles=zeros(length(triplets_to_test),1);
for it=1:length(triplets_to_test)
    
    triplet=indexes_sorted(triplets_to_test(it) ,1:3); 
    im1=triplet(1); im2=triplet(2);  im3=triplet(3);
    
    % Calibration info
    [K1,R1,t1]=readCalibrationOrientation_Strecha(path_to_data,im_names{im1});
    [K2,R2,t2]=readCalibrationOrientation_Strecha(path_to_data,im_names{im2});
    [K3,R3,t3]=readCalibrationOrientation_Strecha(path_to_data,im_names{im3});
    C1=-R1'*t1; C2=-R2'*t2; C3=-R3'*t3;
    a=C2-C1; a=a/norm(a); 
    b=C3-C1; b=b/norm(b);
    c=C3-C2; c=c/norm(c);
    
    angle1=acos(a'*b)*180/pi;
    angle2=acos(-a'*c)*180/pi;
    angle3=acos(b'*c)*180/pi;
    angles(it)=max([angle1,angle2,angle3]);
    %fprintf('Triplet %d/%d (%d,%d,%d) with angle %f.\n',...
    %    it,length(triplets_to_evaluate),im1,im2,im3,angle);
end

figure;
N=histc(angles,90:5:180);
bar(90:5:180,N,'histc')

%%

figure;
scatter(angles,repr_err(:,1,1)-repr_err(:,3,1));
figure;
scatter(angles,repr_err(:,1,1)-repr_err(:,5,1));


figure;
scatter(angles,t_err(:,1,1)-t_err(:,3,1));
figure;
scatter(angles,t_err(:,1,1)-t_err(:,5,1));




