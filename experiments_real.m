% Script to recreate experiments on real data for the PSIVT paper
%  submission "A Critical Review of the Trifocal Tensor Estimation"

clear; close all;


%% Here uncomment the dataset to use.
dataset='fountain-P11';  
% dataset='Herz-Jesu-P8';


%% Some parameters
path_to_data=strcat('Data/',dataset,'/');

switch dataset
    case 'fountain-P11'
        triplets_to_evaluate=1:70;
    case 'Herz-Jesu-P8'
        triplets_to_evaluate=1:50;
end

initial_sample_size=100;
bundle_adj_size=50;

repr_err_th=20;

%% Recover correspondances

corresp_file=matfile(strcat(path_to_data,'Corresp_triplets','.mat'));
indexes_sorted=corresp_file.indexes_sorted;
corresp_by_triplet=corresp_file.Corresp;
im_names=corresp_file.im_names;
clear corresp_file;

%% methods to evaluate

method={...
    @LinearTFTPoseEstimation,...    % 1 - Linear TFT
    @MinimalTFTPoseEstimation,...   % 2 - Minimal TFT (Ressl)
    @NordbergTFTPoseEstimation,...  % 3 - Minimal TFT (Nordberg)
    @PapaFaugTFTPoseEstimation,...  % 4 - Papadopoulo Faugeras TFT
    @PiPoseEstimation,...           % 5 - Pi matrices (Ponce&Hebert)
    @PiColPoseEstimation,...        % 6 - Pi matrices - collinear (Ponce&Hebert)
    @LinearFPoseEstimation,...      % 7 - Linear Fundamental matrices
    @OptimFPoseEstimation};         % 8 - Fundamental matrices

methods_to_test=[1:5,7:8];

%% error vectors
repr_err=zeros(length(triplets_to_evaluate),length(method),2);
rot_err=zeros(length(triplets_to_evaluate),length(method),2);
t_err=zeros(length(triplets_to_evaluate),length(method),2);
iter=zeros(length(triplets_to_evaluate),length(method),2);
time=zeros(length(triplets_to_evaluate),length(method),2);

%% evaluation

for it=1:length(triplets_to_evaluate)
    
    triplet=indexes_sorted(triplets_to_evaluate(it) ,1:3); 
    im1=triplet(1); im2=triplet(2);  im3=triplet(3);
    Corresp=corresp_by_triplet{im1,im2,im3}';
    N=size(Corresp,2);
    
    fprintf('Triplet %d/%d (%d,%d,%d) with %d matching points.\n',...
        it,length(triplets_to_evaluate),im1,im2,im3,N);
    
    % Calibration info
    [K1,R1_true,t1_true,im_size1]=readCalibrationOrientation_Strecha(path_to_data,im_names{im1});
    [K2,R2_true,t2_true]=readCalibrationOrientation_Strecha(path_to_data,im_names{im2});
    [K3,R3_true,t3_true]=readCalibrationOrientation_Strecha(path_to_data,im_names{im3});
    CalM=[K1;K2;K3]; im_size=im_size1;
    R_t0={[R2_true*R1_true.', t2_true-R2_true*R1_true.'*t1_true],...
        [R3_true*R1_true.' t3_true-R3_true*R1_true.'*t1_true]};
    
    %%% Discart correspondances with repr_err > 1 pix
    Reconst0=triangulation3D({K1*eye(3,4),K2*R_t0{1},K3*R_t0{2}},Corresp);
    Reconst0=bsxfun(@rdivide,Reconst0(1:3,:),Reconst0(4,:));
    Corresp_new=project3Dpoints(Reconst0,CalM,[eye(3,4);R_t0{1};R_t0{2}]);
    residuals=Corresp_new-Corresp;
    mask=sum(abs(residuals)>repr_err_th,1)==0;
    Corresp=Corresp(:,mask);
    N=size(Corresp,2);
    REr=ReprError({K1*eye(3,4),K2*R_t0{1},K3*R_t0{2}},Corresp);
    fprintf('%d valid correspondances with reprojection error %f.\n',N,REr);
    
%     % RANSAC with fundamental matrices
%     %a=2*sqrt((2*K1(1,3))^2+(2*K1(2,3))^2)/(4*K1(1,3)*K1(2,3));
%     a=4*norm(im_size)^2/((im_size(1)*im_size(2))^2);
%     inliers21=AC_RANSAC({Corresp(1:2,:),Corresp(3:4,:)},@fundamental_model_ransac,N,8,1,a,1,1,40);
%     inliers31=AC_RANSAC({Corresp(1:2,:),Corresp(5:6,:)},@fundamental_model_ransac,N,8,1,a,1,1,40);
%     inliers32=AC_RANSAC({Corresp(3:4,:),Corresp(5:6,:)},@fundamental_model_ransac,N,8,1,a,1,1,40);
%     inliers=intersect(intersect(inliers21,inliers31),inliers32);
%     fprintf('Ransac found %d inliers out of %d.\n',size(inliers,2),N);
    inliers=1:N; % test without ransac

    % samples for initial estimation and bundle adjustment
    Corresp_inliers=Corresp(:,inliers);
    rng(it);
    init_sample=randsample(inliers,min(initial_sample_size,length(inliers)));
    rng(it);
    ref_sample=randsample(init_sample,min(bundle_adj_size,size(init_sample,2)));
    Corresp_init=Corresp(:,init_sample);
    Corresp_ref=Corresp(:,ref_sample);
    
    fprintf('method... ');
    for m=methods_to_test
        fprintf('%d ',m);
        if (m>6 && N<8) || N<7 % if not enough matches
            repr_err(it,m,:)=inf;    rot_err(it,m,:)=inf;
            t_err(it,m,:)=inf;       iter(it,m,:)=inf;
            time(it,m,:)=inf;
            continue;
        end
        
        % % Pose estimation by method m, measuring time            
        t0=cputime;
        [R_t_2,R_t_3,~,~,nit]=method{m}(Corresp_init,CalM);
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
        
        
        fprintf('(ref)... ');
        % % Apply Bundle Adjustment
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

%save Results/errors_1_to_70_triplets_fountain_100i_50r_40_reprinliers.mat repr_err rot_err t_err iter time...
%    methods_to_test method dataset triplets_to_evaluate initial_sample_size...
%    bundle_adj_size repr_err_th

save Results/errors_1_to_50_triplets_herzjesu_100i_50r_20_reprinliers.mat repr_err rot_err t_err iter time...
    methods_to_test method dataset triplets_to_evaluate initial_sample_size...
    bundle_adj_size repr_err_th

%% Means for all triplets

means_all=zeros(8,5,2);
for m=methods_to_test
    means_all(m,:,1)=[mean(repr_err(:,m,1)), mean(rot_err(:,m,1)),...
        mean(t_err(:,m,1)), mean(iter(:,m,1)), mean(time(:,m,1))];
    means_all(m,:,2)=[mean(repr_err(:,m,2)), mean(rot_err(:,m,2)),...
        mean(t_err(:,m,2)), mean(iter(:,m,2)), mean(time(:,m,2))];
end



%% Means only for valid triplets

rot_err_thresh=5;
trans_err_thresh=10;


indexes_valid=zeros(8,length(triplets_to_evaluate));
num_valid_cases=zeros(8,1);
for m=methods_to_test
    aux=[rot_err(:,m,2)>rot_err_thresh, t_err(:,m,2)>trans_err_thresh];
    indexes_valid(m,:)=(sum(aux,2)==0);
    num_valid_cases(m)=100*sum(indexes_valid(m,:))/length(triplets_to_evaluate);
end

% common valid triplets
indexes_valid_all=sum(~indexes_valid(methods_to_test,:),1)==0;

means_valid=zeros(8,5,2);
for m=methods_to_test
    means_valid(m,:,1)=[mean(repr_err(indexes_valid_all,m,1)),...
        mean(rot_err(indexes_valid_all,m,1)),...
        mean(t_err(indexes_valid_all,m,1)),...
        mean(iter(indexes_valid_all,m,1)),...
        mean(time(indexes_valid_all,m,1))];
    means_valid(m,:,2)=[mean(repr_err(indexes_valid_all,m,2)),...
        mean(rot_err(indexes_valid_all,m,2)),...
        mean(t_err(indexes_valid_all,m,2)),...
        mean(iter(indexes_valid_all,m,2)),...
        mean(time(indexes_valid_all,m,2))];
end

%% latex table

mod=[1:5,7:8];

input.data=zeros(8,5);
input.data=[[means_all(mod,[1:3,5],1),means_all(mod,4,2)];...
    [means_all(1,1:3,2),NaN,NaN]];

input.tableColLabels = {'repr. err. (px)', '\textsf{R} err.','\textsf{t} err.','init. time (s)','iter. BA'};
input.tableRowLabels = {'\textsf{TFT-L}','\textsf{TFT-R}','\textsf{TFT-N}',...
    '\textsf{TFT-PF}','\textsf{TFT-PH}', '\textsf{F-L}', '\textsf{F-O}', '\textsf{BA}'};
input.dataFormat = {'%.3f',3,'%.2f',2}; 
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

% compute angles
angles=zeros(length(triplets_to_evaluate),1);
for it=1:length(triplets_to_evaluate)
    
    triplet=indexes_sorted(triplets_to_evaluate(it) ,1:3); 
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




