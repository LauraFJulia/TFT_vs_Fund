function [inliers,model,ransac_th]=AC_RANSAC(Data,func_model,N,n_sample,d,a,n_out,NFA_th,max_it)
% AC-RANSAC: A contrario random sample consensus. Based on the paper of P.
% Moulon, P. Monasse, R. Marlet: Adaptive Structure from Motion with a
% contrario model estimation.
%
% Input arguments:
% Data       - array/cell of correspondences or image points and their 3D
%              position, or any other necessary data for the model to be
%              computed.
% func_model - function that gives the model (fundamental matrix,
%              orientation, TFT...) with a minimal number of data and also
%              the point-wise rigidness/error of the rest of the data w.r.t. the model.
%              it is of the form [model,vector_errors]=func_model(Data,sample)
% N          - size of the data.
% n_sample   - minimal cardinal of necessary data to compute the model.
% d          - dimention of the error.
% a          - probability of a random data point of having error 1 pixel.
% n_out      - number of possible models that can arise from n_sample data.
% NFA_th     - threshold for NFA that establishes the validity of a
%              method (1 by default).
% max_it     - maximum number of iterations for ransac (100 by default).
%
% Output arguments:
% inliers    - indexes indicating the final inliers in Data.
% model      - final model giving such set of inliers.
% ransac_th  - threshold used for ransac given by the AC method and
%              used in the final model&inliers computation.

% default parameters max_it and NFA_th
if nargin<9 || isempty(max_it)
    max_it=100;
end
if nargin<8 || isempty(NFA_th)
    NFA_th=1;
end

k_inliers=n_sample+1; % max number of inliers found
inliers=[];           % set of inliers
ransac_th=inf;        % ransac threshold

it=0; max_old=max_it;
while it<max_it
    it=it+1;
    
    % pick random sample of n_sample data
    rng(it);
    sample=randsample(N,n_sample);
    
    % compute model and errors from this sample
    % (if the function fails to give a model, we start the next iteration)
    [model_it,vec_errors]=func_model(Data,sample);
    if isempty(model_it)
        if max_it<2*max_old
            max_it=max_it+1;
        end
        continue;
    end
    % data not in the sample used
    nosample=setdiff(1:N,sample);
    
    % sort the list of errors
    [~,ind_sorted]=sort(vec_errors(nosample));
    
    % search for minimum of NFA(model,k)
    NFA_min=NFA_th; k_min=0; err_threshold=inf;
    factor=n_out*prod(N-n_sample:N)/factorial(n_sample);
    for k=n_sample+1:N
        factor=factor*( (N-(k-1))/(k-n_sample) )*a;
        NFA=factor*(vec_errors(nosample(ind_sorted(k-n_sample))))^(d*(k-n_sample));
        if NFA<=NFA_min
            NFA_min=NFA;
            k_min=k;
            err_threshold=vec_errors(nosample(ind_sorted(k-n_sample)));
        end
    end
    
    % If the found model has more inliers or the same number with less
    %  error than the previous we keep it
    if k_min>k_inliers || (k_min==k_inliers && err_threshold<ransac_th)  
        k_inliers=k_min;
        inliers=[sample.' nosample(ind_sorted(1:k_inliers-n_sample))];
        ransac_th=err_threshold;
        model=model_it;
    end 
end

    
    
end





