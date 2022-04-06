function [log_lik,gradient_mean,hessian_mean] = complete_sum_likelihood(times, v,alpha,beta, end_time, old_weights, p, N)  %output: log_likelihoods
% N: the number of monte carlo samples


%This is the Q function of the EM algorithm you need to numerically
%maximise at each iteration.
gradient_mean = -1e10*ones(size([v,alpha,beta]));
hessian_mean = -1e10*ones(numel([v,alpha,beta]),numel([v,alpha,beta]));

gradient_mean = 1e10*ones(size([v,alpha,beta]));
hessian_mean = 1e10*ones(numel([v,alpha,beta]),numel([v,alpha,beta]));
if any(v) <= 0
    log_lik = inf;
    %comp_intensities = 0;
    %disp('v: negative baseline')
    return
end
if min(min(alpha)) < 0 
    log_lik = inf;
    %comp_intensities = 0;
    %disp('A: negative excitation matrix')
    return
end
if any(any(alpha>=beta)) %max(max(A./B)) >= 1 
    log_lik = inf;
    %comp_intensities = 0;
    %disp('A greater than B')
    return
end
if det(eye(p) - alpha./beta) <= 0 %+ 0.01 %|| max(eig(A./B)) >= 1
    log_lik = inf;
    %comp_intensities = 0;
    %disp('1) A./B: determinant less than zero')
    return
end
%Spectral radius of Gamma = \int_0^inf G(u)du is < 1.
if any((eye(p) - alpha./beta)\v < 0)  %any((eye(p) - A./B)\v < 0) %det(eye(p) - A./B) <= 0  %(1-A./B)\v < 0 
    log_lik = inf;
    %comp_intensities = 0;
    %disp('2) A./B: determinant less than zero')
    return
end

%N = size(times,1); % --- make sure this is getting the correct dimension
log_likelihoods = 1:N; 
gradients = zeros(p+2*p^2,N);
hessian = zeros(p+2*p^2,p+2*p^2,N);
if iscell(times)
%    N = size(times,2); % --- make sure this is getting the correct dimension
%tic;
    for i=1:N %can this be vectorised? 
        %log_likelihoods(i) = complete_likelihood(times{i}, v,alpha,beta, end_time,p)*old_weights(i); 
        [log_likelihoods(i),gradients(:,i),hessian(:,:,i)] = complete_likelihood(times{i}, v,alpha,beta, end_time,p); 
        log_likelihoods(i) = log_likelihoods(i)*old_weights(i);
        gradients(:,i) = gradients(:,i)*old_weights(i);
        hessian(:,:,i) = hessian(:,:,i)*old_weights(i);
    end     
%toc
else
    % vectorised version
    [log_likelihoods, gradients] = complete_likelihood(times,v,alpha,beta,end_time,p);  
    if size(log_likelihoods) ~= size(old_weights) % check dimensions first
        old_weights = old_weights';
    end
    log_likelihoods = log_likelihoods.*old_weights;
    if size(gradients) ~= size(old_weights) % check dimensions first
        old_weights = old_weights';
    end
    gradients = gradients.*old_weights;
end
log_lik = mean(log_likelihoods); 
gradient_mean = mean(gradients,2);
hessian_mean = mean(hessian,3);
% disp(gradient_mean)
end 
